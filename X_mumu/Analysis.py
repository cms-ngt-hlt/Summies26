import ROOT
import uproot
import glob
import os
import sys
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import argparse
import shutil, os
from scipy.optimize import curve_fit

ROOT.gROOT.SetBatch(True)
hep.style.use(hep.style.CMS)

from pathlib import Path
import tempfile

REPO_DIR = Path(__file__).resolve().parent

SRC = REPO_DIR / "functions.cc"
HDR_DIR = REPO_DIR / "headers"

tmp_dir = Path(tempfile.mkdtemp(prefix="root_aclic_"))
TMP = tmp_dir / "functions.cc"

shutil.copy(SRC, TMP)
ROOT.gROOT.ProcessLine(f'.L {TMP}+')

#─────────────────────────────────────── Helper functions ────────────────────────────────────────

def flatten(arr):       #flatten arrays, needed because most of the variables are vectors of size 2
    if len(arr) == 0:
        return np.array([])
    first = arr[0]
    if hasattr(first, "__len__"):  
        return np.concatenate([np.asarray(ev) for ev in arr])
    return np.asarray(arr) 

def files_with_events(file_names):
    good_files = []
    for file_name in file_names:
        root_file = ROOT.TFile.Open(file_name)
        if root_file and not root_file.IsZombie() and root_file.Get("Events"):
            good_files.append(file_name)
        if root_file:
            root_file.Close()
    return good_files

def make_chain(file_names):
    chain = ROOT.TChain("Events")
    for file_name in file_names:
        chain.Add(file_name)
    return chain

def gaussian(x, A, mu, sigma):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

def fit_gaussian(counts, edges):
    centers = 0.5 * (edges[:-1] + edges[1:])
    
    A0 = counts.max()
    mu0 = centers[np.argmax(counts)]
    sigma0 = (centers[counts > 0.5 * np.argmax(counts)][-1] - centers[counts > 0.5 * np.argmax(counts)][0])/2
    
    print(A0, mu0, sigma0)

    counts_fit = counts[(centers >= mu0 - 4*sigma0) & (centers <= mu0 + 4*sigma0)]
    centers_fit = centers[(centers >= mu0 - 4*sigma0) & (centers <= mu0 + 4*sigma0)]

    counts_fit_err = [np.sqrt(c) if c != 0 else 1 for c in counts_fit]
    popt, pcov = curve_fit(gaussian, centers_fit, counts_fit, p0=[A0, mu0, sigma0], sigma=counts_fit_err)
    perr = np.sqrt(np.diag(pcov))
    
    A_fit, mu_fit, sigma_fit = popt
    var_mu = pcov[1, 1]
    var_sigma = pcov[2, 2]
    cov_mu_sigma = pcov[1, 2]

    print(popt, perr)
    
    sigma_fit = np.abs(sigma_fit)
    resolution = sigma_fit / mu_fit
    
    dR_dsigma = 1 / mu_fit
    dR_dmu = -sigma_fit / mu_fit**2
    
    var_resolution = (dR_dsigma**2) * var_sigma + (dR_dmu**2) * var_mu + 2 * dR_dsigma * dR_dmu * cov_mu_sigma
    resolution_err = np.sqrt(var_resolution)
    
    return popt, perr, resolution, resolution_err
#─────────────────────────────────────── Define variables ────────────────────────────────────────

def apply_defines(df):
    return (
        df
        .Define("Gen_muon_idx",             "get_muon_idx(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags)")
        .Define("nGen_muons",               "Gen_muon_idx.size()")
        .Define("Gen_muon_all_dR",          "get_dR(Gen_muon_idx, GenPart_eta, GenPart_phi)")

        .Define("Gen_muon_pair_indices",    "deltaR_muon_pairing(Gen_muon_idx, GenPart_eta, GenPart_phi, GenPart_genPartIdxMother)")
        .Define("nGen_muons_pairs",         "Gen_muon_idx.size()")

        .Define("Gen_muon_pt",              "get_pair_variable(Gen_muon_idx, GenPart_pt)")
        .Define("Gen_muon_eta",             "get_pair_variable(Gen_muon_idx, GenPart_eta)")
        .Define("Gen_muon_phi",             "get_pair_variable(Gen_muon_idx, GenPart_phi)")
        .Define("Gen_muon_mass",            "get_pair_variable(Gen_muon_idx, GenPart_mass)")
        .Define("Gen_muon_dR",              "get_dR(Gen_muon_idx, GenPart_eta, GenPart_phi)")

        .Define("Gen_X_pt",                 "get_pair_pt(Gen_muon_idx, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
        .Define("Gen_X_eta",                "get_pair_eta(Gen_muon_idx, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
        .Define("Gen_X_phi",                "get_pair_phi(Gen_muon_idx, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
        .Define("Gen_X_M",                  "get_pair_M(Gen_muon_idx, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")

### L1
        .Define("L1_muon",                  "truth_matching(Gen_muon_idx, GenPart_eta, GenPart_phi, L1GTgmtTkMuon_pt, L1GTgmtTkMuon_eta, L1GTgmtTkMuon_phi)")
        
        .Define("nL1_muon",                 "L1_muon.size()")
        .Define("L1_muon_pt",               "get_reco_variable(L1_muon, L1GTgmtTkMuon_pt)")
        .Define("L1_muon_eta",              "get_reco_variable(L1_muon, L1GTgmtTkMuon_eta)")
        .Define("L1_muon_phi",              "get_reco_variable(L1_muon, L1GTgmtTkMuon_phi)")
        .Define("L1_muon_dR",               "get_dR(L1_muon, L1GTgmtTkMuon_eta, L1GTgmtTkMuon_phi)")

        .Define("L1_X_pt",                  "get_pair_pt(L1_muon, L1GTgmtTkMuon_pt, L1GTgmtTkMuon_eta, L1GTgmtTkMuon_phi, 0.1057f)")
        .Define("L1_X_eta",                 "get_pair_eta(L1_muon, L1GTgmtTkMuon_pt, L1GTgmtTkMuon_eta, L1GTgmtTkMuon_phi, 0.1057f)")
        .Define("L1_X_phi",                 "get_pair_phi(L1_muon, L1GTgmtTkMuon_pt, L1GTgmtTkMuon_eta, L1GTgmtTkMuon_phi, 0.1057f)")
        .Define("L1_X_M",                   "get_pair_M(L1_muon, L1GTgmtTkMuon_pt, L1GTgmtTkMuon_eta, L1GTgmtTkMuon_phi, 0.1057f)")
### HLT
        .Define("Reco_muon",                "truth_matching(Gen_muon_idx, GenPart_eta, GenPart_phi, hltMuon_pt, hltMuon_eta, hltMuon_phi)")
        
        .Define("nReco_muon",               "Reco_muon.size()")
        .Define("Reco_muon_pt",             "get_reco_variable(Reco_muon, hltMuon_pt)")
        .Define("Reco_muon_eta",            "get_reco_variable(Reco_muon, hltMuon_eta)")
        .Define("Reco_muon_phi",            "get_reco_variable(Reco_muon, hltMuon_phi)")
        .Define("Reco_muon_dR",             "get_dR(Reco_muon, hltMuon_eta, hltMuon_phi)")

        .Define("Reco_X_pt",                "get_pair_pt(Reco_muon, hltMuon_pt, hltMuon_eta, hltMuon_phi, 0.1057f)")
        .Define("Reco_X_eta",               "get_pair_eta(Reco_muon, hltMuon_pt, hltMuon_eta, hltMuon_phi, 0.1057f)")
        .Define("Reco_X_phi",               "get_pair_phi(Reco_muon, hltMuon_pt, hltMuon_eta, hltMuon_phi, 0.1057f)")
        .Define("Reco_X_M",                 "get_pair_M(Reco_muon, hltMuon_pt, hltMuon_eta, hltMuon_phi, 0.1057f)")
    )

def get_dataframe(base_dir): #Creates a dataframe for each root file individually
    files = sorted(glob.glob(base_dir + "/step2_*.root"))
    files = files_with_events(files)
    chain = make_chain(files)
    df = ROOT.RDataFrame(chain)
    print(f"total entries for mX = {mass} GeV: {df.Count().GetValue()}")
    return apply_defines(df)

def overall_efficiency(all, trigger):
    all_c = all.Count().GetValue()
    trigger_c = trigger.Count().GetValue()
    eff = trigger_c/all_c
    return eff

#─────────────────────────────────────── Plotting ────────────────────────────────────────

def plot(data_per_mass, var, bins, hist_range, xlabel, save_as, scale):
    path = f"{plot_dir}/Displaced/"
    os.makedirs(path, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 8))
    for mass, arrs in data_per_mass.items():
        counts, edges = np.histogram(arrs[var], bins=bins, range=hist_range, density=False)
        hep.histplot(counts, edges, ax=ax, label=f"m$_X$ = {mass} GeV", color=mass_colors[mass], histtype="step", linewidth=2)
    ax.set_xlabel(xlabel, fontsize=24)
    ax.set_ylabel("Events", fontsize=24)
    ax.legend(loc="upper right", fontsize=20)
    hep.cms.label("Preliminary", data=False, ax=ax, com=14)
    plt.tight_layout()

    if save_as == 'png':            
        if scale == 'log':
            plt.yscale('log')
            fig.savefig(f"{path}{var}_log.png")
        fig.savefig(f"{path}{var}.png")
        plt.close(fig)
        print(f"Saved {var}.png")
    else:
        if scale == 'log':
            plt.yscale('log')
            fig.savefig(f"{path}{var}_log.pdf")
        fig.savefig(f"{path}{var}.pdf")
        plt.close(fig)
        print(f"Saved {var}.pdf")


def plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match, var, bins, hist_range, xlabel):
    path = f"{plot_dir}/Comparisons/"
    os.makedirs(path, exist_ok=True)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    for mass in data_per_mass:
        arrs     = data_per_mass[mass]
        arrs_NGT = data_per_mass_NGT[mass]
        arrs_match = data_per_mass_match[mass]

        counts, edges = np.histogram(arrs[var], bins=bins, range=hist_range, density=False)
        hep.histplot(counts, edges, ax=ax, label=f"m$_X$ = {mass} GeV (All)",
                     color=mass_colors[mass], histtype="step", hatch= "", linewidth=2, linestyle="-")

        counts_NGT, edges_NGT = np.histogram(arrs_NGT[var], bins=bins, range=hist_range, density=False)
        hep.histplot(counts_NGT, edges_NGT, ax=ax, label=f"m$_X$ = {mass} GeV (NGT)",
                     color=mass_colors[mass], histtype="step", hatch= "/", linewidth=2, linestyle="--")

        
        counts_match, edges_match = np.histogram(arrs_match[var], bins=bins, range=hist_range, density=False)
        hep.histplot(counts_match, edges_match, ax=ax, label=f"m$_X$ = {mass} GeV (NGT)",
                     color=mass_colors[mass], histtype="fill", hatch= "", linewidth=2, linestyle="-")

    ax.set_xlabel(xlabel, fontsize=24)
    ax.set_ylabel("Events", fontsize=24)
    ax.legend(loc="upper right", fontsize=14)
    hep.cms.label("Preliminary", data=False, ax=ax, com=14)
    plt.tight_layout()
    fig.savefig(f"{path}{var}.png")
    plt.close(fig)
    print(f"Saved Comparison_{var}.png")

def plot_efficiency(mass, data_per_mass, data_per_mass_NGT, var, bins, hist_range, xlabel, eff_color):
    path = f"{plot_dir}/Efficiencies/"
    os.makedirs(path, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 8))

    arrs     = data_per_mass[mass]
    arrs_NGT = data_per_mass_NGT[mass]

    counts, edges = np.histogram(arrs[var], bins=bins, range=hist_range, density=False)
    hep.histplot(counts, edges, ax=ax, label=f"m$_X$ = {mass} GeV (All)", color=mass_colors[mass], histtype="step", linewidth=2)
    
    counts_NGT, edges_NGT = np.histogram(arrs_NGT[var], bins=bins, range=hist_range, density=False)
    hep.histplot(counts_NGT, edges_NGT, ax=ax, label=f"m$_X$ = {mass} GeV (NGT)", color=mass_colors[mass], histtype="fill", linewidth=2)
    
    Efficiency_NGT = [num/den if den != 0 else 0 for num,den in zip(counts_NGT, counts)]
    
    centers = 0.5*(edges[:-1]+edges[1:])
    bin_widths = np.diff(edges)
    y_error_NGT = [np.sqrt(eff * (1 - eff) / den) if den > 0 else 0 for eff, den in zip(Efficiency_NGT, counts)] 
    
    ax2 = ax.twinx()
    ax2.errorbar(centers, Efficiency_NGT, xerr=0.5*bin_widths, yerr = y_error_NGT, color=eff_color, fmt = '*')
    ax2.set_ylabel("Efficiency", color=eff_color)
    ax2.tick_params(axis='y', labelcolor=eff_color)
    ax2.set_ylim(0,1)
    
    ax.set_xlabel(xlabel, fontsize=24)
    ax.set_ylabel("Events", fontsize=24)
    ax.legend(loc="upper right", fontsize=24)
    hep.cms.label("Preliminary", data=False, ax=ax, com=14)
    plt.tight_layout()
    fig.savefig(f"{path}{mass}_{var}.png")
    plt.close(fig)
    print(f"Saved Efficiency_{mass}_{var}.png")


def plot_genVsReco(mass, data_per_mass_match, var_gen, var_reco, bins, xlabel):
    path = f"{plot_dir}/GenVsReco/"
    os.makedirs(path, exist_ok=True)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    arrs_match = data_per_mass_match[mass]
        
    counts_gen, edges_gen = np.histogram(arrs_match[var_gen], bins=bins, range=(mass-0.5, mass+0.5), density=False)
    hep.histplot(counts_gen, edges_gen, ax=ax, label=f"GEN",
                    color=mass_colors[mass], histtype="step", hatch= "", linewidth=2, linestyle="-")
    
    counts_reco, edges_reco = np.histogram(arrs_match[var_reco], bins=bins, range=(mass-0.5, mass+0.5), density=False)
    hep.histplot(counts_reco, edges_reco, ax=ax, label=f"RECO",
                    color=mass_colors[6], histtype="step", hatch= "", linewidth=2, linestyle="-")

    ax.set_xlabel(xlabel, fontsize=24)
    ax.set_ylabel("Events", fontsize=24)
    ax.legend(loc="upper right", fontsize=16)
    hep.cms.label("Preliminary", data=False, ax=ax, com=14)
    plt.tight_layout()
    fig.savefig(f"{path}{mass}_GenVsReco{var_gen}.png")
    plt.close(fig)
    print(f"Saved genVsReco_{var_gen}.png")

#─────────────────────────────────────── Fitting ────────────────────────────────────────
def plot_fit(mass, data_per_mass_match, var_gen, bins, xlabel, var_reco=None):
    path = f"{plot_dir}/Fit/"
    os.makedirs(path, exist_ok=True)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    arrs_match = data_per_mass_match[mass]

    #GEN
    counts, edges = np.histogram(arrs_match[var_gen], bins=bins, range=(mass-0.5, mass+0.5), density=False)
    hep.histplot(counts, edges, ax=ax, label="GEN",
                 color=mass_colors[mass], histtype="step", hatch="", linewidth=2, linestyle="-")

    popt, perr, resolution, resolution_err = fit_gaussian(counts, edges)
    A_fit, mu_fit, sigma_fit = popt
    mu_err, sigma_err = perr[1], perr[2]
    
    x_fit = np.linspace(edges[0], edges[-1], 500)
    y_fit = gaussian(x_fit, *popt)
    ax.plot(x_fit, y_fit, color="red", linestyle="--", linewidth=2, label=f"$\\mu$={mu_fit:.3f}±{mu_err:.3f}\n$\\sigma$={sigma_fit:.3f}±{sigma_err:.3f}")

    print(f"Mean GEN = {mu_fit}\nsigma = {sigma_fit}")

    #RECO   
    popt2, perr2 = None, None  #
    if var_reco is not None:
        counts2, edges2 = np.histogram(arrs_match[var_reco], bins=bins, range=(mass-0.5, mass+0.5), density=False)
        hep.histplot(counts2, edges2, ax=ax, label="RECO",
                     color=mass_colors[1], histtype="step", hatch="", linewidth=2, linestyle="-")

        popt2, perr2, resolution2, resolution_err2 = fit_gaussian(counts2, edges2)
        A_fit2, mu_fit2, sigma_fit2 = popt2
        mu_err2, sigma_err2 = perr2[1], perr2[2]

        x_fit2 = np.linspace(edges2[0], edges2[-1], 500)
        y_fit2 = gaussian(x_fit2, *popt2)
        ax.plot(x_fit2, y_fit2, color="blue", linestyle="--", linewidth=2,
                label=f"$\\mu$={mu_fit2:.3f}±{mu_err2:.3f}\n$\\sigma$={sigma_fit2:.3f}±{sigma_err2:.3f}")
                
        print(f"Mean RECO = {mu_fit2}\nsigma = {sigma_fit2}")

    ax.set_xlabel(xlabel, fontsize=24)
    ax.set_ylabel("Events", fontsize=24)
    ax.legend(loc="upper right", fontsize=16)
    hep.cms.label("Preliminary", data=False, ax=ax, com=14)
    plt.tight_layout()
    fig.savefig(f"{path}{mass}_fitGenVsReco.png")
    plt.close(fig)
    print(f"Saved {1}_fit_GenVsReco.png")


resolution_gen = {}
resolution_reco = {}
def plot_resolution(data_per_mass_match, var_gen, var_reco=None):
    path = f"{plot_dir}/Resolution/"
    os.makedirs(path, exist_ok=True)
    
    fig, ax = plt.subplots(figsize=(10, 8))

    for mass in data_per_mass:
        arrs_match = data_per_mass_match[mass]

        counts, edges = np.histogram(arrs_match[var_gen], bins=100, range=(mass-0.5, mass+0.5), density=False)
        popt, perr, resolution, resolution_err = fit_gaussian(counts, edges)
        A_fit, mu_fit, sigma_fit = popt

        print(f"Mean GEN = {mu_fit}\nsigma = {sigma_fit}")
        resolution_gen[mass] = resolution

        plt.errorbar(mass, resolution_gen[mass], yerr=resolution_err, color='forestgreen', marker="*", markersize="11", label="GEN")
        
        #RECO   
        popt2, perr2 = None, None  #
        if var_reco is not None:
            counts2, edges2 = np.histogram(arrs_match[var_reco], bins=100, range=(mass-0.5, mass+0.5), density=False)
            popt2, perr2, resolution2, resolution_err2 = fit_gaussian(counts2, edges2)
            A_fit2, mu_fit2, sigma_fit2 = popt2

            print(f"Mean RECO = {mu_fit2}\nsigma = {sigma_fit2}")
            resolution_reco[mass] = resolution2

            plt.errorbar(mass, resolution_reco[mass], yerr=resolution_err2, color='darkmagenta', marker="x", markersize="11", label="RECO")


    ax.set_xlabel(r"m$_X$ [GeV]", fontsize=24)
    ax.set_ylabel(r"Resolution $\sigma / \mu$", fontsize=24)
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))  # duplicate labels overwrite each other, keeping only one
    ax.legend(by_label.values(), by_label.keys(), loc="upper right", fontsize=16)
    hep.cms.label("Preliminary", data=False, ax=ax, com=14)
    plt.tight_layout()
    fig.savefig(f"{path}{var_gen}.png")
    plt.close(fig)
    print(f"Saved resolution_{var_gen}.png")

#─────────────────────────────────────── Main ────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Argument for Analysis_all')
    parser.add_argument('--run-comparison', action='store_true', required=False, help='Run comparison plots across all generated X masses.')
    parser.add_argument('--output-dir', required=False, help='Output directory for plots.', default='./')
    args = parser.parse_args()

    plot_dir = args.output_dir + '/plots/'
    os.makedirs(plot_dir, exist_ok=True)

    mass_dirs = {
        1:  "/eos/user/e/evernazz/Sarah/Xmumu/Xmumu1/",
        2:  "/eos/user/e/evernazz/Sarah/Xmumu/Xmumu2/",
        4:  "/eos/user/e/evernazz/Sarah/Xmumu/Xmumu4/",
        6:  "/eos/user/e/evernazz/Sarah/Xmumu/Xmumu6/",
        8:  "/eos/user/e/evernazz/Sarah/Xmumu/Xmumu8/",
        10: "/eos/user/e/evernazz/Sarah/Xmumu/Xmumu10/",
    }

    mass_colors = {
    1:  "mediumvioletred",
    2:  "purple",
    4:  "royalblue",
    6:  "seagreen",
    8:  "olivedrab",
    10: "yellowgreen",
    }

    to_plot = [ "nGen_muons", "Gen_muon_pt", "Gen_muon_eta", "Gen_muon_phi", "Gen_muon_mass",  "Gen_muon_dR",
                
                "Gen_X_pt", "Gen_X_eta", "Gen_X_phi", "Gen_X_M", 
                ###
                "nL1_muon", "L1_muon_pt", "L1_muon_eta", "L1_muon_phi", 

                "L1_X_pt", "L1_X_eta", "L1_X_phi", "L1_X_M", "L1_muon_dR",
                ###
                "nReco_muon","Reco_muon_pt", "Reco_muon_eta", "Reco_muon_phi", 

                "Reco_X_pt", "Reco_X_eta", "Reco_X_phi", "Reco_X_M", "Reco_muon_dR",
    ]

    data_per_mass = {}
    data_per_mass_NGT = {}
    data_per_mass_match = {}
    for mass, base_dir in mass_dirs.items():
        df = get_dataframe(base_dir)
        df_NGT = df.Filter("DST_PFScouting == true")
        df_match = df.Filter("Reco_muon[0] >= 0 && Reco_muon[1] >= 0")
        print(f"NGT Efficiency for mX = {mass} GeV: ", overall_efficiency(df, df_NGT))
        print(f"Matching Efficiency for mX = {mass} GeV: ", overall_efficiency(df, df_match))
        
        column = df.AsNumpy(to_plot)
        column_NGT = df_NGT.AsNumpy(to_plot)
        column_match = df_match.AsNumpy(to_plot)
        data_per_mass[mass] = {variable: flatten(column[variable]) for variable in to_plot}
        data_per_mass_NGT[mass] = {variable: flatten(column_NGT[variable]) for variable in to_plot}
        data_per_mass_match[mass] = {variable: flatten(column_match[variable]) for variable in to_plot}

    plot(data_per_mass, "nGen_muons",    bins=10, hist_range=(0, 10),          xlabel=r"Gen number of $\mu$", save_as='png', scale='linear')
    plot(data_per_mass, "Gen_muon_pt",   bins=30, hist_range=(0, 80),          xlabel=r"Gen p$_T$ $\mu$ [GeV]", save_as='png', scale='linear')
    plot(data_per_mass, "Gen_muon_eta",  bins=30, hist_range=(-4, 6),          xlabel=r"Gen $\eta(\mu)$", save_as='png', scale='linear')
    plot(data_per_mass, "Gen_muon_phi",  bins=30, hist_range=(-np.pi, np.pi),  xlabel=r"Gen $\phi(\mu)$", save_as='png', scale='linear')
    
    plot(data_per_mass, "Gen_muon_dR",bins=30, hist_range=(0, 1),           xlabel=r"Gen $\Delta R$ X [GeV]", save_as='png', scale='linear')

    plot(data_per_mass, "Gen_X_pt",   bins=30, hist_range=(0, 80),          xlabel=r"Gen p$_T$ X [GeV]", save_as='png', scale='linear')
    plot(data_per_mass, "Gen_X_eta",  bins=30, hist_range=(-4, 6),          xlabel=r"Gen $\eta(X)$", save_as='png', scale='linear')
    plot(data_per_mass, "Gen_X_phi",  bins=30, hist_range=(-np.pi, np.pi),  xlabel=r"Gen $\phi(X)$", save_as='png', scale='linear')
    plot(data_per_mass, "Gen_X_M",    bins=50, hist_range=(0, 12),          xlabel=r"Gen $m_X$ [GeV]", save_as='png', scale='linear')
    ##

    plot(data_per_mass, "nL1_muon", bins=10, hist_range=(0, 10),            xlabel=r"Reco number of $\mu$", save_as='png', scale='linear')

    plot(data_per_mass, "L1_muon_pt",  bins=30, hist_range=(0, 80),         xlabel=r"L1 p$_T$ $\mu$ [GeV]", save_as='png', scale='linear')
    plot(data_per_mass, "L1_muon_eta", bins=30, hist_range=(-3, 3),         xlabel=r"L1 $\eta(\mu)$", save_as='png', scale='linear')
    plot(data_per_mass, "L1_muon_phi", bins=30, hist_range=(-np.pi, np.pi), xlabel=r"L1 $\phi(\mu)$", save_as='png', scale='linear')

    plot(data_per_mass, "L1_X_pt",  bins=30, hist_range=(0, 80),         xlabel=r"L1 p$_T$ X [GeV]", save_as='png', scale='linear')
    plot(data_per_mass, "L1_X_eta", bins=30, hist_range=(-3, 3),         xlabel=r"L1 $\eta(X)$", save_as='png', scale='linear')
    plot(data_per_mass, "L1_X_phi", bins=30, hist_range=(-np.pi, np.pi), xlabel=r"L1 $\phi(X)$", save_as='png', scale='linear')
    plot(data_per_mass, "L1_X_M",   bins=30, hist_range=(0, 15),         xlabel=r"L1 $m_X$ [GeV]", save_as='png', scale='linear')

    ##
    plot(data_per_mass, "nReco_muon", bins=10, hist_range=(0, 10),          xlabel=r"Reco number of $\mu$", save_as='png', scale='linear')

    plot(data_per_mass, "Reco_muon_pt",  bins=30, hist_range=(0, 80),          xlabel=r"Reco p$_T$ $\mu$ [GeV]", save_as='png', scale='linear')
    plot(data_per_mass, "Reco_muon_eta", bins=30, hist_range=(-3, 3),          xlabel=r"Reco $\eta(\mu)$", save_as='png', scale='linear')
    plot(data_per_mass, "Reco_muon_phi", bins=30, hist_range=(-np.pi, np.pi),  xlabel=r"Reco $\phi(\mu)$", save_as='png', scale='linear')

    plot(data_per_mass, "Reco_X_pt",  bins=30, hist_range=(0, 80),          xlabel=r"Reco p$_T$ X [GeV]", save_as='png', scale='linear')
    plot(data_per_mass, "Reco_X_eta", bins=30, hist_range=(-3, 3),          xlabel=r"Reco $\eta(X)$", save_as='png', scale='linear')
    plot(data_per_mass, "Reco_X_phi", bins=30, hist_range=(-np.pi, np.pi),  xlabel=r"Reco $\phi(X)$", save_as='png', scale='linear')
    plot(data_per_mass, "Reco_X_M",   bins=30, hist_range=(0, 15),          xlabel=r"Reco $m_X$ [GeV]", save_as='png', scale='linear')

    plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match, "Gen_muon_dR", bins=30, hist_range=(0, 1),  xlabel=r"Gen $\Delta R$ X [GeV]")
    plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match, "Gen_X_pt",    bins=30, hist_range=(0, 80), xlabel=r"Gen p$_T$ X [GeV]")
    plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match, "Gen_X_eta",   bins=30, hist_range=(-4, 6), xlabel=r"Gen $\eta(X)$")
    plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match, "Gen_X_M",     bins=50, hist_range=(0, 12), xlabel=r"Gen $m_X$ [GeV]")
    plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match, "Gen_X_M",     bins=100, hist_range=(0, 15), xlabel=r"Gen $m_X$ [GeV]")

    plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match, "Reco_muon_dR", bins=30, hist_range=(0, 1),  xlabel=r"Reco $\Delta R$ X [GeV]")
    plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match, "Reco_X_pt",    bins=30, hist_range=(0, 80), xlabel=r"Reco p$_T$ X [GeV]")
    plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match, "Reco_X_eta",   bins=30, hist_range=(-4, 6), xlabel=r"Reco $\eta(X)$")
    plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match, "Reco_X_M",     bins=30, hist_range=(0, 15), xlabel=r"Reco $m_X$ [GeV]")

    plot_efficiency(6, data_per_mass, data_per_mass_NGT, "Gen_X_pt",    bins=30, hist_range=(0, 80),          xlabel=r"Gen p$_T$ X [GeV]",      eff_color="yellowgreen")
    plot_efficiency(6, data_per_mass, data_per_mass_NGT, "Gen_X_eta",   bins=30, hist_range=(-4, 6),          xlabel=r"Gen $\eta(X)$",          eff_color="yellowgreen")
    plot_efficiency(6, data_per_mass, data_per_mass_NGT, "Gen_X_phi",   bins=30, hist_range=(-np.pi, np.pi),  xlabel=r"Gen $\phi(X)$",          eff_color="yellowgreen")
    plot_efficiency(6, data_per_mass, data_per_mass_NGT, "Gen_X_M",     bins=30, hist_range=(2, 8),           xlabel=r"Gen $m_X$ [GeV]",        eff_color="yellowgreen")
    plot_efficiency(6, data_per_mass, data_per_mass_NGT, "Gen_muon_dR", bins=30, hist_range=(0, 4),           xlabel=r"Gen $\Delta R$ X [GeV]", eff_color="yellowgreen")

    plot_genVsReco(1, data_per_mass_match, "Gen_X_M", "Reco_X_M", bins=100, xlabel=r"$m_X$ [GeV]")

    plot_fit(1, data_per_mass_match, "Gen_X_M",  bins=100,  xlabel=r"$m_X$ [GeV]", var_reco="Reco_X_M")
    plot_fit(2, data_per_mass_match, "Gen_X_M",  bins=100,  xlabel=r"$m_X$ [GeV]", var_reco="Reco_X_M")
    plot_fit(4, data_per_mass_match, "Gen_X_M",  bins=100,  xlabel=r"$m_X$ [GeV]", var_reco="Reco_X_M")
    plot_fit(6, data_per_mass_match, "Gen_X_M",  bins=100,  xlabel=r"$m_X$ [GeV]", var_reco="Reco_X_M")
    plot_fit(8, data_per_mass_match, "Gen_X_M",  bins=100,  xlabel=r"$m_X$ [GeV]", var_reco="Reco_X_M")
    plot_fit(10, data_per_mass_match, "Gen_X_M",  bins=100,  xlabel=r"$m_X$ [GeV]", var_reco="Reco_X_M")

    plot_resolution(data_per_mass_match, "Gen_X_M", var_reco="Reco_X_M")
