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
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

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
process = r"$A' \to \mu\mu$"
m = r"m$_{A'}$"

color_all = "#FEA201"
color_L1 =  "#FE9BBA"
color_HLT = "#93D8D5"
color_NGT = "#FF3471"
eff_color = "#FEA201"

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
    
    counts_fit = counts[(centers >= mu0 - 4*sigma0) & (centers <= mu0 + 4*sigma0)]
    centers_fit = centers[(centers >= mu0 - 4*sigma0) & (centers <= mu0 + 4*sigma0)]

    counts_fit_err = [np.sqrt(c) if c != 0 else 1 for c in counts_fit]
    popt, pcov = curve_fit(gaussian, centers_fit, counts_fit, p0=[A0, mu0, sigma0], sigma=counts_fit_err)
    perr = np.sqrt(np.diag(pcov))
    
    A_fit, mu_fit, sigma_fit = popt
    var_mu = pcov[1, 1]
    var_sigma = pcov[2, 2]
    cov_mu_sigma = pcov[1, 2]
    
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
        .Define("Gen_muon_idx",             "get_muon_idx(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags, GenPart_pt)")
        .Define("nGen_muons",               "Gen_muon_idx.size()")
        .Define("Gen_muon_all_dR",          "get_dR(Gen_muon_idx, GenPart_eta, GenPart_phi)")

        .Define("Gen_muon_pair_indices",    "deltaR_muon_pairing(Gen_muon_idx, GenPart_eta, GenPart_phi, GenPart_genPartIdxMother)")
        .Define("nGen_muons_pairs",         "Gen_muon_idx.size()")

        .Define("Gen_muon_leading_pt",      "GenPart_pt[Gen_muon_idx[0]]")
        .Define("Gen_muon_leading_eta",     "GenPart_eta[Gen_muon_idx[0]]")
        .Define("Gen_muon_leading_phi",     "GenPart_phi[Gen_muon_idx[0]]")
        .Define("Gen_muon_leading_mass",    "GenPart_mass[Gen_muon_idx[0]]")

        .Define("Gen_muon_subleading_pt",   "GenPart_pt[Gen_muon_idx[1]]")
        .Define("Gen_muon_subleading_eta",  "GenPart_eta[Gen_muon_idx[1]]")
        .Define("Gen_muon_subleading_phi",  "GenPart_phi[Gen_muon_idx[1]]")
        .Define("Gen_muon_subleading_mass", "GenPart_mass[Gen_muon_idx[1]]")

        .Define("Gen_muon_dR",              "get_dR(Gen_muon_idx, GenPart_eta, GenPart_phi)")

        .Define("Gen_X_pt",                 "get_pair_pt(Gen_muon_idx, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
        .Define("Gen_X_eta",                "get_pair_eta(Gen_muon_idx, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
        .Define("Gen_X_phi",                "get_pair_phi(Gen_muon_idx, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
        .Define("Gen_X_M",                  "get_pair_M(Gen_muon_idx, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")

###-------------------- L1 objects

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

###------------------ HLT objects

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
    # print(f"total entries for mX = {mass} GeV: {df.Count().GetValue()}")
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

    max_count = 0
    for mass, arrs in data_per_mass.items():
        counts, edges = np.histogram(arrs[var], bins=bins, range=hist_range, density=False)
        hep.histplot(counts, edges, ax=ax, label=f"{m} = {mass} GeV", color=mass_colors[mass], histtype="fill", linewidth=2)
        max_count = max(max(counts), max_count)
        
    ax.set_xlabel(xlabel, fontsize=24)
    ax.set_ylabel("Entries", fontsize=24)
    ax.legend(loc="upper right", fontsize=18)

    hep.cms.label("Private work", loc=2, data=True,ax=ax, rlabel = f"{process} (0 PU) | 14 TeV", fontsize = 22)
    ax.grid(True, alpha=0.75, linestyle="dashdot", linewidth=0.75)
    ax.set_ylim(0, max_count * 1.3)

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


def plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match_L1, data_per_mass_match_HLT, var, bins, hist_range, xlabel):
    path = f"{plot_dir}/Comparisons/"
    os.makedirs(path, exist_ok=True)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    max_count = 0
    for mass in data_per_mass:
        arrs     = data_per_mass[mass]
        arrs_NGT = data_per_mass_NGT[mass]
        arrs_match_L1 = data_per_mass_match_L1[mass]
        arrs_match_HLT = data_per_mass_match_HLT[mass]

        counts, edges = np.histogram(arrs[var], bins=bins, range=hist_range, density=False)
        hep.histplot(counts, edges, ax=ax, label=f"m$_X$ = {mass} GeV", color=color_all, histtype="step", hatch="", linewidth=2)

        counts_NGT, edges_NGT = np.histogram(arrs_NGT[var], bins=bins, range=hist_range, density=False)
        hep.histplot(counts_NGT, edges_NGT, ax=ax, label="_nolegend_", color=color_NGT, histtype="fill", hatch=".", linewidth=2)

        counts_match, edges_match = np.histogram(arrs_match_L1[var], bins=bins, range=hist_range, density=False)
        hep.histplot(counts_match, edges_match, ax=ax, label="_nolegend_", color=color_L1, histtype="fill", hatch="//", linewidth=2)

        counts_match_HLT, edges_match_HLT = np.histogram(arrs_match_HLT[var], bins=bins, range=hist_range, density=False)
        hep.histplot(counts_match_HLT, edges_match_HLT, ax=ax, label="_nolegend_", color=color_HLT, histtype="fill", hatch="", linewidth=2)
        max_count = max(max(counts), max_count)

    ax.set_xlabel(xlabel, fontsize=24)
    ax.set_ylabel("Entries", fontsize=24)

    style_handles = [
        Line2D([0], [0], color=color_all, lw=2, linestyle="-",  label="All"),
        Patch(facecolor=color_NGT, edgecolor=color_NGT, label="NGT triggered"),
        Patch(facecolor=color_L1,  edgecolor=color_L1 , label=r"Matched L1 $\mu$'s"),
        Patch(facecolor=color_HLT, edgecolor=color_HLT, label=r"Matched HLT $\mu$'s"),
    ]
    ax.legend(handles=style_handles, loc="upper right", fontsize=18)

    hep.cms.label("Private work", loc=2, data=True,ax=ax, rlabel = f"{process} (0 PU) | 14 TeV", fontsize = 22)
    ax.grid(True, alpha=0.75, linestyle="dashdot", linewidth=0.75)
    ax.set_ylim(0, max_count * 1.3)

    fig.savefig(f"{path}{var}.png")
    plt.close(fig)
    print(f"Saved Comparison_{var}.png")

def plot_recoefficiency(mass, data_per_mass_NGT, data_per_mass_match, var, bins, hist_range, xlabel, eff_color):
    path = f"{plot_dir}/Reco_Efficiencies/"
    os.makedirs(path, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 8))

    arrs_NGT     = data_per_mass_NGT[mass]
    arrs_match   = data_per_mass_match[mass]

    counts_NGT, edges_NGT = np.histogram(arrs_NGT[var], bins=bins, range=hist_range, density=False)
    hep.histplot(counts_NGT, edges_NGT, ax=ax, label=f"m$_X$ = {mass} GeV (NGT)", color=mass_colors[mass], histtype="step", linewidth=2)
    
    counts_match, edges_match = np.histogram(arrs_match[var], bins=bins, range=hist_range, density=False)
    hep.histplot(counts_match, edges_match, ax=ax, label=f"m$_X$ = {mass} GeV (Matched)", color=mass_colors[mass], histtype="fill", linewidth=2)
    
    Efficiency_match = [num/den if den != 0 else 0 for num,den in zip(counts_match, counts_NGT)]
    
    centers = 0.5*(edges_NGT[:-1]+edges_NGT[1:])
    bin_widths = np.diff(edges_NGT)
    y_error_match = [np.sqrt(eff * (1 - eff) / den) if den > 0 else 0 for eff, den in zip(Efficiency_match, counts_NGT)] 
    
    ax2 = ax.twinx()
    ax2.errorbar(centers, Efficiency_match, xerr=0.5*bin_widths, yerr = y_error_match, color=eff_color, fmt = '*')
    ax2.set_ylabel("Efficiency", color=eff_color)
    ax2.tick_params(axis='y', labelcolor=eff_color)
    ax2.set_ylim(0,1)
    
    ax.set_xlabel(xlabel, fontsize=24)
    ax.set_ylabel("Entries", fontsize=24)
    ax.legend(loc="upper right")
    hep.cms.label("Private", data=True, ax=ax, com=14)

    fig.savefig(f"{path}{mass}_{var}.png")
    plt.close(fig)
    print(f"Saved Reco_Efficiency_{mass}_{var}.png")

def plot_efficiency(mass, data_per_mass, data_per_mass_NGT, var, bins, hist_range, xlabel):
    path = f"{plot_dir}/Efficiencies/"
    os.makedirs(path, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 8))

    arrs     = data_per_mass[mass]
    arrs_NGT = data_per_mass_NGT[mass]
    
    areas=[]
    max_count = 0

    counts, edges = np.histogram(arrs[var], bins=bins, range=hist_range, density=False)
    counts_NGT, edges_NGT = np.histogram(arrs_NGT[var], bins=bins, range=hist_range, density=False)

    area     = np.sum(counts * np.diff(edges))
    area_NGT = np.sum(counts_NGT * np.diff(edges_NGT))
    ratio = area_NGT/area*100

    hep.histplot(counts, edges,         ax=ax, label=f"{m} = {mass} GeV (All)" , color=mass_colors[mass], histtype="step", linewidth=2)
    hep.histplot(counts_NGT, edges_NGT, ax=ax, label=f"NGT: {ratio:.2f}%",       color=mass_colors[mass], histtype="fill", linewidth=2)
    
    Efficiency_NGT = [num/den if den != 0 else 0 for num,den in zip(counts_NGT, counts)]


    centers = 0.5*(edges[:-1]+edges[1:])
    bin_widths = np.diff(edges)
    y_error_NGT = [np.sqrt(eff * (1 - eff) / den) if den > 0 else 0 for eff, den in zip(Efficiency_NGT, counts)] 
    
    max_count = max(max(counts), max_count)

    ax2 = ax.twinx()
    ax2.errorbar(centers, Efficiency_NGT, xerr=0.5*bin_widths, yerr = y_error_NGT, color=eff_color, fmt = '*')
    ax2.set_ylabel("Efficiency", color=eff_color)
    ax2.tick_params(axis='y', labelcolor=eff_color)
    ax2.set_ylim(0,1)
    
    ax.set_xlabel(xlabel, fontsize=24)
    ax.set_ylabel("Entries", fontsize=24)
    ax.legend(loc="upper right", fontsize=18)

    hep.cms.label("Private work", loc=2, data=True,ax=ax, rlabel = f"{process} (0 PU) | 14 TeV", fontsize = 22)
    ax.grid(True, alpha=0.75, linestyle="dashdot", linewidth=0.75)
    ax.set_ylim(0, max_count * 1.3)

    fig.savefig(f"{path}{mass}_{var}.png")
    plt.close(fig)
    print(f"Saved Efficiency_{mass}_{var}.png")


def plot_genVsReco(mass, data_per_mass_match, var_gen, var_reco, bins, xlabel):
    path = f"{plot_dir}/GenVsReco/"
    os.makedirs(path, exist_ok=True)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    arrs_match = data_per_mass_match[mass]
        
    counts_gen, edges_gen = np.histogram(arrs_match[var_gen], bins=bins, range=(mass-1, mass+2), density=False)
    hep.histplot(counts_gen, edges_gen, ax=ax, label=f"GEN",
                    color=mass_colors[mass], histtype="fill", hatch= "", linewidth=2, linestyle="-")
    
    counts_reco, edges_reco = np.histogram(arrs_match[var_reco], bins=bins, range=(mass-1, mass+2), density=False)
    hep.histplot(counts_reco, edges_reco, ax=ax, label=f"RECO",
                    color=mass_colors[mass], histtype="step", hatch= "//", linewidth=2, linestyle="-")

    ax.set_xlabel(xlabel, fontsize=24)
    ax.set_ylabel("Entries", fontsize=24)
    ax.legend(loc="upper right", fontsize=16)
    hep.cms.label("Private", data=True, ax=ax, com=14)
    plt.tight_layout()
    fig.savefig(f"{path}{mass}_GenVsReco{var_gen}.png")
    plt.close(fig)
    print(f"Saved genVsReco_{var_gen}.png")

#─────────────────────────────────────── Fitting ────────────────────────────────────────
def plot_fit(mass, data_per_mass_match, var_gen, bins, xlabel, var_recoL1=None, var_recoNGT=None):
    path = f"{plot_dir}/Fit/"
    os.makedirs(path, exist_ok=True)
    
    fig, ax = plt.subplots(figsize=(10, 8))
        
    arrs     = data_per_mass[mass]
    arrs_NGT = data_per_mass_NGT[mass]
    arrs_match_L1 = data_per_mass_match_L1[mass]
    arrs_match_HLT = data_per_mass_match_HLT[mass]

    legend_handles = []
    legend_labels = []

    #GEN
    counts, edges = np.histogram(arrs_match_L1[var_gen], bins=bins, range=(mass-0.5, mass+0.5), density=False)
    gen_hist = hep.histplot(counts, edges, ax=ax, label=f"GEN {m} = {mass} GeV", color=color_all, histtype="step", hatch="", linewidth=2)

    popt, perr, resolution, resolution_err = fit_gaussian(counts, edges)
    A_fit, mu_fit, sigma_fit = popt
    mu_err, sigma_err = perr[1], perr[2]
    
    x_fit = np.linspace(edges[0], edges[-1], 500)
    y_fit = gaussian(x_fit, *popt)
    gen_fit_line, = ax.plot(x_fit, y_fit, color=color_all, linestyle="--", linewidth=2, label=f"$\\mu$={mu_fit:.3f}±{mu_err:.3f}\n$\\sigma$={sigma_fit:.3f}±{sigma_err:.3f}")

    fit_results = {
        "gen_mu": mu_fit, "gen_mu_err": mu_err,
        "gen_sigma": sigma_fit, "gen_sigma_err": sigma_err,
        "gen_res": resolution, "gen_res_err": resolution_err,

    }

    legend_handles += [gen_hist[0], gen_fit_line]
    legend_labels  += [f"GEN {m} = {mass} GeV", f"$\\mu$={mu_fit:.3f}±{mu_err:.3f}\n$\\sigma$={sigma_fit:.3f}±{sigma_err:.3f}"]
    
    #RECO L1   
    if var_recoNGT is None and var_recoL1 is not None:
        counts2, edges2 = np.histogram(arrs_match_L1[var_recoL1], bins=bins, range=(mass-0.5, mass+0.5), density=False)
        l1_hist = hep.histplot(counts2, edges2, ax=ax, label=f"L1 {m} = {mass} GeV", color=color_L1, histtype="step", hatch="", linewidth=2)

        popt2, perr2, resolution2, resolution_err2 = fit_gaussian(counts2, edges2)
        A_fit2, mu_fit2, sigma_fit2 = popt2
        mu_err2, sigma_err2 = perr2[1], perr2[2]

        x_fit2 = np.linspace(edges2[0], edges2[-1], 500)
        y_fit2 = gaussian(x_fit2, *popt2)
        l1_fit_line, = ax.plot(x_fit2, y_fit2, color=color_L1, linestyle="--", linewidth=2, label=f"$\\mu$={mu_fit2:.3f}±{mu_err2:.3f}\n$\\sigma$={sigma_fit2:.3f}±{sigma_err2:.3f}")

        fit_results.update({
            "reco_mu": mu_fit2, "reco_mu_err": mu_err2,
            "reco_sigma": sigma_fit2, "reco_sigma_err": sigma_err2,
            "reco_res": resolution2, "reco_res_err": resolution_err2,

        })

        legend_handles += [l1_hist[0], l1_fit_line]
        legend_labels  += [f"L1 {m} = {mass} GeV", f"$\\mu$={mu_fit2:.3f}±{mu_err2:.3f}\n$\\sigma$={sigma_fit2:.3f}±{sigma_err2:.3f}"]

    #RECO NGT   
    if var_recoNGT is not None and var_recoL1 is None:
        counts3, edges3 = np.histogram(arrs_match_HLT[var_recoNGT], bins=bins, range=(mass-0.5, mass+0.5), density=False)
        ngt_hist = hep.histplot(counts3, edges3, ax=ax, label=f"HLT {m} = {mass} GeV", color=color_NGT, histtype="step", hatch="", linewidth=2)

        popt3, perr3, resolution3, resolution_err3 = fit_gaussian(counts3, edges3)
        A_fit3, mu_fit3, sigma_fit3 = popt3
        mu_err3, sigma_err3 = perr3[1], perr3[2]

        x_fit3 = np.linspace(edges3[0], edges3[-1], 500)
        y_fit3 = gaussian(x_fit3, *popt3)
        ngt_fit_line, = ax.plot(x_fit3, y_fit3, color=color_NGT, linestyle="--", linewidth=2, label=f"$\\mu$={mu_fit3:.3f}±{mu_err3:.3f}\n$\\sigma$={sigma_fit3:.3f}±{sigma_err3:.3f}")

        fit_results.update({
            "reco_mu": mu_fit3, "reco_mu_err": mu_err3,
            "reco_sigma": sigma_fit3, "reco_sigma_err": sigma_err3,
            "reco_res": resolution3, "reco_res_err": resolution_err3,

        })

        legend_handles += [ngt_hist[0], ngt_fit_line]
        legend_labels  += [f"HLT {m} = {mass} GeV", f"$\\mu$={mu_fit3:.3f}±{mu_err3:.3f}\n$\\sigma$={sigma_fit3:.3f}±{sigma_err3:.3f}"]

    #RECO both   
    if var_recoL1 and var_recoNGT is not None:
        counts1, edges1 = np.histogram(arrs_match_L1[var_recoL1], bins=bins, range=(mass-0.5, mass+0.5), density=False)
        l1_hist = hep.histplot(counts1, edges1, ax=ax, label=f"L1 {m} = {mass} GeV", color=color_L1, histtype="step", hatch="", linewidth=2)

        popt1, perr1, resolution1, resolution_err1 = fit_gaussian(counts1, edges1)
        A_fit1, mu_fit1, sigma_fit1 = popt1
        mu_err1, sigma_err1 = perr1[1], perr1[2]

        x_fit1 = np.linspace(edges1[0], edges1[-1], 500)
        y_fit1 = gaussian(x_fit1, *popt1)
        l1_fit_line, = ax.plot(x_fit1, y_fit1, color=color_L1, linestyle="--", linewidth=2, label=f"$\\mu$={mu_fit1:.3f}±{mu_err1:.3f}\n$\\sigma$={sigma_fit1:.3f}±{sigma_err1:.3f}")

        legend_handles += [l1_hist[0], l1_fit_line]
        legend_labels  += [f"L1 {m} = {mass} GeV", f"$\\mu$={mu_fit1:.3f}±{mu_err1:.3f}\n$\\sigma$={sigma_fit1:.3f}±{sigma_err1:.3f}"]

        counts2, edges2 = np.histogram(arrs_match_HLT[var_recoNGT], bins=bins, range=(mass-0.5, mass+0.5), density=False)
        ngt_hist = hep.histplot(counts2, edges2, ax=ax, label=f"HLT {m} = {mass} GeV", color=color_NGT, histtype="step", hatch="", linewidth=2)

        popt2, perr2, resolution2, resolution_err2 = fit_gaussian(counts2, edges2)
        A_fit2, mu_fit2, sigma_fit2 = popt2
        mu_err2, sigma_err2 = perr2[1], perr2[2]

        x_fit2 = np.linspace(edges2[0], edges2[-1], 500)
        y_fit2 = gaussian(x_fit2, *popt2)
        ngt_fit_line, = ax.plot(x_fit2, y_fit2, color=color_NGT, linestyle="--", linewidth=2, label=f"$\\mu$={mu_fit2:.3f}±{mu_err2:.3f}\n$\\sigma$={sigma_fit2:.3f}±{sigma_err2:.3f}")

        fit_results.update({
            "L1_mu": mu_fit1, "L1_mu_err": mu_err1,
            "L1_sigma": sigma_fit1, "L1_sigma_err": sigma_err1,
            "L1_res": resolution1, "L1_res_err": resolution_err1,
            "NGT_mu": mu_fit2, "NGT_mu_err": mu_err2,
            "NGT_sigma": sigma_fit2, "NGT_sigma_err": sigma_err2,
            "NGT_res": resolution2, "NGT_res_err": resolution_err2,

        })

        legend_handles += [ngt_hist[0], ngt_fit_line]
        legend_labels  += [f"HLT {m} = {mass} GeV", f"$\\mu$={mu_fit2:.3f}±{mu_err2:.3f}\n$\\sigma$={sigma_fit2:.3f}±{sigma_err2:.3f}"]

    ax.set_xlabel(xlabel, fontsize=24)
    ax.set_ylabel("Entries", fontsize=24)

    handles, labels = ax.get_legend_handles_labels()
    label_to_handle = dict(zip(labels, handles))

    ordered_labels = [lab for lab in legend_labels if lab in label_to_handle]
    ordered_handles = [label_to_handle[lab] for lab in ordered_labels]

    ax.legend(handles=ordered_handles, labels=ordered_labels, loc="upper right", fontsize=18)
    hep.cms.label("Private work", loc=2, data=True,ax=ax, rlabel = f"{process} (0 PU) | 14 TeV", fontsize = 22)
    ax.grid(True, alpha=0.75, linestyle="dashdot", linewidth=0.75)
    
    fig.savefig(f"{path}{mass}_fitGenVsReco.png")
    plt.close(fig)
    print(f"Saved {mass}_fit_GenVsReco.png")

    return fit_results
 

resolution_gen_L1 = {}
resolution_gen_HLT = {}
resolution_L1 = {}
resolution_HLT = {}
def plot_resolution(data_per_mass_match, var_gen, var_L1, var_NGT):
    path = f"{plot_dir}/Resolution/"
    os.makedirs(path, exist_ok=True)
    
    fig, ax = plt.subplots(figsize=(10, 8))

    for mass in data_per_mass:
        arrs     = data_per_mass[mass]
        arrs_NGT = data_per_mass_NGT[mass]
        arrs_match_L1 = data_per_mass_match_L1[mass]
        arrs_match_HLT = data_per_mass_match_HLT[mass]

        counts, edges = np.histogram(arrs_match_L1[var_gen], bins=100, range=(mass-0.5, mass+0.5), density=False)
        counts1, edges1 = np.histogram(arrs_match_HLT[var_gen], bins=100, range=(mass-0.5, mass+0.5), density=False)
        
        popt, perr, resolution, resolution_err = fit_gaussian(counts, edges)
        popt1, perr1, resolution1, resolution_err1 = fit_gaussian(counts1, edges1)
        
        A_fit, mu_fit, sigma_fit = popt
        A_fit1, mu_fit1, sigma_fit1 = popt1

        print(f"Mean GEN (L1) = {mu_fit}\nsigma = {sigma_fit}")
        print(f"Mean GEN (HLT) = {mu_fit1}\nsigma = {sigma_fit1}")
        resolution_gen_L1[mass] = resolution
        resolution_gen_HLT[mass] = resolution1

        plt.errorbar(mass, resolution_gen_L1[mass], yerr=resolution_err, color=color_all, marker="o", markersize="11", label=f"Gen {m} (L1 matched)")

        plt.errorbar(mass, resolution_gen_HLT[mass], yerr=resolution_err1, color=color_all, marker="s", markersize="11", label=f"Gen {m} (HLT matched)")
        
        popt2, perr2 = None, None
        popt3, perr3 = None, None
        
        counts2, edges2 = np.histogram(arrs_match_L1[var_L1],  bins=100, range=(mass-0.5, mass+0.5), density=False)
        counts3, edges3 = np.histogram(arrs_match_HLT[var_NGT], bins=100, range=(mass-0.5, mass+0.5), density=False)
        
        popt2, perr2, resolution2, resolution_err2 = fit_gaussian(counts2, edges2)
        popt3, perr3, resolution3, resolution_err3 = fit_gaussian(counts3, edges3)
        
        A_fit2, mu_fit2, sigma_fit2 = popt2
        A_fit3, mu_fit3, sigma_fit3 = popt3

        print(f"Mean L1 = {mu_fit2}\nsigma = {sigma_fit2}")
        print(f"Mean HLT = {mu_fit3}\nsigma = {sigma_fit3}")
        resolution_L1[mass] = resolution2
        resolution_HLT[mass] = resolution3

        plt.errorbar(mass, resolution_L1[mass],  yerr=resolution_err2, color=color_L1,  marker="o", markersize="11", label="L1 m$_{A'}$")
        plt.errorbar(mass, resolution_HLT[mass], yerr=resolution_err3, color=color_NGT, marker="o", markersize="11", label="HLT m$_{A'}$")


    ax.set_xlabel(r"m$_{A'}$ [GeV]", fontsize=24)
    ax.set_ylabel(r"Resolution $\sigma / \mu$", fontsize=24)
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))  #Duplicate labels overwrite each other, keeping only one
    ax.legend(by_label.values(), by_label.keys(), loc="upper right", fontsize=18)
    hep.cms.label("Private work", loc=0, data=True,ax=ax, rlabel = f"{process} (0 PU) | 14 TeV", fontsize = 22)
    ax.grid(True, alpha=0.75, linestyle="dashdot", linewidth=0.75)

    fig.savefig(f"{path}{var_gen}.png")
    plt.close(fig)
    print(f"Saved resolution_{var_gen}.png")

#─────────────────────────────────────── Efficiency table ────────────────────────────────────────

def save_efficiency_table(efficiency_table, output_dir):
    """Write a nicely formatted .txt table of total entries / NGT efficiency / matching efficiency per mass."""
    table_path = os.path.join(output_dir, "efficiency_table.txt")

    header = f"{'mX (GeV)':>10} {'Total Entries':>15} {'NGT Efficiency (%)':>20} {'Matching Eff. L1 (%)':>20} {'Matching Eff. HLT (%)':>20}\n"
    sep = "=" * (len(header) - 1) + "\n"

    with open(table_path, "w") as f:
        f.write(sep)
        f.write(header)
        f.write(sep)
        for mass in sorted(efficiency_table):
            row = efficiency_table[mass]
            f.write(f"{mass:>10} {row['total']:>15} {row['ngt_eff']:>20.2f} {row['match_effL1']:>20.2f} {row['match_effHLT']:>20.2f}\n")
        f.write(sep)

    print(f"Saved efficiency table to {table_path}")

#─────────────────────────────────────── Fit table ────────────────────────────────────────

def save_fit_table(fit_table, output_dir):
    table_path = os.path.join(output_dir, "fit_table.txt")

    header = (f"{'m$(A`)$ [GeV]':>10} "
              f"{'Gen mu':>12} {'Gen mu err':>12} {'Gen sigma':>12} {'Gen sigma err':>15} {'Gen res':>15} {'Gen res err':>13} "
              f"{'HLT mu':>12} {'HLT mu err':>12} {'HLT sigma':>12} {'HLT sigma err':>16} {'HLT res':>16} {'HLT res err':>14}\n")
    sep = "-" * (len(header) - 1) + "\n"

    with open(table_path, "w") as f:
        f.write(header)
        f.write(sep)
        for mass in sorted(fit_table):
            row = fit_table[mass]
            f.write(
                f"{mass:>10} "
                f"{row.get('gen_mu', float('nan')):>12.4f} {row.get('gen_mu_err', float('nan')):>12.4f} "
                f"{row.get('gen_sigma', float('nan')):>12.4f} {row.get('gen_sigma_err', float('nan')):>15.4f} "
                f"{row.get('gen_res', float('nan')):>15.4f} {row.get('gen_res_err', float('nan')):>13.4f} "
                f"{row.get('reco_mu', float('nan')):>12.4f} {row.get('reco_mu_err', float('nan')):>12.4f} "
                f"{row.get('reco_sigma', float('nan')):>12.4f} {row.get('reco_sigma_err', float('nan')):>16.4f} "
                f"{row.get('reco_res', float('nan')):>16.4f} {row.get('reco_res_err', float('nan')):>14.4f}\n")
    print(f"Saved fit table to {table_path}")

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

    to_plot = [ "nGen_muons", "Gen_muon_leading_pt", "Gen_muon_leading_eta", "Gen_muon_leading_phi", "Gen_muon_leading_mass",  "Gen_muon_dR",
                "Gen_muon_subleading_pt", "Gen_muon_subleading_eta", "Gen_muon_subleading_phi", "Gen_muon_subleading_mass",

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
    data_per_mass_match_L1 = {}
    data_per_mass_match_HLT = {}
    efficiency_table = {}
    for mass, base_dir in mass_dirs.items():
        df = get_dataframe(base_dir)
        df_NGT = df.Filter("DST_PFScouting == true")
        df_match_L1 = df.Filter("L1_muon[0] >= 0 && L1_muon[1] >= 0")
        df_match_HLT = df.Filter("Reco_muon[0] >= 0 && Reco_muon[1] >= 0")

        total_entries = df.Count().GetValue()
        triggered_entries = df_NGT.Count().GetValue()
        ngt_eff = overall_efficiency(df, df_NGT)
        match_effL1 = overall_efficiency(df, df_match_L1)
        match_effHLT = overall_efficiency(df, df_match_HLT)

        # print(f"Triggered entries for mX = {mass} GeV: ", triggered_entries)
        # print(f"NGT Efficiency for mX = {mass} GeV: ", ngt_eff)
        # print(f"L1 Matching Efficiency for mX = {mass} GeV: ", match_effL1)
        # print(f"HLT Matching Efficiency for mX = {mass} GeV: ", match_effHLT)

        efficiency_table[mass] = {
            "total": total_entries,
            "triggered": triggered_entries,
            "ngt_eff": ngt_eff * 100,
            "match_effL1": match_effL1 * 100,
            "match_effHLT": match_effHLT * 100,
        }

        column = df.AsNumpy(to_plot)
        column_NGT = df_NGT.AsNumpy(to_plot)
        column_match_L1 = df_match_L1.AsNumpy(to_plot)
        column_match_HLT = df_match_HLT.AsNumpy(to_plot)
        
        data_per_mass[mass] = {variable: flatten(column[variable]) for variable in to_plot}
        data_per_mass_NGT[mass] = {variable: flatten(column_NGT[variable]) for variable in to_plot}
        data_per_mass_match_L1[mass] = {variable: flatten(column_match_L1[variable]) for variable in to_plot}
        data_per_mass_match_HLT[mass] = {variable: flatten(column_match_HLT[variable]) for variable in to_plot}

    save_efficiency_table(efficiency_table, args.output_dir)

    # plot(data_per_mass, "nGen_muons",            bins=10, hist_range=(0, 10),          xlabel=r"Gen number of $\mu$",   save_as='png', scale='linear')
    # plot(data_per_mass, "Gen_muon_leading_pt",      bins=30, hist_range=(0, 60),          xlabel=r"Gen p$_T$ Leading $\mu$   [GeV]", save_as='png', scale='linear')
    # plot(data_per_mass, "Gen_muon_subleading_pt",   bins=30, hist_range=(0, 60),          xlabel=r"Gen p$_T$ Subeading $\mu$ [GeV]", save_as='png', scale='linear')
    # plot(data_per_mass, "Gen_muon_leading_eta",  bins=30, hist_range=(-4, 6),          xlabel=r"Gen $\eta(\mu)$", save_as='png', scale='linear')
    # plot(data_per_mass, "Gen_muon_leading_phi",  bins=30, hist_range=(-np.pi, np.pi),  xlabel=r"Gen $\phi(\mu)$", save_as='png', scale='linear')
    
    # plot(data_per_mass, "Gen_muon_dR",   bins=30, hist_range=(0, 1),           xlabel=r"Gen $\Delta R(\mu_1, \mu2)$", save_as='png', scale='linear')

    # plot(data_per_mass, "Gen_X_pt",   bins=30, hist_range=(0, 80),          xlabel=r"Gen p$_T$ X [GeV]", save_as='png', scale='linear')
    # plot(data_per_mass, "Gen_X_eta",  bins=30, hist_range=(-4, 6),          xlabel=r"Gen $\eta(X)$", save_as='png', scale='linear')
    # plot(data_per_mass, "Gen_X_phi",  bins=30, hist_range=(-np.pi, np.pi),  xlabel=r"Gen $\phi(X)$", save_as='png', scale='linear')
    # plot(data_per_mass, "Gen_X_M",    bins=30, hist_range=(0, 15),          xlabel=r"Gen $m_{A'}$ [GeV]", save_as='png', scale='linear')
    # ##

    # plot(data_per_mass, "nL1_muon", bins=10, hist_range=(0, 10),            xlabel=r"Reco number of $\mu$", save_as='png', scale='linear')

    # plot(data_per_mass, "L1_muon_pt",  bins=30, hist_range=(0, 60),         xlabel=r"L1 p$_T$ $\mu$ [GeV]", save_as='png', scale='linear')
    # plot(data_per_mass, "L1_muon_eta", bins=30, hist_range=(-3, 3),         xlabel=r"L1 $\eta(\mu)$", save_as='png', scale='linear')
    # plot(data_per_mass, "L1_muon_phi", bins=30, hist_range=(-np.pi, np.pi), xlabel=r"L1 $\phi(\mu)$", save_as='png', scale='linear')

    # plot(data_per_mass, "L1_X_pt",  bins=30, hist_range=(0, 80),            xlabel=r"L1 p$_T$ X [GeV]", save_as='png', scale='linear')
    # plot(data_per_mass, "L1_X_eta", bins=30, hist_range=(-3, 3),            xlabel=r"L1 $\eta(X)$", save_as='png', scale='linear')
    # plot(data_per_mass, "L1_X_phi", bins=30, hist_range=(-np.pi, np.pi),    xlabel=r"L1 $\phi(X)$", save_as='png', scale='linear')
    # plot(data_per_mass, "L1_X_M",   bins=30, hist_range=(0, 15),            xlabel=r"L1 $m_{A'}$ [GeV]", save_as='png', scale='linear')

    # ##
    # plot(data_per_mass, "nReco_muon", bins=10, hist_range=(0, 10),             xlabel=r"Reco number of $\mu$", save_as='png', scale='linear')

    # plot(data_per_mass, "Reco_muon_pt",  bins=30, hist_range=(0, 60),          xlabel=r"Reco p$_T$ $\mu$ [GeV]", save_as='png', scale='linear')
    # plot(data_per_mass, "Reco_muon_eta", bins=30, hist_range=(-3, 3),          xlabel=r"Reco $\eta(\mu)$", save_as='png', scale='linear')
    # plot(data_per_mass, "Reco_muon_phi", bins=30, hist_range=(-np.pi, np.pi),  xlabel=r"Reco $\phi(\mu)$", save_as='png', scale='linear')

    # plot(data_per_mass, "Reco_X_pt",  bins=30, hist_range=(0, 80),          xlabel=r"Reco p$_T$ X [GeV]", save_as='png', scale='linear')
    # plot(data_per_mass, "Reco_X_eta", bins=30, hist_range=(-3, 3),          xlabel=r"Reco $\eta(X)$", save_as='png', scale='linear')
    # plot(data_per_mass, "Reco_X_phi", bins=30, hist_range=(-np.pi, np.pi),  xlabel=r"Reco $\phi(X)$", save_as='png', scale='linear')
    # plot(data_per_mass, "Reco_X_M",   bins=30, hist_range=(0, 15),          xlabel=r"Reco $m_{A'}$ [GeV]", save_as='png', scale='linear')

    # plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match, "Gen_muon_dR", bins=30, hist_range=(0, 1),  xlabel=r"Gen $\Delta R$ X [GeV]")
    # plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match_L1, data_per_mass_match_HLT, "Gen_muon_leading_pt",    bins=30, hist_range=(0, 60), xlabel=r"Gen p$_T$ Leading $\mu$   [GeV]")
    # plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match_L1, data_per_mass_match_HLT, "Gen_muon_subleading_pt", bins=30, hist_range=(0, 60), xlabel=r"Gen p$_T$ Subleading $\mu$   [GeV]")
    # plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match_L1, data_per_mass_match_HLT, "L1_muon_pt",             bins=30, hist_range=(0, 60), xlabel=r"L1 p$_T$ $\mu$ [GeV]")
    # plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match_L1, data_per_mass_match_HLT, "Reco_muon_pt",           bins=30, hist_range=(0, 60), xlabel=r"NGT p$_T$ $\mu$ [GeV]")

    # plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match, "Gen_X_pt",    bins=30, hist_range=(0, 80), xlabel=r"Gen p$_T$ X [GeV]")
    # plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match, "Gen_X_eta",   bins=30, hist_range=(-4, 6), xlabel=r"Gen $\eta(X)$")
    # plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match_L1, data_per_mass_match_HLT, "Gen_X_M", bins=30, hist_range=(0.5, 10.5), xlabel=r"Gen $m_{A'}$ [GeV]")
    
    plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match_L1, data_per_mass_match_HLT, "Gen_X_M",     bins=100, hist_range=(0, 12), xlabel=r"Gen $m_{A'}$ [GeV]")

    # plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match, "Reco_muon_dR", bins=30, hist_range=(0, 1),  xlabel=r"Reco $\Delta R$ X [GeV]")
    # plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match, "Reco_X_pt",    bins=30, hist_range=(0, 80), xlabel=r"Reco p$_T$ X [GeV]")
    # plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match, "Reco_X_eta",   bins=30, hist_range=(-4, 6), xlabel=r"Reco $\eta(X)$")
    # plot_comparison(data_per_mass, data_per_mass_NGT, data_per_mass_match, "Reco_X_M",     bins=30, hist_range=(0, 15), xlabel=r"Reco $m_{A'}$ [GeV]")

    # for mass in mass_dirs.keys():
        #plot_efficiency(mass, data_per_mass, data_per_mass_NGT, "Gen_muon_dR",               bins=30, hist_range=(0, 4),   xlabel=r"Gen $\Delta R(\mu_1, \mu2)$")
        # plot_efficiency(mass, data_per_mass, data_per_mass_NGT, "Gen_muon_leading_pt",       bins=30, hist_range=(0, 60),  xlabel=r"Gen Leading p$_T (\mu)$ [GeV]")
        # plot_efficiency(mass, data_per_mass, data_per_mass_NGT, "Gen_muon_subleading_pt",    bins=30, hist_range=(0, 60),  xlabel=r"Gen Subleading p$_T (\mu)$ [GeV]")
    #     plot_efficiency(mass, data_per_mass, data_per_mass_NGT, "Gen_muon_leading_eta",      bins=30, hist_range=(-4, 6),          xlabel=r"Gen $\eta($\mu$)$",                eff_color="#CE8D36")
    #     plot_efficiency(mass, data_per_mass, data_per_mass_NGT, "Gen_muon_leading_phi",      bins=30, hist_range=(-np.pi, np.pi),  xlabel=r"Gen $\phi($\mu$)$",                eff_color="#CE8D36")
    #     plot_efficiency(mass, data_per_mass, data_per_mass_NGT, "Gen_muon_leading_mass",     bins=30, hist_range=(mass-1, mass+2), xlabel=r"Gen $m_{\mu}$ [GeV]",    eff_color="#CE8D36")

    #     plot_efficiency(mass, data_per_mass, data_per_mass_NGT, "Gen_X_pt",    bins=30, hist_range=(0, 80),          xlabel=r"Gen p$_T$ X [GeV]",           eff_color="#CE8D36")
    #     plot_efficiency(mass, data_per_mass, data_per_mass_NGT, "Gen_X_eta",   bins=30, hist_range=(-4, 6),          xlabel=r"Gen $\eta(X)$",               eff_color="#CE8D36")
    #     plot_efficiency(mass, data_per_mass, data_per_mass_NGT, "Gen_X_phi",   bins=30, hist_range=(-np.pi, np.pi),  xlabel=r"Gen $\phi(X)$",               eff_color="#CE8D36")
    #     plot_efficiency(mass, data_per_mass, data_per_mass_NGT, "Gen_X_M",     bins=30, hist_range=(mass-1, mass+2),           xlabel=r"Gen $m_{A'}$ [GeV]",   eff_color="#CE8D36")

    #     plot_recoefficiency(mass, data_per_mass_NGT, data_per_mass_match, "Gen_X_M", bins=30, hist_range=(mass-1, mass+2),    xlabel=r"Gen $m_{A'}$ [GeV]",   eff_color="#CE8D36")

    # plot_genVsReco(1, data_per_mass_match, "Gen_X_M", "Reco_X_M", bins=30, xlabel=r"$m_{A'}$ [GeV]")
    # plot_genVsReco(2, data_per_mass_match, "Gen_X_M", "Reco_X_M", bins=30, xlabel=r"$m_{A'}$ [GeV]")
    # plot_genVsReco(4, data_per_mass_match, "Gen_X_M", "Reco_X_M", bins=30, xlabel=r"$m_{A'}$ [GeV]")
    # plot_genVsReco(6, data_per_mass_match, "Gen_X_M", "Reco_X_M", bins=30, xlabel=r"$m_{A'}$ [GeV]")
    # plot_genVsReco(8, data_per_mass_match, "Gen_X_M", "Reco_X_M", bins=30, xlabel=r"$m_{A'}$ [GeV]")
    # plot_genVsReco(10, data_per_mass_match, "Gen_X_M", "Reco_X_M", bins=30, xlabel=r"$m_{A'}$ [GeV]")

    # fit_table = {}
    # for mass in [1, 2, 4, 6, 8, 10]:
        # fit_table[mass] = plot_fit(mass, data_per_mass_match_L1, "Gen_X_M", bins=100, xlabel=r"m$_{A'}$ [GeV]", var_recoL1="L1_X_M", var_recoNGT=None)
    # save_fit_table(fit_table, args.output_dir)

    # plot_resolution(data_per_mass_match_L1, "Gen_X_M", var_L1="L1_X_M", var_NGT="Reco_X_M")