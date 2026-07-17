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

ROOT.gROOT.SetBatch(True)
hep.style.use(hep.style.CMS)

SRC = "/eos/user/s/sbenabde/CERN_Summer_student/X_mumu/functions.cc"
TMP = "/tmp/functions.cc"

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

#─────────────────────────────────────── Define variables ────────────────────────────────────────

def apply_defines(df):
    return (
        df
        .Define("Gen_muon_idx",             "get_muon_idx(GenPart_pdgId, GenPart_statusFlags)")
        .Define("nGen_muons",               "Gen_muon_idx.size()")
        .Define("Gen_muon_all_dR",          "get_dR(Gen_muon_idx, GenPart_eta, GenPart_phi)")

        .Define("Gen_muon_pair_indices",    "deltaR_muon_pairing(Gen_muon_idx, GenPart_eta, GenPart_phi, GenPart_genPartIdxMother)")
        .Define("Gen_muon_pair_M",          "get_pair_M(Gen_muon_pair_indices, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
        .Define("nGen_muons_pairs",         "Gen_muon_pair_indices.size()")

        .Define("Gen_muon_pt",              "get_pair_variable(Gen_muon_pair_indices, GenPart_pt)")
        .Define("Gen_muon_eta",             "get_pair_variable(Gen_muon_pair_indices, GenPart_eta)")
        .Define("Gen_muon_phi",             "get_pair_variable(Gen_muon_pair_indices, GenPart_phi)")
        .Define("Gen_muon_mass",            "get_pair_variable(Gen_muon_pair_indices, GenPart_mass)")
        .Define("Gen_muon_dR",              "get_dR(Gen_muon_pair_indices, GenPart_eta, GenPart_phi)")

        .Define("Gen_X_pt",                 "get_pair_pt(Gen_muon_pair_indices, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
        .Define("Gen_X_eta",                "get_pair_eta(Gen_muon_pair_indices, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
        .Define("Gen_X_phi",                "get_pair_phi(Gen_muon_pair_indices, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
        .Define("Gen_X_M",                  "get_pair_M(Gen_muon_pair_indices, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")

        .Define("Reco_muon",                "truth_matching(Gen_muon_pair_indices, GenPart_eta, GenPart_phi, hltMuon_pt, hltMuon_eta, hltMuon_phi)")
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

def plot(data_per_mass, var, bins, hist_range, xlabel, outname):
    path = "/eos/user/s/sbenabde/CERN_Summer_student/X_mumu/plots/"
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
    fig.savefig(f"{path}{var}.png")
    plt.close(fig)
    print(f"Saved {var}.png")

def plot_comparison(data_per_mass, data_per_mass_NGT, var, bins, hist_range, xlabel):
    path = "/eos/user/s/sbenabde/CERN_Summer_student/X_mumu/plots/Comparisons/"
    os.makedirs(path, exist_ok=True)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    for mass in data_per_mass:
        arrs     = data_per_mass[mass]
        arrs_NGT = data_per_mass_NGT[mass]

        counts, edges = np.histogram(arrs[var], bins=bins, range=hist_range, density=False)
        hep.histplot(counts, edges, ax=ax, label=f"m$_X$ = {mass} GeV (All)",
                     color=mass_colors[mass], histtype="step", hatch= "", linewidth=2, linestyle="-")

        counts_NGT, edges_NGT = np.histogram(arrs_NGT[var], bins=bins, range=hist_range, density=False)
        hep.histplot(counts_NGT, edges_NGT, ax=ax, label=f"m$_X$ = {mass} GeV (NGT)",
                     color=mass_colors[mass], histtype="step", hatch= "/", linewidth=2, linestyle="--")

    ax.set_xlabel(xlabel, fontsize=24)
    ax.set_ylabel("Events", fontsize=24)
    ax.legend(loc="upper right", fontsize=14)
    hep.cms.label("Preliminary", data=False, ax=ax, com=14)
    plt.tight_layout()
    fig.savefig(f"{path}{var}.png")
    plt.close(fig)
    print(f"Saved Comparison_{var}.png")

def plot_efficiency(mass, data_per_mass, data_per_mass_NGT, var, bins, hist_range, xlabel, eff_color):
    path = "/eos/user/s/sbenabde/CERN_Summer_student/X_mumu/plots/Efficiencies/"
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

#─────────────────────────────────────── Main ────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Argument for Analysis_all')
    parser.add_argument('--run-comparison', action='store_true', required=False, help='Run comparison plots across all generated X masses.')
    args = parser.parse_args()

    os.makedirs("plots", exist_ok=True)

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

    to_plot = [ "Gen_muon_pt",
                "Gen_muon_eta", 
                "Gen_muon_phi", 
                "Gen_muon_mass", 
                "Gen_muon_dR",
                
                "Gen_X_pt",
                "Gen_X_eta", 
                "Gen_X_phi", 
                "Gen_X_M", 

                "Reco_muon_pt", 
                "Reco_muon_eta", 
                "Reco_muon_phi", 

                "Reco_X_pt", 
                "Reco_X_eta", 
                "Reco_X_phi", 
                "Reco_X_M"
    ]

    data_per_mass = {}
    data_per_mass_NGT = {}
    for mass, base_dir in mass_dirs.items():
        df = get_dataframe(base_dir)
        df_NGT = df.Filter("DST_PFScouting == true")
        print(f"NGT Efficiency for mX = {mass} GeV: ", overall_efficiency(df, df_NGT))
        
        column = df.AsNumpy(to_plot)
        column_NGT = df_NGT.AsNumpy(to_plot)
        data_per_mass[mass] = {variable: flatten(column[variable]) for variable in to_plot}
        data_per_mass_NGT[mass] = {variable: flatten(column_NGT[variable]) for variable in to_plot}

    plot(data_per_mass, "Gen_muon_dR",bins=30, hist_range=(0, 1),           xlabel= "Generated $\Delta R$ X [GeV]",outname="/eos/user/s/sbenabde/CERN_Summer_student/X_mumu/plots/Gen_muon_dR")

    plot(data_per_mass, "Gen_X_pt",   bins=30, hist_range=(0, 80),          xlabel= "Generated p$_T$ X [GeV]",     outname="/eos/user/s/sbenabde/CERN_Summer_student/X_mumu/plots/Gen_X_pt")
    plot(data_per_mass, "Gen_X_eta",  bins=30, hist_range=(-4, 6),          xlabel= "Generated $\eta(X)$",         outname="/eos/user/s/sbenabde/CERN_Summer_student/X_mumu/plots/Gen_X_eta")
    plot(data_per_mass, "Gen_X_phi",  bins=30, hist_range=(-np.pi, np.pi),  xlabel= "Generated $\phi(X)$",         outname="/eos/user/s/sbenabde/CERN_Summer_student/X_mumu/plots/Gen_X_phi")
    plot(data_per_mass, "Gen_X_M",    bins=30, hist_range=(0, 12),          xlabel= "Generated $m_X$ [GeV]",       outname="/eos/user/s/sbenabde/CERN_Summer_student/X_mumu/plots/Gen_X_M")

    plot(data_per_mass, "Reco_X_pt",  bins=30, hist_range=(0, 80),          xlabel= "Reconstructed p$_T$ X [GeV]", outname="/eos/user/s/sbenabde/CERN_Summer_student/X_mumu/plots/Reco_X_pt")
    plot(data_per_mass, "Reco_X_eta", bins=30, hist_range=(-3, 3),          xlabel= "Reconstructed $\eta(X)$",     outname="/eos/user/s/sbenabde/CERN_Summer_student/X_mumu/plots/Reco_X_eta")
    plot(data_per_mass, "Reco_X_phi", bins=30, hist_range=(-np.pi, np.pi),  xlabel= "Reconstructed $\phi(X)$",     outname="/eos/user/s/sbenabde/CERN_Summer_student/X_mumu/plots/Reco_X_phi")
    plot(data_per_mass, "Reco_X_M",   bins=30, hist_range=(0, 15),          xlabel= "Reconstructed $m_X$ [GeV]",   outname="/eos/user/s/sbenabde/CERN_Summer_student/X_mumu/plots/Reco_X_mass")

    plot_comparison(data_per_mass, data_per_mass_NGT, "Gen_muon_dR", bins=30, hist_range=(0, 1),  xlabel= "Generated $\Delta R$ X [GeV]")
    plot_comparison(data_per_mass, data_per_mass_NGT, "Gen_X_pt",    bins=30, hist_range=(0, 80), xlabel= "Generated p$_T$ X [GeV]")
    plot_comparison(data_per_mass, data_per_mass_NGT, "Gen_X_eta",   bins=30, hist_range=(-4, 6), xlabel= "Generated $\eta(X)$")
    plot_comparison(data_per_mass, data_per_mass_NGT, "Gen_X_M",     bins=30, hist_range=(0, 12), xlabel= "Generated $m_X$ [GeV]")

    plot_efficiency(2, data_per_mass, data_per_mass_NGT, "Gen_X_pt",    bins=30, hist_range=(0, 80),          xlabel= "Generated p$_T$ X [GeV]",      eff_color="yellowgreen")
    plot_efficiency(2, data_per_mass, data_per_mass_NGT, "Gen_X_eta",   bins=30, hist_range=(-4, 6),          xlabel= "Generated $\eta(X)$",          eff_color="yellowgreen")
    plot_efficiency(2, data_per_mass, data_per_mass_NGT, "Gen_X_phi",   bins=30, hist_range=(-np.pi, np.pi),  xlabel= "Generated $\phi(X)$",          eff_color="yellowgreen")
    plot_efficiency(2, data_per_mass, data_per_mass_NGT, "Gen_X_M",     bins=30, hist_range=(1, 4),           xlabel= "Generated $m_X$ [GeV]",        eff_color="yellowgreen")
    plot_efficiency(2, data_per_mass, data_per_mass_NGT, "Gen_muon_dR", bins=30, hist_range=(0, 4),           xlabel= "Generated $\Delta R$ X [GeV]", eff_color="yellowgreen")


