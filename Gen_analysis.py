import ROOT
ROOT.gROOT.SetBatch(True)
import uproot
import glob
import os
import sys
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
hep.style.use(hep.style.CMS)

import shutil, os

SRC = "/eos/user/s/sbenabde/CERN_Summer_student/functions/functions.cc"
HDR = "/eos/user/s/sbenabde/CERN_Summer_student/headers/functions.h"
TMP = "/tmp/functions.cc"

shutil.copy(SRC, TMP)
ROOT.gROOT.ProcessLine(f'.L {TMP}+')

ROOT.ROOT.EnableImplicitMT()

#─────────────────────────────────────── Get files ────────────────────────────────────────

base_dir  = "/eos/user/e/evernazz/Sarah/HHbbtautau/HLTStandard/"
files_ngt = sorted(glob.glob(base_dir + "/step2_*.root"))
 
if not files_ngt:
    raise RuntimeError(f"No ROOT files found in {base_dir}")
 
chain = ROOT.TChain("Events")
for f in files_ngt:
    chain.Add(f)
 
df = ROOT.RDataFrame(chain)
print(f"Total entries in chain: {df.Count().GetValue()}")

#─────────────────────────────────────── Print column names ────────────────────────────────────────

# all_cols = [str(c) for c in df.GetColumnNames()]
# gen_cols = [c for c in all_cols if "GenPart" in c]
# print("\nGenPart columns found:")
# for c in sorted(gen_cols):
#     print(f"  {c}  [{df.GetColumnType(c)}]")


#─────────────────────────────────────── Define variables ────────────────────────────────────────

df = (
    df
    .Define("H_tautau_idx",          f"Get_H_idx(15, GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")
    .Define("H_bb_idx",              f"Get_H_idx(5, GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")

    .Define("H_tautau_mass",         f"GenPart_mass.at(H_tautau_idx)")
    .Define("H_bb_mass",             f"GenPart_mass.at(H_bb_idx)")
    .Define("HH_mass",               f"Get_mHH(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, H_tautau_idx, H_bb_idx)")

    .Define("tau_indices",           "Get_tau_from_H_indices(GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")
    .Define("tau_channel",           "Get_tau_channel(tau_indices, GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")

    .Define("tau_pt",                "Get_tau_pt(tau_indices, GenPart_pt)")
    .Define("tau_pt_lead",           "tau_pt[0]")
    .Define("tau_pt_sublead",        "tau_pt[1]")
)


df_tau_e_tau_h = df.Filter("tau_channel  == 0")
df_tau_m_tau_h = df.Filter("tau_channel == 1")
df_tau_h_tau_h = df.Filter("tau_channel  == 2 ")
df_tau_h_tau_h_triggered = df.Filter("tau_channel  == 2 && HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1 == true")

#─────────────────────────────────────── Define and plot histograms 1D ────────────────────────────────────────

separate_files = False                           #True gets 1 file with all the histograms, false gets a separate file for each histogram


data = {
    # "e_h":  df_tau_e_tau_h.AsNumpy(["HH_mass"])["HH_mass"],
    # "mu_h": df_tau_m_tau_h.AsNumpy(["HH_mass"])["HH_mass"],
    # "h_h":  df_tau_h_tau_h.AsNumpy(["HH_mass"])["HH_mass"],
    "h_h_pt_leading":  df_tau_h_tau_h.AsNumpy(["tau_pt_lead"])["tau_pt_lead"],
    "h_h_pt_leading_triggered":  df_tau_h_tau_h_triggered.AsNumpy(["tau_pt_lead"])["tau_pt_lead"],
    # "h_h_pt_subleading":  df_tau_h_tau_h.AsNumpy(["tau_pt_sublead"])["tau_pt_sublead"],
}

plot_configs = {
    # "e_h":  {"label": r"$\tau \tau \to \tau_e \tau_h$",   "bins": 50, "range": (200, 800), "color": "orchid", "histtype": "fill"},
    # "mu_h": {"label": r"$\tau \tau \to \tau_{\mu} \tau_h$","bins": 50, "range": (200, 800), "color": "darkorchid", "histtype": "fill"},
    # "h_h":  {"label": r"$\tau \tau \to \tau_h \tau_h$",  "bins": 50, "range": (200, 800), "color": "mediumvioletred", "histtype": "fill"},
    "h_h_pt_leading":  {"label": r"All events",  "bins": 50, "range": (0, 300), "color": "orchid", "histtype": "fill"},
    "h_h_pt_leading_triggered":  {"label": r"Events that pass trigger",  "bins": 50, "range": (0, 300), "color": "mediumvioletred", "histtype": "fill"},
    # "h_h_pt_subleading":  {"label": r"$\tau \tau \to \tau_h \tau_h$",  "bins": 50, "range": (0, 200), "color": "orchid", "histtype": "fill"},
}

def make_panel(ax, col, cfg):
    counts, edges = np.histogram(data[col], bins=cfg["bins"], range=cfg["range"])
    hep.histplot(counts, edges, ax=ax, label=cfg["label"], color=cfg["color"], histtype=cfg["histtype"])
    ax.set_xlabel(r"$p_T$ Leading \tau[GeV]")
    ax.set_ylabel("Events")
    ax.legend()
    hep.cms.label("Preliminary", data=False, ax=ax, com=14)

if separate_files:
    for col, cfg in plot_configs.items():
        fig, ax = plt.subplots(figsize=(10, 8))
        make_panel(ax, col, cfg)
        plt.tight_layout()
        plt.savefig(f"plots/{col}.png", dpi=300, bbox_inches="tight")
        plt.close()
else:
    fig, ax = plt.subplots(figsize=(10, 8))
    for col, cfg in plot_configs.items():
        counts, edges = np.histogram(data[col], bins=cfg["bins"], range=cfg["range"])
        hep.histplot(counts, edges, ax=ax, label=cfg["label"], color=cfg["color"], histtype=cfg["histtype"])
    ax.set_xlabel(r"$p_T$ Leading \tau[GeV]")
    ax.set_ylabel("Events")
    ax.legend()
    hep.cms.label("Preliminary", data=False, ax=ax, com=14)
    plt.tight_layout()
    plt.savefig("plots/All_Vs_Trigger.png", dpi=300, bbox_inches="tight")
    plt.close()


#─────────────────────────────────────── Define and plot histograms 2D ────────────────────────────────────────

separate_files_2D = True

data_2D = {
    "leading_vs_subleading_pt": df_tau_h_tau_h.AsNumpy(['tau_pt_lead', 'tau_pt_sublead']),
}

plot_configs_2D = {
    "leading_vs_subleading_pt":  {"label": r"$\tau \tau \to \tau_h \tau_h$",  "x": "tau_pt_lead", "y": "tau_pt_sublead", "bins": 50, "range": [[0, 300], [0, 300]], "cmap": "YlGn"},
}

def make_panel_2d(ax, col, cfg):
    h, xedges, yedges = np.histogram2d(data_2D[col][cfg["x"]], data_2D[col][cfg["y"]], bins=cfg["bins"], range=cfg["range"])
    mesh = ax.pcolormesh(xedges, yedges, h.T, cmap=cfg["cmap"])
    plt.colorbar(mesh, ax=ax, label="Events")
    ax.set_xlabel(r"$p_T Leading$ [GeV]")
    ax.set_ylabel(r"$p_T Subleading}$ [GeV]")
    # ax.set_title(cfg["label"])
    hep.cms.label("Preliminary", data=False, ax=ax, com=14)

if separate_files:
    for col, cfg in plot_configs_2D.items():
        fig, ax = plt.subplots(figsize=(10, 8))
        make_panel_2d(ax, col, cfg)
        plt.tight_layout()
        plt.savefig(f"plots/2d_{col}.png", dpi=300, bbox_inches="tight")
        plt.close()
else:
    fig, axes = plt.subplots(1, 3, figsize=(30, 8))
    for ax, (col, cfg) in zip(axes, plot_configs_2D.items()):
        make_panel_2d(ax, col, cfg)
    plt.tight_layout()
    plt.savefig("plots/2d_3_channels.png", dpi=300, bbox_inches="tight")
    plt.close()
#─────────────────────────────────────── Save root file output ────────────────────────────────────────

# HH_mass = df.Histo1D(("HH_mass", "HH_mass", 100, -10, 200), "HH_mass")

# out = ROOT.TFile("Signal_test.root", "RECREATE")
# for h in [HH_mass]:
#     h.Write()
# out.Close()
# print("\nSaved Signal_test.root")
