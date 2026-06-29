import ROOT
ROOT.gROOT.SetBatch(True)
import uproot
import glob
import os
import sys
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import argparse


hep.style.use(hep.style.CMS)

import shutil, os

SRC = "/eos/user/s/sbenabde/CERN_Summer_student/functions/functions.cc"
HDR = "/eos/user/s/sbenabde/CERN_Summer_student/headers/functions.h"
TMP = "/tmp/functions.cc"

shutil.copy(SRC, TMP)
ROOT.gROOT.ProcessLine(f'.L {TMP}+')

ROOT.ROOT.EnableImplicitMT()
df = ROOT.RDataFrame("Events", "Dataframes/df_ngt_hlt.root")

df = (
    df
    .Define("Gen_tau_idx",          f"Get_tau_from_H_indices(GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")
    
    .Define("Reco_tau",             f"get_4vector(hltHpsPFTau_pt, hltHpsPFTau_eta, hltHpsPFTau_phi, hltHpsPFTau_mass)")
    .Define("Matched_tau_idx",      f"deltaR_matching(Gen_tau_idx, GenPart_eta, GenPart_pt, Reco_tau, hltHpsPFTau_deepTauVSe, hltHpsPFTau_deepTauVSjet, hltHpsPFTau_deepTauVSmu)")
    
    .Define("Reco_tau_pt",          f"Get_tau_pt(Matched_tau_idx, hltHpsPFTau_pt)")
    # .Define("Reco_jet",    f"get_4vector(hltAK4PuppiJet_pt, hltAK4PuppiJet_eta, hltAK4PuppiJet_phi, hltAK4PuppiJet_mass)")
)



separate_files = False      #True gets 1 file with all the histograms, false gets a separate file for each histogram
data = {
    "Reco_tau_pt":  df.AsNumpy(["Reco_tau_pt"])["Reco_tau_pt"],
}

plot_configs = {
    "Reco_tau_pt":  {"label": r"Reco_tau_pt",  "xlabel" : r"$p_T \tau$", "ylabel": "Events", "bins": 30, "range": (0, 300), "color": "forestgreen", "histtype": "step", "hatch": "//"},
}

def make_panel(ax, col, cfg):
    flat_data = np.concatenate(data[col])
    counts, edges = np.histogram(flat_data, bins=cfg["bins"], range=cfg["range"])
    hep.histplot(counts, edges, ax=ax, label=cfg["label"], color=cfg["color"], hatch=cfg["hatch"], histtype=cfg["histtype"])
    ax.set_xlabel(cfg["xlabel"])
    ax.set_ylabel(cfg["ylabel"])
    ax.legend(loc="upper right")
    return counts, edges


fig, ax = plt.subplots(figsize=(10, 8))
for col, cfg in plot_configs.items():
    counts, edges = make_panel(ax, col, cfg)
    hep.cms.label("Preliminary", data=False, ax=ax, com=14)
    plt.tight_layout()
    plt.savefig("plots/Reco/tau_pt.png", dpi=300, bbox_inches="tight")
    plt.close()