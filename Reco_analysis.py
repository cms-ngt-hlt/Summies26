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

# ROOT.ROOT.EnableImplicitMT()
df = ROOT.RDataFrame("Events", "Dataframes/df_ngt_hlt.root")

df = (
    df
    .Define("Gen_tau_idx",          f"Get_tau_from_H_indices(GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")
    .Define("Gen_b_idx",            f"Get_b_from_H_indices(GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")
    
    .Define("Reco_tau",             f"get_4vector(hltHpsPFTau_pt, hltHpsPFTau_eta, hltHpsPFTau_phi, hltHpsPFTau_mass)")
    .Define("Matched_tau_idx",      f"deltaR_matching(Gen_tau_idx, GenPart_eta, GenPart_phi, Reco_tau, hltHpsPFTau_deepTauVSjet)")
    
    .Define("Reco_tau_pt",          f"Get_recotau_pt(Matched_tau_idx, hltHpsPFTau_pt)")

    .Define("DeltaR_tau",           f"Get_dR(Gen_tau_idx, Matched_tau_idx, GenPart_eta, GenPart_phi, hltHpsPFTau_eta, hltHpsPFTau_phi)")
    .Define("Delta_pT_tau",         f"Get_DpT(Gen_tau_idx, Matched_tau_idx, GenPart_pt, hltHpsPFTau_pt)")
    .Define("Delta_eta_tau",        f"Get_Deta(Gen_tau_idx, Matched_tau_idx, GenPart_eta, hltHpsPFTau_eta)")
    .Define("Delta_phi_tau",        f"Get_Dphi(Gen_tau_idx, Matched_tau_idx, GenPart_phi, hltHpsPFTau_phi)")
    
    .Define("Reco_jet",             f"get_4vector(hltAK4PuppiJet_pt, hltAK4PuppiJet_eta, hltAK4PuppiJet_phi, hltAK4PuppiJet_mass)")
    .Define("Matched_jet_idx",      f"deltaR_matching_jets(Gen_b_idx, GenPart_eta, GenPart_phi, Reco_jet, hltAK4PuppiJet_DeepFlavour_prob_b, hltAK4PuppiJet_DeepFlavour_prob_bb, hltAK4PuppiJet_DeepFlavour_prob_c, hltAK4PuppiJet_DeepFlavour_prob_g, hltAK4PuppiJet_DeepFlavour_prob_lepb, hltAK4PuppiJet_DeepFlavour_prob_uds)")
    
    .Define("DeltaR_jet",           f"Get_dR(Gen_b_idx, Matched_jet_idx, GenPart_eta, GenPart_phi, hltAK4PuppiJet_eta, hltAK4PuppiJet_phi)")
    .Define("Delta_pT_jet",         f"Get_DpT(Gen_b_idx, Matched_jet_idx, GenPart_pt, hltAK4PuppiJet_pt)")
    .Define("Delta_eta_jet",        f"Get_Deta(Gen_b_idx, Matched_jet_idx, GenPart_eta, hltAK4PuppiJet_eta)")
    .Define("Delta_phi_jet",        f"Get_Dphi(Gen_b_idx, Matched_jet_idx, GenPart_phi, hltAK4PuppiJet_phi)")

)

separate_files = False      #True gets 1 file with all the histograms, false gets a separate file for each histogram
data = {
    "Delta_phi_jet":  df.AsNumpy(["Delta_phi_jet"])["Delta_phi_jet"],
}

plot_configs = {
    "Delta_phi_jet":  {"label": r"$\Delta \phi$ between gen and reco jets",  "xlabel" : r"$\Delta \phi$", "ylabel": "Events", "bins": 30, "range": (0, 1), "color": "orchid", "histtype": "fill", "hatch": ""},
}

def make_panel(ax, col, cfg):
    flat_data = np.concatenate(data[col])
    flat_data = flat_data[flat_data >= 0]
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
    plt.savefig("plots/Reco/Delta_phi_jets.png", dpi=300, bbox_inches="tight")
    plt.close()

