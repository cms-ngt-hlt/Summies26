import ROOT
import uproot
import glob
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import argparse
import shutil, os
from pathlib import Path
import shutil
import tempfile
import mplhep as hep
hep.style.use("CMS")

REPO_DIR = Path(__file__).resolve().parent

SRC = REPO_DIR / "functions" / "functions.cc"
HDR_DIR = REPO_DIR / "headers"

tmp_dir = Path(tempfile.mkdtemp(prefix="root_aclic_"))
TMP = tmp_dir / "functions.cc"

shutil.copy2(SRC, TMP)
ROOT.gSystem.AddIncludePath(f"-I{HDR_DIR}")
ROOT.gROOT.ProcessLine(f".L {TMP}+")

#───────────────────────────────────────────────────────────────── Functions ────────────────────────────────────────────────────────────────────────────────

def get_leg_arrays(df_in, colnames):
    """Pull columns out of the RDataFrame. Vector-per-event columns are
    split into leg0/leg1 and flattened; scalar-per-event columns are
    passed through as-is."""
    cols = df_in.AsNumpy(colnames)
    out = {}
    for name in colnames:
        col = cols[name]
        if col.dtype == object:
            # vector-per-event column (e.g. RVecI/RVecF)
            leg0 = np.array([ev[0] for ev in col if len(ev) > 0])
            leg1 = np.array([ev[1] for ev in col if len(ev) > 1])
            out[name] = np.concatenate([leg0, leg1])
        else:
            # scalar-per-event column (e.g. nGen_Muon)
            out[name] = np.asarray(col)
    return out
    
def eff_and_err(num, den):
    eff = [n / d if d != 0 else 0 for n, d in zip(num, den)]
    err = np.array([np.sqrt(e * (1 - e) / d) if d > 0 else 0.0 for e, d in zip(eff, den)])
    return eff, err

#────────────────────────────────────────────────────────────────── Get dataframe ────────────────────────────────────────────────────────────────────────────────

df = ROOT.RDataFrame("Events", "/eos/user/s/sbenabde/CERN_Summer_student/Dataframes/df_ngt_hlt.root")

    #────────────────────────────────────────────────────────────── Define variables ────────────────────────────────────────────────────────────────────────────────

df = (
    df
    #─────────────────────────────────────────────────────────────────── Gen ────────────────────────────────────────────────────────────────────────────────

    .Define("Gen_b_idx",      f"Get_b_from_H_indices(GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")
    .Define("Gen_bJets_idx",  f"Match_b_to_GenJet(Gen_b_idx, GenPart_eta, GenPart_phi, GenJet_eta, GenJet_phi, GenJet_pt)")

    .Define("Gen_b_p4",       f"Get_b_p4(Gen_b_idx, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
    .Define("Gen_bjet_p4",    f"Get_b_p4(Gen_bJets_idx, GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass)")

    .Define("Gen_b_leading_jet_pt",       f"(Gen_bjet_p4[0].pt())")
    .Define("Gen_b_leading_jet_eta",      f"(Gen_bjet_p4[0].eta())")
    .Define("Gen_b_leading_jet_phi",      f"(Gen_bjet_p4[0].phi())")
    .Define("Gen_b_leading_jet_mass",     f"(Gen_bjet_p4[0].M())")

    .Define("Gen_b_subleading_jet_pt",    f"(Gen_bjet_p4[1].pt())")
    .Define("Gen_b_subleading_jet_eta",   f"(Gen_bjet_p4[1].eta())")
    .Define("Gen_b_subleading_jet_phi",   f"(Gen_bjet_p4[1].phi())")
    .Define("Gen_b_subleading_jet_mass",  f"(Gen_bjet_p4[1].M())")

    .Define("Gen_H_bb_mass",      f"Get_Gen_H_tautau_mass(Gen_b_p4)") 
    .Define("Gen_H_bbjet_mass",   f"Get_Gen_H_tautau_mass(Gen_bjet_p4)") 



    .Define("Gen_taus_idx",       f"Get_tau_from_H_indices(GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")
    .Define("Gen_tau_idx",        f"Get_taus(Gen_taus_idx, GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")
    
    .Define("tau_channel",        f"Get_tau_channel(Gen_taus_idx, GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")
    .Define("tau_first_pdg",      f"GenPart_pdgId[Gen_tau_idx[0]]")
    .Define("tau_2nd_pdg",        f"GenPart_pdgId[Gen_tau_idx[1]]")
    
    .Define("Gen_tau_p4",         f"Get_visible_tau_p4s(Gen_tau_idx, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
    .Define("Gen_tau_pt",         f"Get_p4_pt(Gen_tau_p4)")
    .Define("Gen_tau_eta",        f"Get_p4_eta(Gen_tau_p4)")
    .Define("Gen_tau_phi",        f"Get_p4_phi(Gen_tau_p4)")
    .Define("Gen_tau_mass",       f"Get_p4_mass(Gen_tau_p4)")

    .Define("Gen_tau_all_p4",         f"Get_all_taus_p4(Gen_tau_idx, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
    .Define("Gen_tau_all_pt",         f"Get_p4_pt(Gen_tau_all_p4)")
    .Define("Gen_tau_all_eta",        f"Get_p4_eta(Gen_tau_all_p4)")
    .Define("Gen_tau_all_phi",        f"Get_p4_phi(Gen_tau_all_p4)")
    .Define("Gen_tau_all_mass",       f"Get_p4_mass(Gen_tau_all_p4)")

    .Define("Gen_tau_leading_pt",         f"Gen_tau_p4[0].pt()")
    .Define("Gen_tau_leading_eta",        f"Gen_tau_p4[0].eta()")
    .Define("Gen_tau_leading_phi",        f"Gen_tau_p4[0].phi()")
    .Define("Gen_tau_leading_mass",       f"Gen_tau_p4[0].M()")

    .Define("Gen_tau_subleading_pt",      f"Gen_tau_p4[1].pt()")
    .Define("Gen_tau_subleading_eta",     f"Gen_tau_p4[1].eta()")
    .Define("Gen_tau_subleading_phi",     f"Gen_tau_p4[1].phi()")
    .Define("Gen_tau_subleading_mass",    f"Gen_tau_p4[1].M()")

    .Define("Gen_Muon_idx",       f"Get_lepton_idx(13, GenPart_pdgId, Gen_taus_idx, GenPart_statusFlags, GenPart_genPartIdxMother)")
    .Define("nGen_Muon",          f"Gen_Muon_idx.size()")
    
    .Define("Gen_muon_pt",        f"Get_variable(Gen_Muon_idx, GenPart_pt)")
    .Define("Gen_muon_eta",       f"Get_variable(Gen_Muon_idx, GenPart_eta)")
    .Define("Gen_muon_phi",       f"Get_variable(Gen_Muon_idx, GenPart_phi)")
    .Define("Gen_muon_mass",      f"Get_variable(Gen_Muon_idx, GenPart_mass)")

    .Define("Gen_Electron_idx",   f"Get_lepton_idx(11, GenPart_pdgId, Gen_taus_idx, GenPart_statusFlags, GenPart_genPartIdxMother)")
    .Define("nGen_Electron",      f"Gen_Electron_idx.size()")

    .Define("Gen_electron_pt",    f"Get_variable(Gen_Electron_idx, GenPart_pt)")
    .Define("Gen_electron_eta",   f"Get_variable(Gen_Electron_idx, GenPart_eta)")
    .Define("Gen_electron_phi",   f"Get_variable(Gen_Electron_idx, GenPart_phi)")
    .Define("Gen_electron_mass",  f"Get_variable(Gen_Electron_idx, GenPart_mass)")
   
    .Define("Gen_H_tautau_mass",  f"Get_Gen_H_tautau_mass(Gen_tau_all_p4)")

    .Define("Gen_HH_mass",        f"Get_gen_mHH(Gen_tau_p4, Gen_b_p4)") 
    
   # ───────────────────────────────────────────────────────────────────── Reco ────────────────────────────────────────────────────────────────────────────────

    .Define("Matched_tau_idx",              f"deltaR_matching(Gen_tau_idx, GenPart_pdgId, Gen_tau_pt, Gen_tau_eta, Gen_tau_phi, hltHpsPFTau_pt, hltHpsPFTau_eta, hltHpsPFTau_phi, hltElectron_pt, hltElectron_eta, hltElectron_phi, hltMuon_pt, hltMuon_eta, hltMuon_phi, hltHpsPFTau_deepTauVSjet, false)")
    .Define("Matched_tagged_tau_idx",       f"deltaR_matching(Gen_tau_idx, GenPart_pdgId, Gen_tau_pt, Gen_tau_eta, Gen_tau_phi, hltHpsPFTau_pt, hltHpsPFTau_eta, hltHpsPFTau_phi, hltElectron_pt, hltElectron_eta, hltElectron_phi, hltMuon_pt, hltMuon_eta, hltMuon_phi, hltHpsPFTau_deepTauVSjet, true)")

    .Define("Reco_tau_pt",           "Get_variable(Matched_tau_idx, hltHpsPFTau_pt)")
    .Define("Reco_tau_eta",          "Get_variable(Matched_tau_idx, hltHpsPFTau_eta)")
    .Define("Reco_tau_phi",          "Get_variable(Matched_tau_idx, hltHpsPFTau_phi)")
    .Define("Reco_tau_mass",         "Get_variable(Matched_tau_idx, hltHpsPFTau_mass)")
    .Define("Reco_tau_DR",           "Get_dR_tau(Gen_tau_p4, Matched_tau_idx, hltHpsPFTau_eta, hltHpsPFTau_phi)")

    .Define("Reco_leading_tau_pt",    "Reco_tau_pt[0]")
    .Define("Reco_leading_tau_eta",   "Reco_tau_eta[0]")
    .Define("Reco_leading_tau_phi",   "Reco_tau_phi[0]")   
    .Define("Reco_leading_tau_mass",  "Reco_tau_mass[0]")   
    .Define("Reco_leading_tau_DR",    "Reco_tau_DR[0]") 

    .Define("Reco_subleading_tau_pt",    "Reco_tau_pt[1]")
    .Define("Reco_subleading_tau_eta",   "Reco_tau_eta[1]")
    .Define("Reco_subleading_tau_phi",   "Reco_tau_phi[1]")   
    .Define("Reco_subleading_tau_mass",  "Reco_tau_mass[1]")   
    .Define("Reco_subleading_tau_DR",    "Reco_tau_DR[1]") 

    .Define("Reco_tagged_tau_pt",           "Get_variable(Matched_tagged_tau_idx, hltHpsPFTau_pt)")
    .Define("Reco_tagged_tau_eta",          "Get_variable(Matched_tagged_tau_idx, hltHpsPFTau_eta)")
    .Define("Reco_tagged_tau_phi",          "Get_variable(Matched_tagged_tau_idx, hltHpsPFTau_phi)")
    .Define("Reco_tagged_tau_mass",         "Get_variable(Matched_tagged_tau_idx, hltHpsPFTau_mass)")
    .Define("Reco_tagged_tau_DR",           "Get_dR_tau(Gen_tau_p4, Matched_tagged_tau_idx, hltHpsPFTau_eta, hltHpsPFTau_phi)")

    .Define("Reco_leading_tagged_tau_pt",    "Reco_tagged_tau_pt[0]")
    .Define("Reco_leading_tagged_tau_eta",   "Reco_tagged_tau_eta[0]")
    .Define("Reco_leading_tagged_tau_phi",   "Reco_tagged_tau_phi[0]")   
    .Define("Reco_leading_tagged_tau_mass",  "Reco_tagged_tau_mass[0]")   
    .Define("Reco_leading_tagged_tau_DR",    "Reco_tagged_tau_DR[0]") 

    .Define("Reco_subleading_tagged_tau_pt",    "Reco_tagged_tau_pt[1]")
    .Define("Reco_subleading_tagged_tau_eta",   "Reco_tagged_tau_eta[1]")
    .Define("Reco_subleading_tagged_tau_phi",   "Reco_tagged_tau_phi[1]")   
    .Define("Reco_subleading_tagged_tau_mass",  "Reco_tagged_tau_mass[1]")   
    .Define("Reco_subleading_tagged_tau_DR",    "Reco_tagged_tau_DR[1]")    

    .Define("Matched_jet_idx",       "deltaR_matching_jets(Gen_bJets_idx, GenJet_eta, GenJet_phi, hltAK4PuppiJet_eta, hltAK4PuppiJet_phi, hltAK4PuppiJet_DeepFlavour_prob_b, hltAK4PuppiJet_DeepFlavour_prob_bb, hltAK4PuppiJet_DeepFlavour_prob_c, hltAK4PuppiJet_DeepFlavour_prob_g, hltAK4PuppiJet_DeepFlavour_prob_lepb, hltAK4PuppiJet_DeepFlavour_prob_uds, false)")
    .Define("Matched_bjet_idx",      "deltaR_matching_jets(Gen_bJets_idx, GenJet_eta, GenJet_phi, hltAK4PuppiJet_eta, hltAK4PuppiJet_phi, hltAK4PuppiJet_DeepFlavour_prob_b, hltAK4PuppiJet_DeepFlavour_prob_bb, hltAK4PuppiJet_DeepFlavour_prob_c, hltAK4PuppiJet_DeepFlavour_prob_g, hltAK4PuppiJet_DeepFlavour_prob_lepb, hltAK4PuppiJet_DeepFlavour_prob_uds, true)")
    
    .Define("Matched_jet_bscore",    "Get_bscore(Gen_bJets_idx, hltAK4PuppiJet_DeepFlavour_prob_b, hltAK4PuppiJet_DeepFlavour_prob_bb, hltAK4PuppiJet_DeepFlavour_prob_c, hltAK4PuppiJet_DeepFlavour_prob_g, hltAK4PuppiJet_DeepFlavour_prob_lepb, hltAK4PuppiJet_DeepFlavour_prob_uds)")
    
    .Define("Reco_jet_pT",           "Get_variable(Matched_jet_idx, hltAK4PuppiJet_pt)")
    .Define("Reco_jet_eta",          "Get_variable(Matched_jet_idx, hltAK4PuppiJet_eta)")
    .Define("Reco_jet_phi",          "Get_variable(Matched_jet_idx, hltAK4PuppiJet_phi)")
    .Define("Reco_jet_mass",         "Get_variable(Matched_jet_idx, hltAK4PuppiJet_mass)")
    .Define("Reco_jet_DR",           "Get_dR_jet(Gen_bJets_idx, Matched_jet_idx, GenJet_eta, GenJet_phi, hltAK4PuppiJet_eta, hltAK4PuppiJet_phi)")

    .Define("Reco_leading_jet_pT",   "Reco_jet_pT[0]")
    .Define("Reco_leading_jet_eta",  "Reco_jet_eta[0]")
    .Define("Reco_leading_jet_phi",  "Reco_jet_phi[0]")
    .Define("Reco_leading_jet_mass", "Reco_jet_mass[0]")
    
    .Define("Reco_subleading_jet_pT",   "Reco_jet_pT[1]")
    .Define("Reco_subleading_jet_eta",  "Reco_jet_eta[1]")
    .Define("Reco_subleading_jet_phi",  "Reco_jet_phi[1]")
    .Define("Reco_subleading_jet_mass", "Reco_jet_mass[1]")


    .Define("Reco_bjet_pT",           "Get_variable(Matched_bjet_idx, hltAK4PuppiJet_pt)")
    .Define("Reco_bjet_eta",          "Get_variable(Matched_bjet_idx, hltAK4PuppiJet_eta)")
    .Define("Reco_bjet_phi",          "Get_variable(Matched_bjet_idx, hltAK4PuppiJet_phi)")
    .Define("Reco_bjet_mass",         "Get_variable(Matched_bjet_idx, hltAK4PuppiJet_mass)")
    .Define("Reco_bjet_DR",           "Get_dR_jet(Gen_bJets_idx, Matched_bjet_idx, GenJet_eta, GenJet_phi, hltAK4PuppiJet_eta, hltAK4PuppiJet_phi)")

    .Define("Reco_leading_bjet_pT",   "Reco_bjet_pT[0]")
    .Define("Reco_leading_bjet_eta",  "Reco_bjet_eta[0]")
    .Define("Reco_leading_bjet_phi",  "Reco_bjet_phi[0]")
    .Define("Reco_leading_bjet_mass", "Reco_bjet_mass[0]")
    
    .Define("Reco_subleading_bjet_pT",   "Reco_bjet_pT[1]")
    .Define("Reco_subleading_bjet_eta",  "Reco_bjet_eta[1]")
    .Define("Reco_subleading_bjet_phi",  "Reco_bjet_phi[1]")
    .Define("Reco_subleading_bjet_mass", "Reco_bjet_mass[1]")


    .Define("Reco_H_tautau_p4",      "Build_Higgs_p4(Matched_tau_idx, hltHpsPFTau_pt, hltHpsPFTau_eta, hltHpsPFTau_phi, hltHpsPFTau_mass)")
    .Define("Reco_H_tautau_mass",    "Get_Higgs_mass(Reco_H_tautau_p4)")

    .Define("Reco_H_bb_p4",          "Build_Higgs_p4(Matched_jet_idx, hltAK4PuppiJet_pt, hltAK4PuppiJet_eta, hltAK4PuppiJet_phi, hltAK4PuppiJet_mass)")
    .Define("Reco_H_bb_mass",        "Get_Higgs_mass(Reco_H_bb_p4)")

    .Define("Reco_HH_mass",          "Get_mHH(Reco_H_tautau_p4, Reco_H_bb_p4)")   

    .Define("Matched_muon_idx",     f"deltaR_matching_lepton(Gen_muon_pt, Gen_muon_eta, Gen_muon_phi, hltMuon_pt, hltMuon_eta, hltMuon_phi)")
    .Define("Matched_electron_idx", f"deltaR_matching_lepton(Gen_electron_pt, Gen_electron_eta, Gen_electron_phi, hltElectron_pt, hltElectron_eta, hltElectron_phi)")
 
    .Define("Reco_electron_pt",      "Get_variable(Matched_electron_idx, hltElectron_pt)")
    .Define("Reco_muon_pt",          "Get_variable(Matched_muon_idx, hltMuon_pt)")

)

#─────────────────────────────────────────────────────────────────── Print column names ────────────────────────────────────────────────────────────────────────────────

# all_cols = [str(c) for c in df.GetColumnNames()]
# L1_cols = [c for c in all_cols if "L1_p" in c]
# print("\nGenPart columns found:")
# for c in sorted(L1_cols):
    # print(f"  {c}  [{df.GetColumnType(c)}]")

df.Display(["tau_channel", "tau_first_pdg", "tau_2nd_pdg"], 50).Print()

#──────────────────────────────────────────────────────────────────── Filter dataframe ────────────────────────────────────────────────────────────────────────────────
channels = {
    0: {"filter": "tau_channel == 0", "suffix": "e", "title": r"$\tau \tau \to \tau_e \tau_h$"},
    1: {"filter": "tau_channel == 1", "suffix": "m", "title": r"$\tau \tau \to \tau_{\mu} \tau_h$"},
    2: {"filter": "tau_channel == 2", "suffix": "h", "title": r"$\tau \tau \to \tau_h \tau_h$"},
}
 
channel = {idx: df.Filter(cfg["filter"]) for idx, cfg in channels.items()}
 
streams = {
"HLT_e":    channel[0].Filter("HLT_Ele32_WPTight_L1Seeded == true && HLT_Ele30_WPTight_L1Seeded_LooseDeepTauPFTauHPS30_eta2p1_CrossL1 == true"),
"HLT_m":    channel[1].Filter("HLT_IsoMu24_FromL1TkMuon == true && HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1 == true"),
"HLT_h":    channel[2].Filter("HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1 == true"),
 
"NGT_e":    channel[0].Filter("DST_PFScouting == true"),
"NGT_m":    channel[1].Filter("DST_PFScouting == true"),
"NGT_h":    channel[2].Filter("DST_PFScouting == true"),
}

streams_reco = {
    "Reco_tau_e": channel[0].Filter("DST_PFScouting == true").Filter("(Matched_tau_idx[0] >= 0 || Matched_tau_idx[1] >= 0) && Matched_electron_idx[0] >= 0"),
    "Reco_tau_m": channel[1].Filter("DST_PFScouting == true").Filter("(Matched_tau_idx[0] >= 0 || Matched_tau_idx[1] >= 0) && Matched_muon_idx[0] >= 0"),
    "Reco_tau_h": channel[2].Filter("DST_PFScouting == true").Filter("Matched_tau_idx[0] >= 0 && Matched_tau_idx[1] >= 0"),
    

    "Reco_jet_e": channel[0].Filter("DST_PFScouting == true").Filter("Matched_jet_idx[0] >= 0 && Matched_jet_idx[1] >= 0"),
    "Reco_jet_m": channel[1].Filter("DST_PFScouting == true").Filter("Matched_jet_idx[0] >= 0 && Matched_jet_idx[1] >= 0"),
    "Reco_jet_h": channel[2].Filter("DST_PFScouting == true").Filter("Matched_jet_idx[0] >= 0 && Matched_jet_idx[1] >= 0"),

    "Reco_both_e": channel[0].Filter("DST_PFScouting == true").Filter("(Matched_tau_idx[0] >= 0 || Matched_tau_idx[1] >= 0) && Matched_electron_idx[0] >= 0 && Matched_jet_idx[0] >= 0 && Matched_jet_idx[1] >= 0"),
    "Reco_both_m": channel[1].Filter("DST_PFScouting == true").Filter("(Matched_tau_idx[0] >= 0 || Matched_tau_idx[1] >= 0) && Matched_muon_idx[0] >= 0 && Matched_jet_idx[0] >= 0 && Matched_jet_idx[1] >= 0"),
    "Reco_both_h": channel[2].Filter("DST_PFScouting == true").Filter("Matched_tau_idx[0] >= 0 && Matched_tau_idx[1] >= 0 && Matched_jet_idx[0] >= 0 && Matched_jet_idx[1] >= 0"),

##Tagging

    "Reco_tagged_tau_e": channel[0].Filter("DST_PFScouting == true").Filter("(Matched_tagged_tau_idx[0] >= 0 || Matched_tagged_tau_idx[1] >= 0) && Matched_electron_idx[0] >= 0"),
    "Reco_tagged_tau_m": channel[1].Filter("DST_PFScouting == true").Filter("(Matched_tagged_tau_idx[0] >= 0 || Matched_tagged_tau_idx[1] >= 0) && Matched_muon_idx[0] >= 0"),
    "Reco_tagged_tau_h": channel[2].Filter("DST_PFScouting == true").Filter("Matched_tagged_tau_idx[0] >= 0 && Matched_tagged_tau_idx[1] >= 0"),
    

    "Reco_bjet_e": channel[0].Filter("DST_PFScouting == true").Filter("Matched_bjet_idx[0] >= 0 && Matched_bjet_idx[1] >= 0"),
    "Reco_bjet_m": channel[1].Filter("DST_PFScouting == true").Filter("Matched_bjet_idx[0] >= 0 && Matched_bjet_idx[1] >= 0"),
    "Reco_bjet_h": channel[2].Filter("DST_PFScouting == true").Filter("Matched_bjet_idx[0] >= 0 && Matched_bjet_idx[1] >= 0"),

    "Reco_both_tagged_e": channel[0].Filter("DST_PFScouting == true").Filter("(Matched_tagged_tau_idx[0] >= 0 || Matched_tagged_tau_idx[1] >= 0) && Matched_electron_idx[0] >= 0 && Matched_bjet_idx[0] >= 0 && Matched_bjet_idx[1] >= 0"),
    "Reco_both_tagged_m": channel[1].Filter("DST_PFScouting == true").Filter("(Matched_tagged_tau_idx[0] >= 0 || Matched_tagged_tau_idx[1] >= 0) && Matched_muon_idx[0] >= 0 && Matched_bjet_idx[0] >= 0 && Matched_bjet_idx[1] >= 0"),
    "Reco_both_tagged_h": channel[2].Filter("DST_PFScouting == true").Filter("Matched_tagged_tau_idx[0] >= 0 && Matched_tagged_tau_idx[1] >= 0 && Matched_bjet_idx[0] >= 0 && Matched_bjet_idx[1] >= 0"),

}

# ────────────────────────────────────────────────────────────── Define histograms  ────────────────────────────────────────────────────────────────────────────────    
 
var_configs = {
    # "Matched_jet_bscore": {"branch": "Matched_jet_bscore", "xlabel": r"b-score", "bins": 30, "range": (0, 1)},
    
    # "Reco_tau_DR": {"branch": "Reco_tau_DR", "xlabel": r"$\Delta$ R", "bins": 30, "range": (0, 0.5)},
    # "Reco_jet_DR": {"branch": "Reco_jet_DR", "xlabel": r"$\Delta$ R", "bins": 30, "range": (0, 0.5)},
    
    # "Gen_tau_pt":      {"branch": "Gen_tau_pt",      "xlabel": r"Gen $p_{T}(\tau's)$ [GeV]", "bins": 30, "range": (0, 300)},
    # "Gen_tau_all_pt":  {"branch": "Gen_tau_all_pt",  "xlabel": r"Gen $p_{T}(\tau's)$ [GeV]", "bins": 30, "range": (0, 300)},

    "Gen_tau_leading_pt":     {"branch": "Gen_tau_leading_pt",     "xlabel": r"Gen $p_{T}(\mu)$ [GeV]",    "bins": 30, "range": (0, 300)},
    "Gen_tau_subleading_pt":  {"branch": "Gen_tau_subleading_pt",  "xlabel": r"Gen $p_{T}(\tau)$ [GeV]", "bins": 30, "range": (0, 300)},
    
    "Gen_tau_leading_eta":    {"branch": "Gen_tau_leading_eta",    "xlabel": r"Gen $\eta(e)$ Leading",    "bins": 30, "range": (-3, 3)},
    "Gen_tau_subleading_eta": {"branch": "Gen_tau_subleading_eta", "xlabel": r"Gen $\eta(\tau)$", "bins": 30, "range": (-3, 3)},
    
    # "bscore":          {"branch": "Matched_jet_bscore", "xlabel": r"$b-score$",    "bins": 30, "range": (0, 1)},

    # "nGen_Muon":       {"branch": "nGen_Muon",     "xlabel": r"n Gen $\mu$",       "bins": 4, "range": (-0.5, 3.5)},
    # "nGen_Electron":   {"branch": "nGen_Electron", "xlabel": r"n Gen $e$",         "bins": 4, "range": (-0.5, 3.5)},
    
    "Gen_b_leading_jet_pt":     {"branch": "Gen_b_leading_jet_pt",     "xlabel": r"Gen $p_{T}$ Leading jet [GeV]",   "bins": 30, "range": (0, 300)},
    "Gen_b_subleading_jet_pt":  {"branch": "Gen_b_subleading_jet_pt",  "xlabel": r"Gen $p_{T}$ Subleading jet [GeV]","bins": 30, "range": (0, 300)},
    
    "Gen_b_leading_jet_eta":    {"branch": "Gen_b_leading_jet_eta",    "xlabel": r"Gen $\eta$ Leading jet",          "bins": 30, "range": (-3, 3)},
    "Gen_b_subleading_jet_eta": {"branch": "Gen_b_subleading_jet_eta", "xlabel": r"Gen $\eta$ Subleading jet",       "bins": 30, "range": (-3, 3)},
    
    # "Gen_electron_pt":  {"branch": "Gen_electron_pt",  "xlabel": r"Gen $p_{T}(e)$ [GeV]",   "bins": 30, "range": (0, 200)},
    # "Gen_muon_pt":      {"branch": "Gen_muon_pt",      "xlabel": r"Gen $p_{T}(\mu)$ [GeV]", "bins": 30, "range": (0, 200)},

    # "Gen_H_tautau_mass":       {"branch": "Gen_H_tautau_mass",     "xlabel": r"Gen $m_H$ [GeV]",  "bins": 30, "range": (110, 140)},



    # "Reco_leading_tau_pt":     {"branch": "Reco_leading_tau_pt",     "xlabel": r"Reco $p_{T}(\tau)$ Leading [GeV]",    "bins": 30, "range": (0, 300)},
    # "Reco_subleading_tau_pt":  {"branch": "Reco_subleading_tau_pt",  "xlabel": r"Reco $p_{T}(\tau)$ Subleading [GeV]", "bins": 30, "range": (0, 300)},

    # "Reco_leading_jet_pT":     {"branch": "Reco_leading_jet_pT",    "xlabel": r"Reco $p_{T}$ Leading jet [GeV]",    "bins": 30, "range": (0, 300)},
    # "Reco_subleading_jet_pT":  {"branch": "Reco_subleading_jet_pT", "xlabel": r"Reco $p_{T}$ Subleading jet [GeV]", "bins": 30, "range": (0, 300)},
    
    # "Reco_leading_jet_eta":    {"branch": "Reco_leading_jet_eta",    "xlabel": r"Reco $\eta$ Leading jet",    "bins": 30, "range": (-3, 3)},
    # "Reco_subleading_jet_eta": {"branch": "Reco_subleading_jet_eta", "xlabel": r"Reco $\eta$ Subleading jet", "bins": 30, "range": (-3, 3)},
    
    # "Reco_electron_pt": {"branch": "Reco_electron_pt", "xlabel": r"Reco $p_{T}(e)$ [GeV]",    "bins": 30, "range": (0, 200)},
    # "Reco_muon_pt":     {"branch": "Reco_muon_pt",     "xlabel": r"Reco $p_{T}(\mu)$ [GeV]",  "bins": 30, "range": (0, 200)},

}
    # "all": {"label": r"All events",          "color": "mediumvioletred",  "histtype": "step", "hatch": ""},
    # "all": {"label": r"Gen-Reco $\tau$/e",           "color": "#C39047",  "histtype": "fill", "hatch": ""},

selection_style = {
    "all": {"label": r"All events",          "color": "mediumvioletred",  "histtype": "step", "hatch": ""},
    "HLT": {"label": "Events that pass HLT", "color": "purple",           "histtype": "fill", "hatch": ""},
    "NGT": {"label": "Events that pass NGT", "color": "darkorchid",       "histtype": "step", "hatch": "//"},
}


branches = [vcfg["branch"] for vcfg in var_configs.values()]

data = {}
for ch, cfg in channels.items():
    suf = cfg["suffix"]
 
    all_arrs = get_leg_arrays(channel[ch],           branches)
    hlt_arrs = get_leg_arrays(streams[f"HLT_{suf}"], branches)
    ngt_arrs = get_leg_arrays(streams[f"NGT_{suf}"], branches)
    
    reco_tau_arrs  = get_leg_arrays(streams_reco[f"Reco_tau_{suf}"], branches)
    reco_jet_arrs  = get_leg_arrays(streams_reco[f"Reco_jet_{suf}"], branches)
    reco_both_arrs = get_leg_arrays(streams_reco[f"Reco_both_{suf}"], branches)

    reco_tagged_tau_arrs  = get_leg_arrays(streams_reco[f"Reco_tagged_tau_{suf}"], branches)
    reco_bjet_arrs        = get_leg_arrays(streams_reco[f"Reco_bjet_{suf}"], branches)
    reco_both_tagged_arrs = get_leg_arrays(streams_reco[f"Reco_both_tagged_{suf}"], branches)

    for var, vcfg in var_configs.items():
        data[f"{var}_all_{suf}"] = all_arrs[vcfg["branch"]]
        data[f"{var}_HLT_{suf}"] = hlt_arrs[vcfg["branch"]]
        data[f"{var}_NGT_{suf}"] = ngt_arrs[vcfg["branch"]]

        data[f"{var}_Reco_tau_{suf}"]  = reco_tau_arrs[vcfg["branch"]]
        data[f"{var}_Reco_jet_{suf}"]  = reco_jet_arrs[vcfg["branch"]]
        data[f"{var}_Reco_both_{suf}"] = reco_both_arrs[vcfg["branch"]]

        data[f"{var}_Reco_tagged_tau_{suf}"]  = reco_tagged_tau_arrs[vcfg["branch"]]
        data[f"{var}_Reco_bjet_{suf}"]  = reco_bjet_arrs[vcfg["branch"]]
        data[f"{var}_Reco_both_tagged_{suf}"] = reco_both_tagged_arrs[vcfg["branch"]]

 
#──────────────────────────────────────────────────────────────────────── Plotting ────────────────────────────────────────────────────────────────────────────────

def plot_kinematic(var, ch):
    path = "/eos/user/s/sbenabde/CERN_Summer_student/HH_bbtautau/plots/"
    os.makedirs(path, exist_ok=True)
 
    fig, ax = plt.subplots(figsize=(10, 8))
    suf   = channels[ch]["suffix"]
    vcfg  = var_configs[var]
    bins, rng = vcfg["bins"], vcfg["range"]
 
    for sel, cfg in selection_style.items():
        arr = data[f"{var}_{sel}_{suf}"]
        c, e = np.histogram(arr, bins=bins, range=rng)
        hep.histplot(c, e, ax=ax, label=cfg["label"], color=cfg["color"], histtype=cfg["histtype"], hatch=cfg["hatch"])
 
    ax.set_xlabel(vcfg["xlabel"])
    ax.set_ylabel("Events")
    ax.legend(loc="upper right", fontsize=20)
    ax.text(0.125, -0.08, channels[ch]["title"], transform=ax.transAxes, fontsize=20, ha="right", va="top")
    hep.cms.label("Preliminary", data=False, ax=ax, com=14)
 
    outname = f"plots/kinematics/{var}_{suf}.png"
    plt.tight_layout()
    plt.savefig(outname, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved {outname}")

        
def plot_efficiency(var, ch, save_as, scale):
    path = "/eos/user/s/sbenabde/CERN_Summer_student/HH_bbtautau/plots/Efficiency"
    os.makedirs(path, exist_ok=True)
 
    fig, ax = plt.subplots(figsize=(10, 8))
    suf   = channels[ch]["suffix"]
    vcfg  = var_configs[var]
    bins, rng = vcfg["bins"], vcfg["range"]
 
    hist_data = {}
    for sel, cfg in selection_style.items():
        arr = data[f"{var}_{sel}_{suf}"]
        counts, edges = np.histogram(arr, bins=bins, range=rng, density=False)
        hist_data[sel] = counts
        hep.histplot(counts, edges, ax=ax, label=cfg["label"], color=cfg["color"], histtype=cfg["histtype"], hatch=cfg["hatch"])
 
    ax.set_xlabel(vcfg["xlabel"], fontsize=24)
    ax.set_ylabel("Events", fontsize=24)
    hep.cms.label("Preliminary", data=False, ax=ax, com=14)
    ax.text(0.125, -0.08, channels[ch]["title"], transform=ax.transAxes, fontsize=20, ha="right", va="top")
 
    centers     = 0.5 * (edges[:-1] + edges[1:])
    bin_widths  = np.diff(edges)
    den         = hist_data["all"]
  
    eff_HLT, err_HLT = eff_and_err(hist_data["HLT"], den)
    eff_NGT, err_NGT = eff_and_err(hist_data["NGT"], den)
 
    eff_color = "olivedrab"
    ax2 = ax.twinx()
    ax2.errorbar(centers, eff_HLT, xerr=0.5 * bin_widths, yerr=err_HLT, color='olivedrab', fmt='*',   label="HLT efficiency")
    ax2.errorbar(centers, eff_NGT, xerr=0.5 * bin_widths, yerr=err_NGT, color='yellowgreen', fmt='*', label="NGT efficiency")
    ax2.set_ylabel("Efficiency", color=eff_color)
    ax2.tick_params(axis='y', labelcolor=eff_color)
    ax2.set_ylim(0, 1)
    plt.tight_layout()
    ax2.legend(loc="upper left", fontsize=20)
    ax.legend(loc="upper right", fontsize=20)
    if scale == 'log':
        ax.set_yscale('log')
 
    outname = f"{path}/{var}_{suf}_{scale}.{save_as}"
    fig.savefig(outname)
    plt.close(fig)
    print(f"Saved {outname}")
 

def plot_reco_efficiency(var, ch, save_as, scale, No_tag = False, Tau=False, Jet=False, Both=False):
    fig, ax = plt.subplots(figsize=(10, 8))
    suf   = channels[ch]["suffix"]
    vcfg  = var_configs[var]
    bins, rng = vcfg["bins"], vcfg["range"]

    hist_data = {}
    arr_all  = data[f"{var}_all_{suf}"]
    
    if Tau==True:
        arr_reco = data[f"{var}_Reco_tau_{suf}"]
        arr_reco_tagged = data[f"{var}_Reco_tagged_tau_{suf}"]
        lb   = r"Matched+$\tau$-tagged"
        path = "/eos/user/s/sbenabde/CERN_Summer_student/HH_bbtautau/plots/Reco_Efficiency_noBTagging_TEST"

    if Jet==True:
        arr_reco = data[f"{var}_Reco_jet_{suf}"]
        arr_reco_tagged = data[f"{var}_Reco_bjet_{suf}"]
        lb   = r"Matched+b-tagged"
        path = "/eos/user/s/sbenabde/CERN_Summer_student/HH_bbtautau/plots/Reco_Efficiency_noTauTagging_TEST"
    
    if Both==True:
        arr_reco = data[f"{var}_Reco_both_{suf}"]
        arr_reco_tagged = data[f"{var}_Reco_both_tagged_{suf}"]
        lb   = r"Matched+$\tau$&b-tagged"
        path = "/eos/user/s/sbenabde/CERN_Summer_student/HH_bbtautau/plots/Reco_Efficiency_withTagging_TEST"
   
    if No_tag==True: 
        arr_reco        = data[f"{var}_NGT_{suf}"]
        arr_reco_tagged = data[f"{var}_NGT_{suf}"]
        lb   = r"Matched"
        path = "/eos/user/s/sbenabde/CERN_Summer_student/HH_bbtautau/plots/Reco_Efficiency_noTagging_TEST"
    
    os.makedirs(path, exist_ok=True)
    
    counts, edges = np.histogram(arr_all, bins=bins, range=rng, density=False)
    counts2, edges2 = np.histogram(arr_reco, bins=bins, range=rng, density=False)
    counts3, edges3 = np.histogram(arr_reco_tagged, bins=bins, range=rng, density=False)
    hist_data['all']  = counts
    hist_data['Reco'] = counts2
    hist_data['Reco_tagged'] = counts3
    hep.histplot(counts, edges, ax=ax,   label="Pass NGT",          color="mediumblue",         histtype="step", hatch= ""),
    hep.histplot(counts2, edges2, ax=ax, label="Matched",           color="mediumslateblue",    histtype="step", hatch= "//"),
    hep.histplot(counts3, edges3, ax=ax, label="Matched+tagged",    color="steelblue",          histtype="fill", hatch= ""),

    ax.set_xlabel(vcfg["xlabel"])
    ax.set_ylabel("Events")
    hep.cms.label("Preliminary", data=False, ax=ax, com=14)
    ax.text(0.125, -0.08, channels[ch]["title"], transform=ax.transAxes, fontsize=20, ha="right", va="top")

    centers     = 0.5 * (edges[:-1] + edges[1:])
    bin_widths  = np.diff(edges)
    den         = hist_data["all"]

    eff_NGT, err_NGT = eff_and_err(hist_data["Reco"], den)
    eff_NGT2, err_NGT2 = eff_and_err(hist_data["Reco_tagged"], den)

    eff_color = "deeppink"
    ax2 = ax.twinx()
    #ax2.errorbar(centers, eff_NGT, xerr=0.5 * bin_widths,  yerr=err_NGT,  color=eff_color, fmt='*', label="Reco efficiency")
    ax2.errorbar(centers, eff_NGT2, xerr=0.5 * bin_widths, yerr=err_NGT2, color=eff_color, fmt='*', label="Reco efficiency (tagged)")
    ax2.set_ylabel("Efficiency", color=eff_color)
    ax2.tick_params(axis='y',    labelcolor=eff_color)
    ax2.set_ylim(0, 1)
    ax2.legend(loc="upper left", fontsize=20)
    ax.legend(loc="upper right", fontsize=20)
    plt.tight_layout()

    if scale == 'log':
        ax.set_yscale('log')

    outname = f"{path}/{var}_{suf}_{scale}.{save_as}"
    fig.savefig(outname)
    plt.close(fig)
    print(f"Saved /{outname}")


# for var in var_configs:
#     for ch in channels:
#         plot_kinematic(var, ch)

# for var in var_configs:
#     for ch in channels:
#         plot_efficiency(var, ch, save_as='png', scale='lin')

for ch in channels:
    plot_reco_efficiency("Gen_b_leading_jet_pt",     ch, save_as='png', scale='lin', Tau=True)
    plot_reco_efficiency("Gen_b_subleading_jet_pt",  ch, save_as='png', scale='lin', Tau=True)
    plot_reco_efficiency("Gen_b_leading_jet_eta",    ch, save_as='png', scale='lin', Jet=True)
    plot_reco_efficiency("Gen_b_subleading_jet_eta", ch, save_as='png', scale='lin', Jet=True)

    plot_reco_efficiency("Gen_tau_leading_pt",      ch, save_as='png', scale='lin',   Tau=True)
    plot_reco_efficiency("Gen_tau_subleading_pt",   ch, save_as='png', scale='lin',   Tau=True)

    plot_reco_efficiency("Gen_tau_leading_pt",      ch, save_as='png', scale='lin', Jet=True)
    plot_reco_efficiency("Gen_tau_subleading_pt",   ch, save_as='png', scale='lin', Jet=True)
    plot_reco_efficiency("Gen_b_leading_jet_pt",    ch, save_as='png', scale='lin', Jet=True)
    plot_reco_efficiency("Gen_b_subleading_jet_pt", ch, save_as='png', scale='lin', Jet=True)
    plot_reco_efficiency("Gen_tau_leading_eta",     ch, save_as='png', scale='lin', Jet=True)
    plot_reco_efficiency("Gen_tau_subleading_eta",  ch, save_as='png', scale='lin', Jet=True)




#     # data = {

#     #     "gen_b_massH_NGT":  streams.get("NGT_h").AsNumpy(["Gen_H_bb_mass"])["Gen_H_bb_mass"],  
#     #     "gen_bJet_massH_NGT":  streams.get("NGT_h").AsNumpy(["Gen_H_bbjet_mass"])["Gen_H_bbjet_mass"],  

#     #     "gen_mH_NGT_e":  streams.get("NGT_e").AsNumpy(["Gen_H_tautau_mass"])["Gen_H_tautau_mass"],  
#     #     "gen_mH_NGT_m":  streams.get("NGT_m").AsNumpy(["Gen_H_tautau_mass"])["Gen_H_tautau_mass"],  
#     #     "gen_mH_NGT_h":  streams.get("NGT_h").AsNumpy(["Gen_H_tautau_mass"])["Gen_H_tautau_mass"], 

#     #     "mH_NGT_e":  streams.get("NGT_e").AsNumpy(["Reco_H_tautau_mass"])["Reco_H_tautau_mass"],  
#     #     "mH_NGT_m":  streams.get("NGT_m").AsNumpy(["Reco_H_tautau_mass"])["Reco_H_tautau_mass"],  
#     #     "mH_NGT_h":  streams.get("NGT_h").AsNumpy(["Reco_H_tautau_mass"])["Reco_H_tautau_mass"],  

#     #     "mHb_all":  channel.get(0).AsNumpy(["Reco_H_bb_mass"])["Reco_H_bb_mass"],
#     #     "mHb_HLT":  streams.get("HLT_m").AsNumpy(["Reco_H_bb_mass"])["Reco_H_bb_mass"],
#     #     "mHb_NGT":  streams.get("NGT_m").AsNumpy(["Reco_H_bb_mass"])["Reco_H_bb_mass"],    

#     #     "mHH_all":  channel.get(0).AsNumpy(["Reco_HH_mass"])["Reco_HH_mass"],
#     #     "mHH_HLT":  streams.get("HLT_m").AsNumpy(["Reco_HH_mass"])["Reco_HH_mass"],
#     #     "mHH_NGT":  streams.get("NGT_m").AsNumpy(["Reco_HH_mass"])["Reco_HH_mass"], 
        
#     #     "gen_mH_NGT":  df_NGT.AsNumpy(["Gen_H_tautau_mass"])["Gen_H_tautau_mass"],  
#     #     "mH_NGT":      df_NGT.AsNumpy(["Reco_H_tautau_mass"])["Reco_H_tautau_mass"],  

#     #     "gen_mHb_NGT":  df_NGT.AsNumpy(["Gen_H_bb_mass"])["Gen_H_bb_mass"],  
#     #     "mHb_NGT":      df_NGT.AsNumpy(["Reco_H_bb_mass"])["Reco_H_bb_mass"], 

#     #     "gen_mHH_NGT":  df_NGT.AsNumpy(["Gen_HH_mass"])["Gen_HH_mass"],  
#     #     "mHH_NGT":      df_NGT.AsNumpy(["Reco_HH_mass"])["Reco_HH_mass"],  

#     # }


#     # plot_configs = {       

#     #     "gen_b_massH_NGT":         {"label": r"$H \to bb$ from gen b quarks",   "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events that pass NGT Scouting", "bins": 30, "range": (0, 180), "color": "purple", "histtype": "fill", "hatch": ""},
#     #     "gen_bJet_massH_NGT":      {"label": r"$H \to bb$ from gen b jets",     "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events that pass NGT Scouting", "bins": 30, "range": (0, 180), "color": "mediumvioletred", "histtype": "step", "hatch": "/"},

#     #     "mH_NGT_e":      {"label": r"$\tau \tau \to \tau_e \tau_h$",   "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 150), "color": "purple", "histtype": "step", "hatch": ""},
#     #     "mH_NGT_m":      {"label": r"$\tau \tau \to \tau_\mu \tau_h$",   "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 150), "color": "mediumvioletred", "histtype": "step", "hatch": ""},
#     #     "mH_NGT_h":      {"label": r"$\tau \tau \to \tau_h \tau_h$",   "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 150), "color": "orchid", "histtype": "step", "hatch": ""},

#     #     "gen_mH_NGT_e":      {"label": r"GEN $\tau \tau \to \tau_e \tau_h$",   "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 150), "color": "deepskyblue", "histtype": "fill", "hatch": ""},
#     #     "gen_mH_NGT_m":      {"label": r"GEN $\tau \tau \to \tau_\mu \tau_h$", "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 150), "color": "darkcyan", "histtype": "fill", "hatch": ""},
#     #     "gen_mH_NGT_h":      {"label": r"GEN $\tau \tau \to \tau_h \tau_h$",   "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 150), "color": "royalblue", "histtype": "fill", "hatch": ""},

#     #     "gen_mH_NGT":  {"label": r"GEN $H \to \tau \tau$",   "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 150), "color": "purple", "histtype": "fill", "hatch": ""},
#     #     "mH_NGT":      {"label": r"RECO $H \to \tau \tau$",  "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 150), "color": "mediumvioletred", "histtype": "fill", "hatch": ""},

#     #     "gen_mHb_NGT":  {"label": r"Gen $H \to bb$",   "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events that pass NGT Scouting", "bins": 30, "range": (0, 180), "color": "purple", "histtype": "fill", "hatch": ""},
#     #     "mHb_NGT":      {"label": r"Reco $H \to bb$",  "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events that pass NGT Scouting", "bins": 30, "range": (0, 180), "color": "mediumvioletred", "histtype": "fill", "hatch": ""},

#     #     "gen_mHH_NGT":  {"label": r"Gen. $HH \to bb \tau \tau$",   "xlabel": r"$m_{HH}$ [GeV]", "ylabel": "Events that pass NGT Scouting", "bins": 30, "range": (50, 800), "color": "purple", "histtype": "fill", "hatch": ""},
#     #     "mHH_NGT":      {"label": r"Reco. $HH \to bb \tau \tau$",  "xlabel": r"$m_{HH}$ [GeV]", "ylabel": "Events that pass NGT Scouting", "bins": 30, "range": (50, 800), "color": "mediumvioletred", "histtype": "fill", "hatch": ""},

#     #     "mHb_all":      {"label": r"All events", "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 150), "color": "lightcyan", "histtype": "fill", "hatch": ""},
#     #     "mHb_HLT":      {"label": r"Pass HLT",   "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 150), "color": "deepskyblue", "histtype": "step", "hatch": "//"},
#     #     "mHb_NGT":      {"label": r"Pass NGT",   "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 150), "color": "royalblue", "histtype": "step", "hatch": "/"},

#     #     "mHH_all":      {"label": r"All events", "xlabel": r"$m_{HH}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 200), "color": "lightcyan", "histtype": "fill", "hatch": ""},
#     #     "mHH_HLT":      {"label": r"Pass HLT",   "xlabel": r"$m_{HH}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 200), "color": "deepskyblue", "histtype": "step", "hatch": "//"},
#     #     "mHH_NGT":      {"label": r"Pass NGT",   "xlabel": r"$m_{HH}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 200), "color": "royalblue", "histtype": "step", "hatch": "/"},
#     #     }

#     # separate_files = False   # True: one file per histogram. False: overlay all in plot_configs onto one figure.

#     # if separate_files:
#     #     for col, cfg in plot_configs.items():
#     #         fig, ax = plt.subplots(figsize=(10, 8))
#     #         make_panel(ax, col, cfg)
#     #         plt.tight_layout()
#     #         plt.savefig(f"plots/mHH/muon_channel/{col}.png", dpi=300, bbox_inches="tight")
#     #         plt.close()
#     # else:
#     #     fig, ax = plt.subplots(figsize=(10, 8))
#     #     for col, cfg in plot_configs.items():
#     #         make_panel(ax, col, cfg)
#     #     plt.tight_layout()
#     #     plt.yscale('linear')
#     #     plt.savefig("plots/mHH/mH_GEN_b_vs_jet.png", dpi=300, bbox_inches="tight")
#     #     plt.close()

#     # separate_files = False      #True gets 1 file with all the histograms, false gets a separate file for each histogram


#     # data = {      
#     # }

#     # plot_configs = {
#     #     "e_h":  {"label": r"$\tau \tau \to \tau_e \tau_h$",   "bins": 50, "range": (100, 150), "color": "orchid", "histtype": "fill"},
#     #     "mu_h": {"label": r"$\tau \tau \to \tau_{\mu} \tau_h$","bins": 50, "range": (100, 150), "color": "darkorchid", "histtype": "fill"},
#     #     "h_h":  {"label": r"$\tau \tau \to \tau_h \tau_h$",  "bins": 50, "range": (100, 150), "color": "mediumvioletred", "histtype": "fill"},
        
#     #     "h_h_mHH":  {"label": r"All events",  "bins": 50, "range": (100, 150), "color": "orchid", "histtype": "fill"},
#     #     "h_h_mHH_HLT":  {"label": r"Events that pass trigger",  "bins": 50, "range": (100, 150), "color": "mediumvioletred", "histtype": "fill"},

#     #     "hadr_pt_leading":      {"label": r"All events",            "xlabel" : r"$p_T \tau_1$", "ylabel": "Events", "bins": 30, "range": (0, 300), "color": "palevioletred", "histtype": "step", "hatch": ""},
#     #     "hadr_pt_leading_HLT":  {"label": r"Events that pass HLT",  "xlabel" : r"$p_T \tau_1$", "ylabel": "Events", "bins": 30, "range": (0, 300), "color": "purple", "histtype": "fill", "hatch": "//"},
#     #     "hadr_pt_leading_NGT":  {"label": r"Events that pass NGT",  "xlabel" : r"$p_T \tau_1$", "ylabel": "Events", "bins": 30, "range": (0, 300), "color": "magenta", "histtype": "step", "hatch": ""},
        
#     #     "hadr_pt_subleading":      {"label": r"All events",            "xlabel" : r"$p_T \tau_2$", "ylabel": "Events", "bins": 30, "range": (0, 300), "color": "violet", "histtype": "step", "hatch": ""},
#     #     "hadr_pt_subleading_HLT":  {"label": r"Events that pass HLT",  "xlabel" : r"$p_T \tau_2$", "ylabel": "Events", "bins": 30, "range": (0, 300), "color": "orchid", "histtype": "fill", "hatch": ""},
#     #     "hadr_pt_subleading_NGT":  {"label": r"Events that pass NGT",  "xlabel" : r"$p_T \tau_2$", "ylabel": "Events", "bins": 30, "range": (0, 300), "color": "rebeccapurple", "histtype": "step", "hatch": "//"},

#     #     "hadr_eta1":      {"label": r"All events",            "xlabel" : r"$\eta \tau_1$", "ylabel": "Events", "bins": 30, "range": (-3, 3), "color": "olivedrab", "histtype": "step", "hatch": ""},
#     #     "hadr_eta1_HLT":  {"label": r"Events that pass HLT",  "xlabel" : r"$\eta \tau_1$", "ylabel": "Events", "bins": 30, "range": (-3, 3), "color": "yellowgreen", "histtype": "fill", "hatch": ""},
#     #     "hadr_eta1_NGT":  {"label": r"Events that pass NGT",  "xlabel" : r"$\eta \tau_1$", "ylabel": "Events", "bins": 30, "range": (-3, 3), "color": "yellowgreen", "histtype": "step", "hatch": "//"},

#     #     "hadr_pt":      {"label": r"All events",              "xlabel" : r"$p_T$ [GeV]", "ylabel": "Events", "bins": 30, "range": (0, 300), "color": "palevioletred", "histtype": "step", "hatch": ""},

#     #     "hadr_eta2":      {"label": r"All events",            "xlabel" : r"$\eta \tau_2$", "ylabel": "Events", "bins": 30, "range": (-3, 3), "color": "greenyellow", "histtype": "step", "hatch": ""},
#     #     "hadr_eta2_HLT":  {"label": r"Events that pass HLT",  "xlabel" : r"$\eta \tau_2$", "ylabel": "Events", "bins": 30, "range": (-3, 3), "color": "forestgreen",  "histtype": "fill", "hatch": ""},
#     #     "hadr_eta2_NGT":  {"label": r"Events that pass NGT",  "xlabel" : r"$\eta \tau_2$", "ylabel": "Events", "bins": 30, "range": (-3, 3), "color": "forestgreen", "histtype": "step", "hatch": "//"},
#     # }

#     # def make_panel(ax, col, cfg):
#     #     counts, edges = np.histogram(data[col], bins=cfg["bins"], range=cfg["range"])
#     #     hep.histplot(counts, edges, ax=ax, label=cfg["label"], color=cfg["color"], hatch = cfg["hatch"], histtype=cfg["histtype"])
#     #     ax.set_xlabel(cfg["xlabel"])
#     #     ax.set_ylabel(cfg["ylabel"])
#     #     ax.legend(loc="upper right")
#     #     hep.cms.label("Preliminary", data=False, ax=ax, com=14)
#     #     return counts, edges

#     # eff_color = "orchid"
#     # efficiency_counts = []
#     # if separate_files:
#     #     for col, cfg in plot_configs.items():
#     #         fig, ax = plt.subplots(figsize=(10, 8))
#     #         make_panel(ax, col, cfg)
#     #         plt.tight_layout()
#     #         plt.savefig(f"plots/{col}.png", dpi=300, bbox_inches="tight")
#     #         plt.close()
#     # else:
#     #     fig, ax = plt.subplots(figsize=(10, 8))
#     #     for col, cfg in plot_configs.items():
#     #         counts, edges = make_panel(ax, col, cfg)
#     #         efficiency_counts.append(counts)
            
#     #     Efficiency_HLT = [num/den if den != 0 else 0 for num,den in zip(efficiency_counts[1], efficiency_counts[0])]
#     #     Efficiency_NGT = [num/den if den != 0 else 0 for num,den in zip(efficiency_counts[2], efficiency_counts[0])]
#     #     centers = 0.5*(edges[:-1]+edges[1:])
#     #     bin_widths = np.diff(edges)
#     #     y_error_HLT = [np.sqrt(eff * (1 - eff) / den) if den > 0 else 0 for eff, den in zip(Efficiency_HLT, efficiency_counts[0])] 
#     #     y_error_NGT = [np.sqrt(eff * (1 - eff) / den) if den > 0 else 0 for eff, den in zip(Efficiency_NGT, efficiency_counts[0])] 
#     #     ax2 = ax.twinx()

#     #     ax2.errorbar(centers, Efficiency_HLT, xerr=0.5*bin_widths, yerr = y_error_HLT, color=eff_color, fmt = 'o')
#     #     ax2.errorbar(centers, Efficiency_NGT, xerr=0.5*bin_widths, yerr = y_error_NGT, color='purple', fmt = '*')
#     #     ax2.set_ylabel("Efficiency", color=eff_color)
#     #     ax2.tick_params(axis='y', labelcolor=eff_color)
#     #     ax2.set_ylim(0,1)
#     #     hep.cms.label("Preliminary", data=False, ax=ax, com=14)
#     #     plt.tight_layout()
#     #     ax2.legend(loc='upper right')
#     #     plt.savefig("plots/eta/eta2+eff.png", dpi=300, bbox_inches="tight")
#     #     plt.close()

#     #─────────────────────────────────────── Define and plot histograms 2D ────────────────────────────────────────

#     # separate_files_2D = True

#     # data_2D = {

#     # }

#     # plot_configs_2D = {
#     #     "pTLeading_Vs_Subleading_all": {"x": "tau_pt_lead", "y": "tau_pt_sublead", "bins": 50, "range": [[0, 300], [0, 300]], "cmap": "YlGn", "cbar_label": "All events"},
#     #     "pTLeading_Vs_Subleading_HLT": {"x": "tau_pt_lead", "y": "tau_pt_sublead", "bins": 50, "range": [[0, 300], [0, 300]], "cmap": "YlGn", "cbar_label": "Events that pass HLT"},
#     #     "pTLeading_Vs_Subleading_NGT": {"x": "tau_pt_lead", "y": "tau_pt_sublead", "bins": 50, "range": [[0, 300], [0, 300]], "cmap": "YlGn", "cbar_label": "Events that pass NGT"},
        
#     #     "pT_Efficiency_HLT": {"x": "tau_pt_lead", "y": "tau_pt_sublead", "bins": 50, "range": [[0, 300], [0, 300]], "cmap": "YlGn", "cbar_label": "Efficiency HLT"},
#     #     "pT_Efficiency_NGT": {"x": "tau_pt_lead", "y": "tau_pt_sublead", "bins": 50, "range": [[0, 300], [0, 300]], "cmap": "YlGn", "cbar_label": "Efficiency NGT"},
#     # }

#     # bins = plot_configs_2D["eta_Efficiency_HLT"]["bins"]
#     # range_ = plot_configs_2D["eta_Efficiency_HLT"]["range"]

#     # hist_all, xedges, yedges = np.histogram2d(
#     #     data_2D["eta1_Vs_eta2_all"]["eta1"],
#     #     data_2D["eta1_Vs_eta2_all"]["eta2"],
#     #     bins=bins, range=range_
#     # )
#     # hist_HLT, xedges, yedges = np.histogram2d(
#     #     data_2D["eta1_Vs_eta2_HLT"]["eta1"],
#     #     data_2D["eta1_Vs_eta2_HLT"]["eta2"],
#     #     bins=bins, range=range_
#     # )
#     # hist_ngt, xedges, yedges = np.histogram2d(
#     #     data_2D["eta1_Vs_eta2_NGT"]["eta1"],
#     #     data_2D["eta1_Vs_eta2_NGT"]["eta2"],
#     #     bins=bins, range=range_
#     # )
#     # hist_efficiency_HLT = np.divide(hist_HLT, hist_all,
#     #     out=np.zeros_like(hist_HLT, dtype=float),
#     #     where=hist_all > 0
#     # )
#     # hist_efficiency_NGT = np.divide(hist_ngt, hist_all,
#     #     out=np.zeros_like(hist_ngt, dtype=float),
#     #     where=hist_all > 0
#     # )
#     # data_2D["eta_Efficiency_HLT"] = (hist_efficiency_HLT, xedges, yedges) #Appends Efficiency tuple in data_2D dict
#     # data_2D["eta_Efficiency_NGT"] = (hist_efficiency_NGT, xedges, yedges) #Appends Efficiency NGT tuple in data_2D dict

#     # def make_panel_2d(ax, col, cfg):
#     #     if isinstance(data_2D[col], tuple): #For Efficiency histogram, because that one is a tuple in the dictionnary
#     #         h, xedges, yedges = data_2D[col]
#     #     else:                               #The other histograms are stored as dictionnaries in data, not tuples, so different treatement
#     #         h, xedges, yedges = np.histogram2d(
#     #             data_2D[col][cfg["x"]], data_2D[col][cfg["y"]],
#     #             bins=cfg["bins"], range=cfg["range"]
#     #         )

#     #     mesh = ax.pcolormesh(xedges, yedges, h.T, cmap=cfg["cmap"])
#     #     plt.colorbar(mesh, ax=ax, label=cfg["cbar_label"])
#     #     ax.set_xlabel(r"$\eta^{\mathrm{leading}}$")
#     #     ax.set_ylabel(r"$\eta^{\mathrm{subleading}}$")
#     #     hep.cms.label("Preliminary", data=False, ax=ax, com=14)

#     # if separate_files_2D:
#     #     for col, cfg in plot_configs_2D.items():
#     #         fig, ax = plt.subplots(figsize=(10, 8))
#     #         make_panel_2d(ax, col, cfg)
#     #         plt.tight_layout()
#     #         plt.savefig(f"plots/eta/2D_{col}.png", dpi=300, bbox_inches="tight")
#     #         plt.close()
#     # else:
#     #     fig, axes = plt.subplots(1, 4, figsize=(30, 8))
#     #     for ax, (col, cfg) in zip(axes, plot_configs_2D.items()):
#     #         make_panel_2d(ax, col, cfg)
#     #     plt.tight_layout()
#     #     plt.savefig("plots/eta/2D_eta_All+Triggered+EffHLT+EFFNGT.png", dpi=300, bbox_inches="tight")
#     #     plt.close()

#     # import pdb
#     # pdb.set_trace()
#     # all_events = df_tau_h.Histo2D(("h_all_events", "", 50, 0, 300, 50, 0, 300), "tau_pt_lead", "tau_pt_sublead")
#     # all_events_HLT = df_tau_h_HLT.Histo2D(("h_all_events_HLT", "", 50, 0, 300, 50, 0, 300), "tau_pt_lead", "tau_pt_sublead")

#     # h_tot = all_events.GetValue()
#     # h_pas = all_events_HLT.GetValue()

#     # h_eff = h_pas.Clone("h_eff")
#     # h_eff.Reset("ICES")

#     # h_tot.SetTitle(r"Events;$p_T\tau_{1}$;$p_T\tau_{2}$")
#     # h_pas.SetTitle(r"Events;$p_T\tau_{1}$;$p_T\tau_{2}$")
#     # h_eff.SetTitle(r"Efficiency;$p_T\tau_{1}$;$p_T\tau_{2}$")

#     # h_eff.Divide(h_pas, h_tot, 1.0, 1.0, "B")

#     # for h, name in zip([h_tot, h_pas, h_eff], ["2D_pTLeading_Vs_Subleading_total", "2D_pTLeading_Vs_Subleading_HLT", "2D_Efficiency"]):
#     #     xedges = np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1, h.GetNbinsX() + 2)])
#     #     yedges = np.array([h.GetYaxis().GetBinLowEdge(i) for i in range(1, h.GetNbinsY() + 2)])

#     #     #Bin contents
#     #     z = np.zeros((h.GetNbinsY(), h.GetNbinsX()))

#     #     for ix in range(1, h.GetNbinsX() + 1):
#     #         for iy in range(1, h.GetNbinsY() + 1):
#     #             z[iy-1, ix-1] = h.GetBinContent(ix, iy)

#     #     #Plot
#     #     plt.figure(figsize=(10, 8))
#     #     plt.pcolormesh(xedges, yedges, z, cmap='RdPu')
#     #     plt.colorbar(label="Efficiency")
#     #     hep.cms.label("Preliminary", data=False, com=14)
#     #     plt.xlabel(h.GetXaxis().GetTitle())
#     #     plt.ylabel(h.GetYaxis().GetTitle())
#     #     plt.tight_layout()
#     #     plt.savefig(f"/eos/user/s/sbenabde/CERN_Summer_student/plots/{name}")
# #     #─────────────────────────────────────── Save root file output ────────────────────────────────────────

# #     # HH_mass = df.Histo1D(("HH_mass", "HH_mass", 100, -10, 200), "HH_mass")

# #     # out = ROOT.TFile("Signal_test.root", "RECREATE")
# #     # for h in [HH_mass]:
# #     #     h.Write()
# #     # out.Close()
# #     # print("\nSaved Signal_test.root")

