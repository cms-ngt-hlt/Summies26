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
from pathlib import Path
import shutil
import tempfile

REPO_DIR = Path(__file__).resolve().parent

SRC = REPO_DIR / "functions" / "functions.cc"
HDR_DIR = REPO_DIR / "headers"

tmp_dir = Path(tempfile.mkdtemp(prefix="root_aclic_"))
TMP = tmp_dir / "functions.cc"

shutil.copy2(SRC, TMP)
ROOT.gSystem.AddIncludePath(f"-I{HDR_DIR}")
ROOT.gROOT.ProcessLine(f".L {TMP}+")

# ROOT.ROOT.EnableImplicitMT()

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

def overall_efficiency(all, trigger):
    all_c = all.Count().GetValue()
    trigger_c = trigger.Count().GetValue()
    eff = trigger_c/all_c
    return eff

#─────────────────────────────────────── Get files ────────────────────────────────────────
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Argument for Gen_Analysis')
    parser.add_argument('--run-merging', action='store_true', required=False, help='Run merging of NGT and HLT input files.')
    args = parser.parse_args()

    if args.run_merging == True:
        os.makedirs("Dataframes", exist_ok=True)

        base_dir_hlt  = "/eos/user/e/evernazz/Sarah/HHbbtautau/HLTStandard/"
        base_dir_ngt  = "/eos/user/e/evernazz/Sarah/HHbbtautau/NGTScouting/"
        
        files_hlt = sorted(glob.glob(base_dir_hlt + "/step2_*.root"))
        files_ngt = sorted(glob.glob(base_dir_ngt + "/step2_*.root"))

        files_hlt = files_with_events(files_hlt)
        files_ngt = files_with_events(files_ngt)
    
        chain_hlt = make_chain(files_hlt)
        chain_ngt = make_chain(files_ngt)
        
        df_hlt = ROOT.RDataFrame(chain_hlt)
        df_ngt = ROOT.RDataFrame(chain_ngt)

        print(f"Total entries in chain HLT: {df_hlt.Count().GetValue()}")
        print(f"Total entries in chain NGT: {df_ngt.Count().GetValue()}")

        ########################################################################
        # Merge the two dataframes and save (HLT and NGT)
        ########################################################################

        # Keep only non-HLT columns from NGT
        cols_ngt_filtered = [str(c) for c in df_ngt.GetColumnNames() if not str(c).startswith("HLT_")]
        print(f"Keeping {len(cols_ngt_filtered)} NGT columns for merging.")
        # Keep only HLT columns from HLT
        cols_hlt_filtered = [str(c) for c in df_hlt.GetColumnNames() if str(c).startswith("HLT_") and not str(c).endswith("_pHLT")]
        print(f"Keeping {len(cols_hlt_filtered)} HLT columns for merging.")

        chain_ngt.AddFriend(chain_hlt)
        df = ROOT.RDataFrame(chain_ngt)

        print("Saving merged dataframe...")

        os.makedirs("inputs", exist_ok=True)
        snapshot_columns = ROOT.vector(str)()
        for column in cols_ngt_filtered + cols_hlt_filtered:
            snapshot_columns.push_back(column)

        snapshot_options = ROOT.RDF.RSnapshotOptions()
        # snapshot_options.fCompressionAlgorithm = ROOT.kLZ4
        snapshot_options.fCompressionLevel = 4
        df.Snapshot("Events", "/eos/user/s/sbenabde/CERN_Summer_student/Dataframes/df_ngt_hlt.root", snapshot_columns, snapshot_options)
    else:
        df = ROOT.RDataFrame("Events", "/eos/user/s/sbenabde/CERN_Summer_student/Dataframes/df_ngt_hlt.root")

    #─────────────────────────────────────── Define variables ────────────────────────────────────────

    df = (
        df
        #-------------- Gen -------------------

        .Define("Gen_b_idx",             f"Get_b_from_H_indices(GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")
        .Define("Gen_bJets_idx",         f"Match_b_to_GenJet(Gen_b_idx, GenPart_eta, GenPart_phi, GenJet_eta, GenJet_phi)")

        .Define("Gen_b_p4",              f"Get_b_p4(Gen_b_idx, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
        .Define("Gen_bjet_p4",           f"Get_b_p4(Gen_bJets_idx, GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass)")

        .Define("Gen_H_bb_mass",         f"Get_Gen_H_tautau_mass(Gen_b_p4)") 
        .Define("Gen_H_bbjet_mass",      f"Get_Gen_H_tautau_mass(Gen_bjet_p4)") 

        .Define("Gen_tau_idx",           f"Get_tau_from_H_indices(GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")
        .Define("tau_channel",           f"Get_tau_channel(Gen_tau_idx, GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")

        .Define("Gen_tau_p4",            f"Get_visible_tau_p4s(Gen_tau_idx, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
        .Define("Gen_tau_pt",            f"Get_p4_pt(Gen_tau_p4)")
        .Define("Gen_tau_eta",           f"Get_p4_eta(Gen_tau_p4)")
        .Define("Gen_tau_phi",           f"Get_p4_phi(Gen_tau_p4)")
        .Define("Gen_tau_mass",          f"Get_p4_mass(Gen_tau_p4)")

        .Define("Gen_H_tautau_mass",     f"Get_Gen_H_tautau_mass(Gen_tau_p4)")

        .Define("Gen_HH_mass",           f"Get_gen_mHH(Gen_tau_p4, Gen_b_p4)") 
        
        #-------------- Reco -------------------

        .Define("Matched_tau_idx",       f"deltaR_matching(Gen_tau_pt, Gen_tau_eta, Gen_tau_phi, hltHpsPFTau_pt, hltHpsPFTau_eta, hltHpsPFTau_phi, hltHpsPFTau_deepTauVSjet)")

        .Define("Reco_tau_pt",           "Get_variable(Matched_tau_idx, hltHpsPFTau_pt)")
        .Define("Reco_tau_eta",          "Get_variable(Matched_tau_idx, hltHpsPFTau_eta)")
        .Define("Reco_tau_phi",          "Get_variable(Matched_tau_idx, hltHpsPFTau_phi)")
        .Define("Reco_tau_mass",         "Get_variable(Matched_tau_idx, hltHpsPFTau_mass)")
        .Define("Reco_tau_DR",           "Get_dR_tau(Gen_tau_p4, Matched_tau_idx, hltHpsPFTau_eta, hltHpsPFTau_phi)")

        .Define("Matched_jet_idx",       "deltaR_matching_jets(Gen_bJets_idx, GenJet_eta, GenJet_phi, hltAK4PuppiJet_eta, hltAK4PuppiJet_phi, hltAK4PuppiJet_DeepFlavour_prob_b, hltAK4PuppiJet_DeepFlavour_prob_bb, hltAK4PuppiJet_DeepFlavour_prob_c, hltAK4PuppiJet_DeepFlavour_prob_g, hltAK4PuppiJet_DeepFlavour_prob_lepb, hltAK4PuppiJet_DeepFlavour_prob_uds)")

        .Define("Matched_jet_bscore",    "Get_bscore(Gen_bJets_idx, hltAK4PuppiJet_DeepFlavour_prob_b, hltAK4PuppiJet_DeepFlavour_prob_bb, hltAK4PuppiJet_DeepFlavour_prob_c, hltAK4PuppiJet_DeepFlavour_prob_g, hltAK4PuppiJet_DeepFlavour_prob_lepb, hltAK4PuppiJet_DeepFlavour_prob_uds)")
        .Define("Reco_jet_pT",           "Get_variable(Matched_jet_idx, hltAK4PuppiJet_pt)")
        .Define("Reco_jet_eta",          "Get_variable(Matched_jet_idx, hltAK4PuppiJet_eta)")
        .Define("Reco_jet_phi",          "Get_variable(Matched_jet_idx, hltAK4PuppiJet_phi)")
        .Define("Reco_jet_mass",         "Get_variable(Matched_jet_idx, hltAK4PuppiJet_mass)")
        .Define("Reco_jet_DR",           "Get_dR_jet(Gen_bJets_idx, Matched_jet_idx, GenJet_eta, GenJet_phi, hltAK4PuppiJet_eta, hltAK4PuppiJet_phi)")

        .Define("Reco_H_tautau_p4",      "Build_Higgs_p4(Matched_tau_idx, hltHpsPFTau_pt, hltHpsPFTau_eta, hltHpsPFTau_phi, hltHpsPFTau_mass)")
        .Define("Reco_H_tautau_mass",    "Get_Higgs_mass(Reco_H_tautau_p4)")

        .Define("Reco_H_bb_p4",          "Build_Higgs_p4(Matched_jet_idx, hltAK4PuppiJet_pt, hltAK4PuppiJet_eta, hltAK4PuppiJet_phi, hltAK4PuppiJet_mass)")
        .Define("Reco_H_bb_mass",        "Get_Higgs_mass(Reco_H_bb_p4)")

        .Define("Reco_HH_mass",          "Get_mHH(Reco_H_tautau_p4, Reco_H_bb_p4)")    
)
    #─────────────────────────────────────── Print column names ────────────────────────────────────────

    # all_cols = [str(c) for c in df.GetColumnNames()]
    # L1_cols = [c for c in all_cols if "L1_p" in c]
    # print("\nGenPart columns found:")
    # for c in sorted(L1_cols):
        # print(f"  {c}  [{df.GetColumnType(c)}]")

    #─────────────────────────────────────── Filter dataframe ────────────────────────────────────────
    df_NGT =    df.Filter("tau_channel >= 0 && DST_PFScouting == true")

    channel = {
        0: df.Filter("tau_channel == 0"),
        1: df.Filter("tau_channel == 1"),
        2: df.Filter("tau_channel == 2"),
    }
    
    streams = {
    "HLT_e":    channel[0].Filter("HLT_Ele32_WPTight_L1Seeded == true && HLT_Ele30_WPTight_L1Seeded_LooseDeepTauPFTauHPS30_eta2p1_CrossL1 == true"),
    "HLT_m":    channel[1].Filter("HLT_IsoMu24_FromL1TkMuon == true && HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1 == true"),
    "HLT_h":    channel[2].Filter("HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1 == true"),

    "NGT_e":    channel[0].Filter("DST_PFScouting == true"),
    "NGT_m":    channel[1].Filter("DST_PFScouting == true"),
    "NGT_h":    channel[2].Filter("DST_PFScouting == true"),
    }

    results = []
    for name, df_s in streams.items():
        n_trig    = df_s.Count()    #Get amount of events that pass each trigger
        n_tau_ok  = df_s.Filter("Matched_tau_idx[0] >= 0 && Matched_tau_idx[1] >= 0").Count()   #Get amount of events that pass each trigger and also have 2 tau's recontructed
        n_jet_ok  = df_s.Filter("Matched_jet_idx[0] >= 0 && Matched_jet_idx[1] >= 0").Count()
        n_both_ok = df_s.Filter("Matched_tau_idx[0] >= 0 && Matched_tau_idx[1] >= 0 && Matched_jet_idx[0] >= 0 && Matched_jet_idx[1] >= 0").Count()
        results.append((name, n_trig, n_tau_ok, n_jet_ok, n_both_ok))

    print(f"{'stream':8s} {'#Triggered':>13s} {'Trig. Eff':>13s} "
        f"{'#Reco tau':>13s} {'#Reco+trig. tau':>20s} "
        f"{'#Reco Jet':>13s} {'#Reco+trig. Jet':>20s} "
        f"{'#Reco both':>13s} {'#Reco+trig.both':>20s}")

    for name, n_trig, n_tau_ok, n_jet_ok, n_both_ok in results:
        ch = 0 if "_e" in name else (1 if "_m" in name else 2)
        denom = channel[ch].Count().GetValue()
        trig_count = n_trig.GetValue()

        trig_eff = trig_count / denom if denom > 0 else 0.0

        tau_abs  = n_tau_ok.GetValue()  / denom if denom > 0 else 0.0
        jet_abs  = n_jet_ok.GetValue()  / denom if denom > 0 else 0.0
        both_abs = n_both_ok.GetValue() / denom if denom > 0 else 0.0

        tau_cond  = n_tau_ok.GetValue()  / trig_count if trig_count > 0 else 0.0
        jet_cond  = n_jet_ok.GetValue()  / trig_count if trig_count > 0 else 0.0
        both_cond = n_both_ok.GetValue() / trig_count if trig_count > 0 else 0.0

        print(f"{name:8s} {trig_count:13d} {trig_eff:13f} "
            f"{tau_abs:15f} {tau_cond:15f} "
            f"{jet_abs:15f} {jet_cond:15f} "
            f"{both_abs:15f} {both_cond:15f}")