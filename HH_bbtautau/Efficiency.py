import ROOT
import glob
import os
import argparse
import shutil
import tempfile
from pathlib import Path
import re

REPO_DIR = Path(__file__).resolve().parent
SRC = REPO_DIR  / "functions.cc"
HDR_DIR = REPO_DIR / "headers"

tmp_dir = Path(tempfile.mkdtemp(prefix="root_aclic_"))
TMP = tmp_dir / "functions.cc"
shutil.copy2(SRC, TMP)
ROOT.gSystem.AddIncludePath(f"-I{HDR_DIR}")
ROOT.gROOT.ProcessLine(f".L {TMP}+")


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

import re

def filter_by_index(file_list, min_idx=0, max_idx=19):
    filtered = []
    for f in file_list:
        m = re.search(r"step2_(\d+)\.root$", f)
        if m and min_idx <= int(m.group(1)) <= max_idx:
            filtered.append(f)
    filtered.sort(key=lambda f: int(re.search(r"step2_(\d+)\.root$", f).group(1)))
    return filtered

#─────────────────────────────────────── Get files ────────────────────────────────────────
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Argument for Gen_Analysis')
    parser.add_argument('--run-merging', action='store_true', required=False, help='Run merging of NGT and HLT input files.')
    parser.add_argument('--out', default='efficiency_table.txt', help='Output .txt file for the efficiency table.')
    args = parser.parse_args()

    if args.run_merging:
        os.makedirs("Dataframes", exist_ok=True)

        base_dir_hlt = "/eos/user/e/evernazz/Sarah/HHbbtautau/HLTStandard/"
        base_dir_ngt = "/eos/user/e/evernazz/Sarah/HHbbtautau/NGTScouting_ImprovedPixelTracking"

        files_hlt_all = glob.glob(base_dir_hlt + "/step2_*.root")
        files_ngt_all = glob.glob(base_dir_ngt + "/step2_*.root")

        files_hlt = files_with_events(filter_by_index(files_hlt_all, 0, 19))
        files_ngt = files_with_events(filter_by_index(files_ngt_all, 0, 19))

        chain_hlt = make_chain(files_hlt)
        chain_ngt = make_chain(files_ngt)

        df_hlt = ROOT.RDataFrame(chain_hlt)
        df_ngt = ROOT.RDataFrame(chain_ngt)

        print(f"Total entries in chain HLT: {df_hlt.Count().GetValue()}")
        print(f"Total entries in chain NGT: {df_ngt.Count().GetValue()}")

        cols_ngt_filtered = [str(c) for c in df_ngt.GetColumnNames() if not str(c).startswith("HLT_")]
        print(f"Keeping {len(cols_ngt_filtered)} NGT columns for merging.")
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
        snapshot_options.fCompressionLevel = 4
        df.Snapshot("Events", "/eos/user/s/sbenabde/CERN_Summer_student/Dataframes/df_pixel_tracking.root", snapshot_columns, snapshot_options)
    else:
        df = ROOT.RDataFrame("Events", "/eos/user/s/sbenabde/CERN_Summer_student/Dataframes/df_pixel_tracking.root")

    #─────────────────────────────────────── Define variables ────────────────────────────────────────
    df = (
        df
        #─────────────────────────────────────────────────────────────────── Gen ────────────────────────────────────────────────────────────────────────────────

        .Define("Gen_taus_idx",       f"Get_tau_from_H_indices(GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")
            .Define("Gen_tau_idx",        f"Get_taus(Gen_taus_idx, GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother, GenPart_pt)")
            
            .Define("tau_channel",        f"Get_tau_channel(Gen_taus_idx, GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")
            .Define("tau_first_pdg",      f"GenPart_pdgId[Gen_tau_idx[0]]")
            .Define("tau_2nd_pdg",        f"GenPart_pdgId[Gen_tau_idx[1]]")
            
            .Define("Gen_tau_p4",         f"Get_visible_tau_p4s(Gen_tau_idx, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
            .Define("Gen_tau_pt",         f"Get_p4_pt(Gen_tau_p4)")
            .Define("Gen_tau_eta",        f"Get_p4_eta(Gen_tau_p4)")
            .Define("Gen_tau_phi",        f"Get_p4_phi(Gen_tau_p4)")
            .Define("Gen_tau_mass",       f"Get_p4_mass(Gen_tau_p4)")
        
            .Define("Gen_tau_all_p4",         f"Get_all_taus_p4(Gen_taus_idx, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
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
            
            #─────────────────────────────── Jets ────────────────────────────────
        
            .Define("Gen_b_idx_first", f"Get_b_from_H_indices(GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother, true)")
            .Define("Gen_b_idx",       f"Get_b_from_H_indices(GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother, false)")
            # .Define("Daughters",      f"GetPDG(Gen_b_idx, GenPart_pdgId)")
            # .Define("Mothers",        f"GetPDG_mother(Gen_b_idx, GenPart_pdgId,GenPart_genPartIdxMother)")
        
            .Define("n_Gen_b",        f"Gen_b_idx.size()")
        
            .Define("Gen_bJets_idx",  f"Match_b_to_GenJet(Gen_b_idx, GenPart_eta, GenPart_phi, GenJet_eta, GenJet_phi)")
            .Define("Gen_bJets_flavour", "Get_Flavour_for_Jets(Gen_bJets_idx, GenJet_hadronFlavour)")
            
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
        
        
        #─────────────────────────────── Muons ────────────────────────────────
        
            .Define("Gen_Muon_idx",       f"Get_lepton_idx(13, GenPart_pdgId, Gen_taus_idx, GenPart_statusFlags, GenPart_genPartIdxMother)")
            .Define("nGen_Muon",          f"Gen_Muon_idx.size()")
            
            .Define("Gen_muon_pt",        f"Get_variable(Gen_Muon_idx, GenPart_pt)")
            .Define("Gen_muon_eta",       f"Get_variable(Gen_Muon_idx, GenPart_eta)")
            .Define("Gen_muon_phi",       f"Get_variable(Gen_Muon_idx, GenPart_phi)")
            .Define("Gen_muon_mass",      f"Get_variable(Gen_Muon_idx, GenPart_mass)")
        
        #─────────────────────────────── Electrons ────────────────────────────────
        
            .Define("Gen_Electron_idx",   f"Get_lepton_idx(11, GenPart_pdgId, Gen_taus_idx, GenPart_statusFlags, GenPart_genPartIdxMother)")
            .Define("nGen_Electron",      f"Gen_Electron_idx.size()")
        
            .Define("Gen_electron_pt",    f"Get_variable(Gen_Electron_idx, GenPart_pt)")
            .Define("Gen_electron_eta",   f"Get_variable(Gen_Electron_idx, GenPart_eta)")
            .Define("Gen_electron_phi",   f"Get_variable(Gen_Electron_idx, GenPart_phi)")
            .Define("Gen_electron_mass",  f"Get_variable(Gen_Electron_idx, GenPart_mass)")
           
           #─────────────────────────────── Higgs ────────────────────────────────
        
            .Define("Gen_H_bb",             f"Get_H_fromdecay_idx(Gen_b_idx_first, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)") 
            .Define("Gen_H_bb_mass",        f"Gen_H_bb.M()") 
        
            .Define("Gen_H_bbjets",         f"Get_H_fromdecay(Gen_bjet_p4)") 
            .Define("Gen_H_bbjets_mass",    f"Gen_H_bbjets.M()") 
        
            .Define("Gen_H_tt_all",         f"Get_H_fromdecay(Gen_tau_all_p4)")
            .Define("Gen_H_tt",             f"Get_H_fromdecay(Gen_tau_p4)")
            .Define("Gen_H_tt_all_mass",    f"Gen_H_tt_all.M()")
            .Define("Gen_H_tt_mass",        f"Gen_H_tt.M()")
        
            .Define("Gen_HH_firstb",        f"Get_HH(Gen_H_bb, Gen_H_tt)") 
            .Define("Gen_HH",               f"Get_HH(Gen_H_bbjets, Gen_H_tt)") 
            .Define("Gen_mHH_firstb",       f"Gen_HH_firstb.M()") 
            .Define("Gen_mHH",              f"Gen_HH.M()") 
            
        #────────────────────────────────────────────────────────────── Reco HLT ────────────────────────────────────────────────────────────────────────────────
        #─────────────────────────────── Taus ────────────────────────────────
        
            .Define("Matched_tau_idx",          f"deltaR_matching(Gen_tau_idx, GenPart_pdgId, Gen_tau_pt, Gen_tau_eta, Gen_tau_phi, hltHpsPFTau_pt, hltHpsPFTau_eta, hltHpsPFTau_phi, hltElectron_pt, hltElectron_eta, hltElectron_phi, hltMuon_pt, hltMuon_eta, hltMuon_phi, hltHpsPFTau_deepTauVSjet, false)")
            .Define("Matched_tagged_tau_idx",   f"deltaR_matching(Gen_tau_idx, GenPart_pdgId, Gen_tau_pt, Gen_tau_eta, Gen_tau_phi, hltHpsPFTau_pt, hltHpsPFTau_eta, hltHpsPFTau_phi, hltElectron_pt, hltElectron_eta, hltElectron_phi, hltMuon_pt, hltMuon_eta, hltMuon_phi, hltHpsPFTau_deepTauVSjet, true)")
        
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
        
        #─────────────────────────────── Tagged taus ────────────────────────────────
        
            .Define("Reco_tagged_tau_pt",               "Get_variable(Matched_tagged_tau_idx, hltHpsPFTau_pt)")
            .Define("Reco_tagged_tau_eta",              "Get_variable(Matched_tagged_tau_idx, hltHpsPFTau_eta)")
            .Define("Reco_tagged_tau_phi",              "Get_variable(Matched_tagged_tau_idx, hltHpsPFTau_phi)")
            .Define("Reco_tagged_tau_mass",             "Get_variable(Matched_tagged_tau_idx, hltHpsPFTau_mass)")
            .Define("Reco_tagged_tau_DR",               "Get_dR_tau(Gen_tau_p4, Matched_tagged_tau_idx, hltHpsPFTau_eta, hltHpsPFTau_phi)")
        
            .Define("Reco_leading_tagged_tau_pt",       "Reco_tagged_tau_pt[0]")
            .Define("Reco_leading_tagged_tau_eta",      "Reco_tagged_tau_eta[0]")
            .Define("Reco_leading_tagged_tau_phi",      "Reco_tagged_tau_phi[0]")   
            .Define("Reco_leading_tagged_tau_mass",     "Reco_tagged_tau_mass[0]")   
            .Define("Reco_leading_tagged_tau_DR",       "Reco_tagged_tau_DR[0]") 
        
            .Define("Reco_subleading_tagged_tau_pt",    "Reco_tagged_tau_pt[1]")
            .Define("Reco_subleading_tagged_tau_eta",   "Reco_tagged_tau_eta[1]")
            .Define("Reco_subleading_tagged_tau_phi",   "Reco_tagged_tau_phi[1]")   
            .Define("Reco_subleading_tagged_tau_mass",  "Reco_tagged_tau_mass[1]")   
            .Define("Reco_subleading_tagged_tau_DR",    "Reco_tagged_tau_DR[1]")    
        
        #─────────────────────────────── Jets ────────────────────────────────
        
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
        
        #─────────────────────────────── b-jets ────────────────────────────────
        
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
        
        #─────────────────────────────── Muons & Electrons ────────────────────────────────
        
            .Define("Matched_muon_idx",     f"deltaR_matching_lepton(Gen_muon_pt, Gen_muon_eta, Gen_muon_phi, hltMuon_pt, hltMuon_eta, hltMuon_phi)")
            .Define("Matched_electron_idx", f"deltaR_matching_lepton(Gen_electron_pt, Gen_electron_eta, Gen_electron_phi, hltElectron_pt, hltElectron_eta, hltElectron_phi)")
         
            .Define("Reco_electron_pt",      "Get_variable(Matched_electron_idx, hltElectron_pt)")
            .Define("Reco_muon_pt",          "Get_variable(Matched_muon_idx, hltMuon_pt)")
        
        #─────────────────────────────── Higgs ────────────────────────────────
        
            .Define("Reco_H_tt",        "Build_Higgs_p4(Matched_tau_idx, hltHpsPFTau_pt, hltHpsPFTau_eta, hltHpsPFTau_phi, hltHpsPFTau_mass)")
            .Define("Reco_H_bb",        "Build_Higgs_p4(Matched_jet_idx, hltAK4PuppiJet_pt, hltAK4PuppiJet_eta, hltAK4PuppiJet_phi, hltAK4PuppiJet_mass)")
            
            .Define("Reco_H_tt_mass",   "Reco_H_tt.M()")
            .Define("Reco_H_bb_mass",   "Reco_H_bb.M()")
        
            .Define("Reco_HH",          "Get_HH(Reco_H_tt, Reco_H_bb)")   
            .Define("Reco_HH_mass",     "Reco_HH.M()")   
        
        )

    #─────────────────────────────────────── Channel + trigger + leg definitions ────────────────────────────────────────
    JET_MATCH_EXPR = {
        "NoTag": "Matched_jet_idx[0] >= 0 && Matched_jet_idx[1] >= 0",
        "bTag":  "Matched_bjet_idx[0] >= 0 && Matched_bjet_idx[1] >= 0",
    }

    def leg_exprs(channel_label, tau_var):
        if channel_label == "electron":
            leg1_name = "e"
            leg2_name = "τ"
        elif channel_label == "muon":
            leg1_name = "μ"
            leg2_name = "τ"
        else:  # hadronic
            leg1_name = "τ1"
            leg2_name = "τ2"

        leg1_expr = f"{tau_var}[0] >= 0"
        leg2_expr = f"{tau_var}[1] >= 0"
        combined_expr = f"({leg1_expr}) && ({leg2_expr})"
        return leg1_expr, leg1_name, leg2_expr, leg2_name, combined_expr

    channels = {
        0: {"label": "electron", "filter": "tau_channel == 0", "trig": "HLT_Ele32_WPTight_L1Seeded == true || HLT_Ele30_WPTight_L1Seeded_LooseDeepTauPFTauHPS30_eta2p1_CrossL1 == true", "trig_L1": "L1_pSingleEGEle51_final == true || L1_pSingleTkEle36_final == true || L1_pSingleIsoTkEle28_final == true || L1_pPuppiTauTkIsoEle45_22_final == true"},
        1: {"label": "muon", "filter": "tau_channel == 1", "trig": "HLT_IsoMu24_FromL1TkMuon == true || HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1 == true", "trig_L1": "L1_pSingleTkMuon22_final == true || L1_pPuppiTauTkMuon42_18_final == true"},
        2: {"label": "hadronic", "filter": "tau_channel == 2", "trig": "HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1 == true", "trig_L1": "L1_pDoublePuppiTau52_52_final == true"},
    }

    scenarios = {
        "NoTag":   ("Matched_tau_idx",        "NoTag"), #No tau tag, no b tag
        "TauTag":  ("Matched_tagged_tau_idx", "NoTag"), #tau tag, no b tag
        "bTag":    ("Matched_tau_idx",        "bTag"),  #no tau tag, b tag
        "BothTag": ("Matched_tagged_tau_idx", "bTag"),  #tau tag, b tag
    }

    channel_df = {idx: df.Filter(cfg["filter"]) for idx, cfg in channels.items()}
    channel_order = ["electron", "muon", "hadronic"] 

    #─────────────────────────────────────── Book all actions (lazy) ────────────────────────────────────────
    booked_rows = []
    for idx, cfg in channels.items():
        base = channel_df[idx]
        total_n = base.Count()

        for stream_label, trig_filter in [("L1", cfg["trig_L1"]), ("HLT", cfg["trig"]), ("NGT", "DST_PFScouting == true")]:
            trig_df = base.Filter(trig_filter)
            triggered_n = trig_df.Count()

            for tag_label, (tau_var, jet_key) in scenarios.items():
                leg1_expr, leg1_name, leg2_expr, leg2_name, lep_expr = leg_exprs(cfg["label"], tau_var)
                jet_expr = JET_MATCH_EXPR[jet_key]

                leg1_reco_n = trig_df.Filter(leg1_expr).Count() #number of triggered e/mu/tau
                leg2_reco_n = trig_df.Filter(leg2_expr).Count() #number of triggered hadronic tau
                lep_reco_n  = trig_df.Filter(lep_expr).Count()  #number of triggered e/mu/tau AND tau h
                jet_reco_n  = trig_df.Filter(jet_expr).Count()  #number of triggered jet
                both_reco_n = trig_df.Filter(f"({lep_expr}) && ({jet_expr})").Count() #number of both triggered

                booked_rows.append({
                    "stream": stream_label,
                    "channel": cfg["label"],
                    "tagging": tag_label,
                    "leg1_name": leg1_name,
                    "leg2_name": leg2_name,
                    "total": total_n,
                    "triggered": triggered_n,
                    "leg1_reco": leg1_reco_n,
                    "leg2_reco": leg2_reco_n,
                    "lep_reco": lep_reco_n,
                    "jet_reco": jet_reco_n,
                    "both_reco": both_reco_n,
                })

    #─────────────────────────────────────── Resolve values + compute fractions ────────────────────────────────────────
    def safe_div(a, b):
        return a / b if b > 0 else 0.0

    rows = []
    for r in booked_rows:
        total     = r["total"].GetValue()
        triggered = r["triggered"].GetValue()
        leg1_reco = r["leg1_reco"].GetValue()
        leg2_reco = r["leg2_reco"].GetValue()
        lep_reco  = r["lep_reco"].GetValue()
        jet_reco  = r["jet_reco"].GetValue()
        both_reco = r["both_reco"].GetValue()

        rows.append({
            "tagging":   r["tagging"],
            "stream":    r["stream"],
            "channel":   r["channel"],
            "leg1_name": r["leg1_name"],
            "leg2_name": r["leg2_name"],
            "total":     total,
            "triggered": triggered,
            "trig_eff":      safe_div(triggered, total),
            "leg1_reco_eff": safe_div(leg1_reco, triggered),
            "leg2_reco_eff": safe_div(leg2_reco, triggered),
            "lep_reco_eff":  safe_div(lep_reco,  triggered),
            "jet_reco_eff":  safe_div(jet_reco,  triggered),
            "both_reco_eff": safe_div(both_reco, triggered),
        })

     #─────────────────────────────────────── Write merged LaTeX table ────────────────────────────────────────
    def get_row(channel, stream, tagging):
        return next(r for r in rows if r["channel"] == channel and r["stream"] == stream and r["tagging"] == tagging)

    def pct(x):
        return f"{x * 100:.2f} \\%"

    channel_tex_labels = {
        "electron": r"$\tau\tau \to e\tau_h$",
        "muon":     r"$\tau\tau \to \mu\tau_h$",
        "hadronic": r"$\tau\tau \to \tau_h\tau_h$",
    }
    stream_order_tex = ["L1", "HLT", "NGT"]

    tex_lines = []
    tex_lines.append(r"\begingroup")
    tex_lines.append(r"\renewcommand{\arraystretch}{1.3}")
    tex_lines.append(r"\begin{table}[H]")
    tex_lines.append(r"    \caption{Total number of triggered events, $\Delta R$-matching efficiency (untagged and tagged) for L1, HLT and NGT streams, for the 3 different $\tau$ decay channels.}")
    tex_lines.append(r"    \centering")
    tex_lines.append(r"    \resizebox{175mm}{!}{")
    tex_lines.append(r"    \begin{tabular}{c c|c||c|c|c|c||c||c|c|c}")
    tex_lines.append(r"     &  &  & \multicolumn{5}{c|}{Matched} & \multicolumn{2}{c|}{Tagged} & Efficiency \\ ")
    tex_lines.append(r"     Channel & Triggered & Stream & $e$/$\mu$/$\tau$ & $\tau_h$ & both leptons & $b$-jet & Both & lepton & $b$-jet & \\ \hline")

    for i, ch_label in enumerate(channel_order):
        tex_lines.append(r"    \multirow{3}{*}{%s}" % channel_tex_labels[ch_label])
        for j, stream in enumerate(stream_order_tex):
            r_notag = get_row(ch_label, stream, "NoTag")
            r_btag  = get_row(ch_label, stream, "bTag")
            r_tautag = get_row(ch_label, stream, "TauTag")
            r_both  = get_row(ch_label, stream, "BothTag")

            triggered = r_notag["triggered"]
            e_mu_tau  = pct(r_notag["leg1_reco_eff"])
            tau_h     = pct(r_notag["leg2_reco_eff"])
            both_lep  = pct(r_notag["lep_reco_eff"])
            jet       = pct(r_notag["jet_reco_eff"])       # plain jet match (was r_btag)
            both_m    = pct(r_notag["both_reco_eff"])      # untagged lep + plain jet (was r_btag)
            tag_lep   = pct(r_tautag["lep_reco_eff"])
            tag_bjet  = pct(r_btag["jet_reco_eff"])         # b-tagged jet match (was "-")
            eff       = pct(r_both["both_reco_eff"])

            row_str = (f"       & {triggered}  & {stream} & {e_mu_tau} & {tau_h} & {both_lep} "
                       f"& {jet} & {both_m} & {tag_lep} & {tag_bjet} & {eff} \\\\")
            tex_lines.append(row_str)

        if i < len(channel_order) - 1:
            tex_lines.append(r"    \hline")

    tex_lines.append(r"    \end{tabular}")
    tex_lines.append(r"    }")
    tex_lines.append(r"    \label{tab:HH_recoeff}")
    tex_lines.append(r"\end{table}")
    tex_lines.append(r"\endgroup")

    tex_text = "\n".join(tex_lines)
    print(tex_text)

    tex_out = args.out.rsplit(".", 1)[0] + ".tex"
    with open(tex_out, "w") as f:
        f.write(tex_text + "\n")

    print(f"\nLaTeX table written to {tex_out}")
