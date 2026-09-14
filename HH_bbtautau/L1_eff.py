import ROOT
import os
from pathlib import Path
import shutil
import mplhep as hep
import tempfile
hep.style.use("CMS")

REPO_DIR = Path(__file__).resolve().parent

SRC = Path("/eos/user/s/sbenabde/CERN_Summer_student/HH_bbtautau/functions.cc")
HDR_DIR = SRC.parent / "headers"

SVFIT_IFACE_CC = HDR_DIR / "Tools" / "Tools" / "src" / "SVfitInterface.cc"
CMSSW_SRC = Path("/eos/home-s/sbenabde/CERN_Summer_student/CMSSW_9_4_4/src")

tmp_dir = Path(tempfile.mkdtemp(prefix="root_aclic_"))
TMP = tmp_dir / "functions.cc"

shutil.copy2(SRC, TMP)

ROOT.gSystem.AddIncludePath(f"-I{HDR_DIR}")
ROOT.gSystem.AddIncludePath(f"-I{CMSSW_SRC}")

CMSSW_LIB = Path("/eos/home-s/sbenabde/CERN_Summer_student/CMSSW_9_4_4/lib/slc7_amd64_gcc630")
ROOT.gSystem.AddDynamicPath(str(CMSSW_LIB))
ROOT.gSystem.Load("libTauAnalysisClassicSVfit.so")

ROOT.gROOT.ProcessLine(f".L {SVFIT_IFACE_CC}+")
ROOT.gROOT.ProcessLine(f".L {TMP}+")

output = "/eos/user/s/sbenabde/CERN_Summer_student/HH_bbtautau/Tables/L1_efficiencies.txt"
output2 = "/eos/user/s/sbenabde/CERN_Summer_student/HH_bbtautau/Tables/L1_list_efficiencies.txt"

#Functions
def overall_efficiency(all, trigger):
    all_c = all.Count().GetValue()
    trigger_c = trigger.Count().GetValue()
    eff = trigger_c/all_c
    return eff 

def exclusive_efficiency(df, lists_dict, target_key):
    target_paths = lists_dict[target_key] #the L1 pathss we're checking, at least one of does must fire
    
    other_paths = []    #the other ones, we don't want these to fire too
    for key, paths in lists_dict.items():
        if key != target_key: #If another list than the one we are investigating,
            other_paths += paths #take paths in that list and add them to a single list of all the other paths we don't want

    
    any_target = " || ".join(target_paths)
    none_other = " && ".join(f"!{p}" for p in other_paths) #adding ! in front of condition makes it false in Boolean logic. !A && !B means A must be false and B must be false

    condition = f"({any_target}) && ({none_other})"

    df_triggered = df.Filter(condition)
    return overall_efficiency(df, df_triggered), condition


df = ROOT.RDataFrame("Events", "/eos/user/s/sbenabde/CERN_Summer_student/Dataframes/df_pixel_tracking.root")
print(f"Total entries in df: {df.Count().GetValue()}")
    #────────────────────────────────────────────────────────────── Define variables ────────────────────────────────────────────────────────────────────────────────

df = (
    df

#───────────────────────────────────────────────────────── Gen ────────────────────────────────────────────────────────────────────────────────

#─────────────────────────────── Taus ────────────────────────────────

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

    .Define("Gen_bJets_idx",      f"Match_b_to_GenJet(Gen_b_idx, GenPart_eta, GenPart_phi, GenJet_eta, GenJet_phi)")
    .Define("Gen_bJets_flavour",  "Get_Flavour_for_Jets(Gen_bJets_idx, GenJet_hadronFlavour)")
    
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

#Filters
df_L1 = df.Filter("tau_channel == 1")  

L1_electron_mix = ['L1_pDoubleEGEle37_24_final', 'L1_pDoubleIsoTkPho22_12_final', 'L1_pDoubleTkEle25_12_final', 'L1_pIsoTkEleEGEle22_12_final', 'L1_pSingleEGEle51_final', 'L1_pSingleIsoTkEle28_final', 'L1_pSingleIsoTkPho36_final', 'L1_pSingleTkEle36_final']
L1_jet_mix = ['L1_pDoublePuppiJet112_112_final', 'L1_pDoublePuppiJet160_35_mass620_final', 'L1_pSinglePuppiJet230_final']
L1_tau_mix = ['L1_pDoublePuppiTau52_52_final']
L1_muon_mix = ['L1_pDoubleTkMuon15_7_final', 'L1_pDoubleTkMuon_4_4_OS_Dr1p2_final', 'L1_pDoubleTkMuon_4p5_4p5_OS_Er2_Mass7to18_final', 'L1_pDoubleTkMuon_OS_Er1p5_Dr1p4_final', 'L1_pSingleTkMuon22_final', 'L1_pTripleTkMuon5_3_3_final', 'L1_pTripleTkMuon_5_3_0_DoubleTkMuon_5_3_OS_MassTo9_final', 'L1_pTripleTkMuon_5_3p5_2p5_OS_Mass5to17_final']
L1_met_mix = ['L1_pPuppiHT450_final', 'L1_pPuppiMET200_final', 'L1_pPuppiMHT140_final']
L1_cross_mix = ['L1_pDoubleTkElePuppiHT_8_8_390_final', 'L1_pNNPuppiTauPuppiMet_55_190_final', 'L1_pDoubleTkMuPuppiJetPuppiMet_3_3_60_130_final', 'L1_pDoubleTkMuonTkEle5_5_9_final', 'L1_pDoubleTkMuPuppiHT_3_3_300_final', 'L1_pPuppiTauTkIsoEle45_22_final', 'L1_pPuppiTauTkMuon42_18_final', 'L1_pTkEleIsoPuppiHT_26_190_final', 'L1_pTkElePuppiJet_28_40_MinDR_final', 'L1_pTkEleTkMuon10_20_final', 'L1_pTkMuPuppiJetPuppiMet_3_110_120_final', 'L1_pTkMuTriPuppiJet_12_40_dRMax_DoubleJet_dEtaMax_final', 'L1_pTkMuonDoubleTkEle6_17_17_final', 'L1_pTkMuonPuppiHT6_320_final', 'L1_pTkMuonTkEle7_23_final', 'L1_pTkMuonTkIsoEle7_20_final']

# L1_pQuadJet70_55_40_40_final and L1_pPuppiHT400_final not available in our files
data = {}
all_cols = [str(c) for c in df.GetColumnNames()]
L1_cols  = [c for c in all_cols if "L1_p" in c]

print("\nComputing the efficiencies...")
for c in sorted(L1_cols):
    column_name = f"{c}  [{df.GetColumnType(c)}]"
    df_L1_filtered = df_L1.Filter(f"{c} == true")
    efficiency_L1 = str(overall_efficiency(df_L1, df_L1_filtered))
    data[column_name] = efficiency_L1

#Generate output table with efficiencies for each L1 path
sorted_data = {k: v for k, v in sorted(data.items(), key=lambda item: item[1])}
with open(output, "w") as out:
    for key, value in sorted_data.items():
        v = float(value) * 100
        out.write(f"{key:<70} {v:.2f} %\n")
print(f"Output saved in {output}")

# filters = {}
# for s_e, s_jet, s_tau, s_mu, s_met, s_mix in L1_electron_mix, L1_jet_mix, L1_tau_mix, L1_muon_mix, L1_met_mix, L1_cross_mix:
#     name = df_{string} 
#     L1_condition_e = f"{string} == true"
#     L1_condition_e = f"{string} == true"
#     triggered_df = df.Filter(L1_condition)

L1_lists = {
    "electron": L1_electron_mix,
    "jet": L1_jet_mix,
    "tau": L1_tau_mix,
    "muon": L1_muon_mix,
    "MET": L1_met_mix,
    # "cross": L1_cross_mix,
}

results = {}
for key in L1_lists:
    eff, condition = exclusive_efficiency(df_L1, L1_lists, key)
    results[key] = eff
    print(f"{key:10s}:eff = {eff*100:.2f}%\n")

with open(output2, "a") as f:
    for key, eff in results.items():
        f.Write("Muon channel")
        f.write(f"Exclusive efficiency ({key} only): {eff*100:.2f}% \n")
print(f"Output saved in {output2}")
