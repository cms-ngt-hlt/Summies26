import ROOT
import glob
import os
import argparse
import shutil
import tempfile
from pathlib import Path

REPO_DIR = Path(__file__).resolve().parent
SRC = REPO_DIR / "functions" / "functions.cc"
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


#─────────────────────────────────────── Get files ────────────────────────────────────────
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Argument for Gen_Analysis')
    parser.add_argument('--run-merging', action='store_true', required=False,
                         help='Run merging of NGT and HLT input files.')
    parser.add_argument('--out', default='efficiency_table.txt',
                         help='Output .txt file for the efficiency table.')
    args = parser.parse_args()

    if args.run_merging:
        os.makedirs("Dataframes", exist_ok=True)

        base_dir_hlt = "/eos/user/e/evernazz/Sarah/HHbbtautau/HLTStandard/"
        base_dir_ngt = "/eos/user/e/evernazz/Sarah/HHbbtautau/NGTScouting/"

        files_hlt = files_with_events(sorted(glob.glob(base_dir_hlt + "/step2_*.root")))
        files_ngt = files_with_events(sorted(glob.glob(base_dir_ngt + "/step2_*.root")))

        chain_hlt = make_chain(files_hlt)
        chain_ngt = make_chain(files_ngt)

        df_hlt = ROOT.RDataFrame(chain_hlt)
        df_ngt = ROOT.RDataFrame(chain_ngt)

        print(f"Total entries in chain HLT: {df_hlt.Count().GetValue()}")
        print(f"Total entries in chain NGT: {df_ngt.Count().GetValue()}")

        cols_ngt_filtered = [str(c) for c in df_ngt.GetColumnNames() if not str(c).startswith("HLT_")]
        print(f"Keeping {len(cols_ngt_filtered)} NGT columns for merging.")
        cols_hlt_filtered = [str(c) for c in df_hlt.GetColumnNames()
                              if str(c).startswith("HLT_") and not str(c).endswith("_pHLT")]
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
        df.Snapshot("Events", "/eos/user/s/sbenabde/CERN_Summer_student/Dataframes/df_ngt_hlt.root",
                    snapshot_columns, snapshot_options)
    else:
        df = ROOT.RDataFrame("Events", "/eos/user/s/sbenabde/CERN_Summer_student/Dataframes/df_ngt_hlt.root")

    #─────────────────────────────────────── Define variables ────────────────────────────────────────
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

        leg1_expr = f"{tau_var}[0] > 0"
        leg2_expr = f"{tau_var}[1] > 0"
        combined_expr = f"({leg1_expr}) && ({leg2_expr})"
        return leg1_expr, leg1_name, leg2_expr, leg2_name, combined_expr

    channels = {
        0: {"label": "electron", 
            "filter": "tau_channel == 0",
            "trig": "HLT_Ele32_WPTight_L1Seeded == true && HLT_Ele30_WPTight_L1Seeded_LooseDeepTauPFTauHPS30_eta2p1_CrossL1 == true"},
        
        1: {"label": "muon", 
            "filter": "tau_channel == 1",
            "trig": "HLT_IsoMu24_FromL1TkMuon == true && HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1 == true"},
        
        2: {"label": "hadronic", 
            "filter": "tau_channel == 2",
            "trig": "HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1 == true"},
    }

    scenarios = {
        "NoTag":   ("Matched_tau_idx",        "NoTag"),
        "TauTag":  ("Matched_tagged_tau_idx", "NoTag"),
        "bTag":    ("Matched_tau_idx",        "bTag"),
        "BothTag": ("Matched_tagged_tau_idx", "bTag"),
    }

    channel_df = {idx: df.Filter(cfg["filter"]) for idx, cfg in channels.items()}

    #─────────────────────────────────────── Book all actions (lazy) ────────────────────────────────────────
    booked_rows = []
    for idx, cfg in channels.items():
        base = channel_df[idx]
        total_n = base.Count()

        for stream_label, trig_filter in [("HLT", cfg["trig"]), ("NGT", "DST_PFScouting == true")]:
            trig_df = base.Filter(trig_filter)
            triggered_n = trig_df.Count()

            for tag_label, (tau_var, jet_key) in scenarios.items():
                leg1_expr, leg1_name, leg2_expr, leg2_name, lep_expr = leg_exprs(cfg["label"], tau_var)
                jet_expr = JET_MATCH_EXPR[jet_key]

                leg1_reco_n = trig_df.Filter(leg1_expr).Count()
                leg2_reco_n = trig_df.Filter(leg2_expr).Count()
                lep_reco_n  = trig_df.Filter(lep_expr).Count()
                jet_reco_n  = trig_df.Filter(jet_expr).Count()
                both_reco_n = trig_df.Filter(f"({lep_expr}) && ({jet_expr})").Count()

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

    #─────────────────────────────────────── Write formatted table to .txt ────────────────────────────────────────
    col_gap = "  "

    def fmt_cell(value, width, kind):
        if kind == "s":
            return f"{value:<{width}}"
        elif kind == "d":
            return f"{value:>{width}d}"
        elif kind == "pct":
            return f"{value * 100:>{width - 1}.2f}%"
        return f"{value:>{width}}"

    tag_order = ["NoTag", "TauTag", "bTag", "BothTag"]
    channel_order = ["electron", "muon", "hadronic"]
    stream_order = ["HLT", "NGT"]

    all_blocks = []
    for tag in tag_order:
        lines = [f"# Tagging: {tag}"]

        for ch_label in channel_order:
            ch_rows = [r for r in rows if r["tagging"] == tag and r["channel"] == ch_label]
            leg1_name = ch_rows[0]["leg1_name"]
            leg2_name = ch_rows[0]["leg2_name"]

            columns = [
                ("Channel",             "channel",       9,  "s"),
                ("Stream",              "stream",        6,  "s"),
                ("Total",               "total",         9,  "d"),
                ("Triggered",           "triggered",     10, "d"),
                ("Acceptance",          "trig_eff",      11, "pct"),
                (f"{leg1_name} Matched", "leg1_reco_eff", 13, "pct"),
                (f"{leg2_name} Matched", "leg2_reco_eff", 13, "pct"),
                (f"{leg1_name}+{leg2_name} Matched",      "lep_reco_eff",  15, "pct"),
                ("Jet Matched",          "jet_reco_eff",  12, "pct"),
                ("Both Matched",          "both_reco_eff", 12, "pct"),
            ]
            header_line = col_gap.join(f"{name:<{w}}" if kind == "s" else f"{name:>{w}}" for name, _, w, kind in columns)
            sep_line = "-" * len(header_line)
            lines += [header_line, sep_line]

            for i, stream in enumerate(stream_order):
                row = next(r for r in ch_rows if r["stream"] == stream)
                row_for_print = dict(row)
                row_for_print["channel"] = ch_label.capitalize() if i == 0 else ""
                lines.append(col_gap.join(fmt_cell(row_for_print[key], w, kind) for _, key, w, kind in columns))

            lines.append("")

        all_blocks.append("\n".join(lines))

    table_text = "\n\n".join(all_blocks)
    print(table_text)

    with open(args.out, "w") as f:
        f.write(table_text + "\n")

    print(f"\nTable written to {args.out}")

    df_check = (
    df.Filter("tau_channel == 2")
      .Define("old_idx0", "Get_tau_from_H_indices(GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)[0]")
      .Define("new_idx0", "Gen_tau_idx[0]")
      .Define("old_idx1", "Get_tau_from_H_indices(GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)[1]")
      .Define("new_idx1", "Gen_tau_idx[1]")
      .Define("idx_diff0", "old_idx0 - new_idx0")
      .Define("idx_diff1", "old_idx1 - new_idx1")
    )

    h0 = df_check.Histo1D(("d0", "", 21, -10, 10), "idx_diff0")
    h1 = df_check.Histo1D(("d1", "", 21, -10, 10), "idx_diff1")
    print("Mean index diff leg0:", h0.GetMean(), "  nonzero fraction:", (h0.Integral() - h0.GetBinContent(h0.FindBin(0))) / h0.Integral())
    print("Mean index diff leg1:", h1.GetMean())