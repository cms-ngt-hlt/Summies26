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

SRC = "/eos/user/s/sbenabde/CERN_Summer_student/HH_bbtautau/functions/functions.cc"
HDR = "/eos/user/s/sbenabde/CERN_Summer_student/HH_bbtautau/headers/functions.h"
TMP = "/tmp/functions.cc"

shutil.copy(SRC, TMP)
ROOT.gROOT.ProcessLine(f'.L {TMP}+')

ROOT.ROOT.EnableImplicitMT()

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

def make_panel(ax, col, cfg):
    counts, edges = np.histogram(data[col], bins=cfg["bins"], range=cfg["range"])
    hep.histplot(counts, edges, ax=ax, label=cfg["label"], color=cfg["color"], hatch=cfg["hatch"], histtype=cfg["histtype"])
    ax.set_xlabel(cfg["xlabel"])
    ax.set_ylabel(cfg["ylabel"])
    ax.legend(loc="upper right")
    hep.cms.label("Preliminary", data=False, ax=ax, com=14)
    return counts, edges

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

    streams = {
    "HLT_e":    df.Filter("tau_channel == 0 && HLT_Ele32_WPTight_L1Seeded == true && HLT_Ele30_WPTight_L1Seeded_LooseDeepTauPFTauHPS30_eta2p1_CrossL1 == true"),
    "HLT_m":    df.Filter("tau_channel == 1 && HLT_IsoMu24_FromL1TkMuon == true && HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1 == true"),
    "HLT_h":    df.Filter("tau_channel == 2 && HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1 == true"),

    "NGT_e":    df.Filter("tau_channel == 0 && DST_PFScouting == true"),
    "NGT_m":    df.Filter("tau_channel == 1 && DST_PFScouting == true"),
    "NGT_h":    df.Filter("tau_channel == 2 && DST_PFScouting == true"),
    }

    df_NGT =    df.Filter("tau_channel >= 0 && DST_PFScouting == true")

    results = []
    for name, df_s in streams.items():
        n_trig    = df_s.Count()
        n_tau_ok  = df_s.Filter("Matched_tau_idx[0] >= 0 && Matched_tau_idx[1] >= 0").Count()
        n_jet_ok  = df_s.Filter("Matched_jet_idx[0] >= 0 && Matched_jet_idx[1] >= 0").Count()
        n_both_ok = df_s.Filter("Matched_tau_idx[0] >= 0 && Matched_tau_idx[1] >= 0 && "
                                "Matched_jet_idx[0] >= 0 && Matched_jet_idx[1] >= 0").Count()
        results.append((name, n_trig, n_tau_ok, n_jet_ok, n_both_ok))

    channel = {
        0: df.Filter("tau_channel == 0"),
        1: df.Filter("tau_channel == 1"),
        2: df.Filter("tau_channel == 2"),
    }

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

# ─────────────────────────────────────── Define and plot histograms 1D ────────────────────────────────────────

    os.makedirs("plots/kinematics", exist_ok=True)
    
    channel_suffix = {0: "e", 1: "m", 2: "h"}
    channel_title = {
        0: r"$\tau \tau \to \tau_e \tau_h$",
        1: r"$\tau \tau \to \tau_{\mu} \tau_h$",
        2: r"$\tau \tau \to \tau_h \tau_h$",
    }
 
    var_configs = {
        # "pt":     {"branch": "Gen_tau_pt",  "xlabel": r"$p_{T}(\tau)$ [GeV]",  "bins": 30, "range": (0, 300)},
        # "eta":    {"branch": "Gen_tau_eta", "xlabel": r"$\eta(\tau)$",         "bins": 30, "range": (-3, 3)},
        "bscore": {"branch": "Matched_jet_bscore", "xlabel": r"$b-score$",  "bins": 30, "range": (0, 1)},
    }

    def get_leg_array(df_in, colname):
        """Pull a per-event 2-vector column out of the RDataFrame and flatten
        leg 0 (e.g. leading) + leg 1 (e.g. subleading) into one 1D array."""
        col = df_in.AsNumpy([colname])[colname]
        leg0 = np.array([ev[0] for ev in col if len(ev) > 0])
        leg1 = np.array([ev[1] for ev in col if len(ev) > 1])
        return np.concatenate([leg0, leg1])
 
    data = {}
    for var, vcfg in var_configs.items():
        branch = vcfg["branch"]
        for ch, suf in channel_suffix.items():
            data[f"{var}_all_{suf}"] = get_leg_array(channel[ch],           branch)
            data[f"{var}_HLT_{suf}"] = get_leg_array(streams[f"HLT_{suf}"], branch)
            data[f"{var}_NGT_{suf}"] = get_leg_array(streams[f"NGT_{suf}"], branch)
            data[f"{var}_NGT_{suf}"] = get_leg_array(streams[f"NGT_{suf}"], branch)
 
    selection_style = {
        "all": {"label": "All events",           "color": "mediumvioletred",  "histtype": "step", "hatch": ""},
        "HLT": {"label": "Events that pass HLT", "color": "purple",           "histtype": "fill", "hatch": ""},
        "NGT": {"label": "Events that pass NGT", "color": "darkorchid",       "histtype": "step", "hatch": "//"},
    }
 
    def plot_kinematic(var, ch):
        suf   = channel_suffix[ch]
        vcfg  = var_configs[var]
        bins, rng = vcfg["bins"], vcfg["range"]
 
        fig, ax = plt.subplots(figsize=(10, 8))
 
        counts = {}
        edges = None
        for sel, cfg in selection_style.items():
            arr = data[f"{var}_{sel}_{suf}"]
            c, e = np.histogram(arr, bins=bins, range=rng)
            counts[sel] = c
            edges = e
            hep.histplot(c, e, ax=ax, label=cfg["label"], color=cfg["color"], histtype=cfg["histtype"], hatch=cfg["hatch"])
 
        ax.set_xlabel(vcfg["xlabel"])
        ax.set_ylabel("Events")
        ax.legend(loc="upper right")
        ax.text(0.125, -0.08, channel_title[ch], transform=ax.transAxes,
                 fontsize=20, ha="right", va="top")
        hep.cms.label("Preliminary", data=False, ax=ax, com=14)
 
        centers     = 0.5 * (edges[:-1] + edges[1:])
        bin_widths  = np.diff(edges)
        den         = counts["all"]
 
        def eff_and_err(num, den):
            eff = np.divide(num, den, out=np.zeros_like(num, dtype=float), where=den > 0)
            err = np.array([np.sqrt(e * (1 - e) / d) if d > 0 else 0.0 for e, d in zip(eff, den)])
            return eff, err
 
        eff_NGT, err_NGT = eff_and_err(counts["NGT"], den)
 
        eff_color = "olivedrab"
        ax2 = ax.twinx()
        ax2.errorbar(centers, eff_NGT, xerr=0.5 * bin_widths, yerr=err_NGT, color='yellowgreen', fmt='*', label="NGT efficiency")
        ax2.set_ylabel("Efficiency", color=eff_color)
        ax2.tick_params(axis='y', labelcolor=eff_color)
        ax2.set_ylim(0, 1)
        ax2.legend(loc="upper right")
 
        plt.tight_layout()
        outname = f"plots/kinematics/{var}_{suf}.png"
        plt.savefig(outname, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"Saved {outname}")
 
    for var in var_configs:
        for ch in channel_suffix:
            plot_kinematic(var, ch)


    # data = {

    #     "gen_b_massH_NGT":  streams.get("NGT_h").AsNumpy(["Gen_H_bb_mass"])["Gen_H_bb_mass"],  
    #     "gen_bJet_massH_NGT":  streams.get("NGT_h").AsNumpy(["Gen_H_bbjet_mass"])["Gen_H_bbjet_mass"],  

    #     "gen_mH_NGT_e":  streams.get("NGT_e").AsNumpy(["Gen_H_tautau_mass"])["Gen_H_tautau_mass"],  
    #     "gen_mH_NGT_m":  streams.get("NGT_m").AsNumpy(["Gen_H_tautau_mass"])["Gen_H_tautau_mass"],  
    #     "gen_mH_NGT_h":  streams.get("NGT_h").AsNumpy(["Gen_H_tautau_mass"])["Gen_H_tautau_mass"], 

    #     "mH_NGT_e":  streams.get("NGT_e").AsNumpy(["Reco_H_tautau_mass"])["Reco_H_tautau_mass"],  
    #     "mH_NGT_m":  streams.get("NGT_m").AsNumpy(["Reco_H_tautau_mass"])["Reco_H_tautau_mass"],  
    #     "mH_NGT_h":  streams.get("NGT_h").AsNumpy(["Reco_H_tautau_mass"])["Reco_H_tautau_mass"],  

    #     "mHb_all":  channel.get(0).AsNumpy(["Reco_H_bb_mass"])["Reco_H_bb_mass"],
    #     "mHb_HLT":  streams.get("HLT_m").AsNumpy(["Reco_H_bb_mass"])["Reco_H_bb_mass"],
    #     "mHb_NGT":  streams.get("NGT_m").AsNumpy(["Reco_H_bb_mass"])["Reco_H_bb_mass"],    

    #     "mHH_all":  channel.get(0).AsNumpy(["Reco_HH_mass"])["Reco_HH_mass"],
    #     "mHH_HLT":  streams.get("HLT_m").AsNumpy(["Reco_HH_mass"])["Reco_HH_mass"],
    #     "mHH_NGT":  streams.get("NGT_m").AsNumpy(["Reco_HH_mass"])["Reco_HH_mass"], 
        
    #     "gen_mH_NGT":  df_NGT.AsNumpy(["Gen_H_tautau_mass"])["Gen_H_tautau_mass"],  
    #     "mH_NGT":      df_NGT.AsNumpy(["Reco_H_tautau_mass"])["Reco_H_tautau_mass"],  

    #     "gen_mHb_NGT":  df_NGT.AsNumpy(["Gen_H_bb_mass"])["Gen_H_bb_mass"],  
    #     "mHb_NGT":      df_NGT.AsNumpy(["Reco_H_bb_mass"])["Reco_H_bb_mass"], 

    #     "gen_mHH_NGT":  df_NGT.AsNumpy(["Gen_HH_mass"])["Gen_HH_mass"],  
    #     "mHH_NGT":      df_NGT.AsNumpy(["Reco_HH_mass"])["Reco_HH_mass"],  

    # }


    # plot_configs = {       

        # "gen_b_massH_NGT":         {"label": r"$H \to bb$ from gen b quarks",   "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events that pass NGT Scouting", "bins": 30, "range": (0, 180), "color": "purple", "histtype": "fill", "hatch": ""},
        # "gen_bJet_massH_NGT":      {"label": r"$H \to bb$ from gen b jets",     "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events that pass NGT Scouting", "bins": 30, "range": (0, 180), "color": "mediumvioletred", "histtype": "step", "hatch": "/"},

        # "mH_NGT_e":      {"label": r"$\tau \tau \to \tau_e \tau_h$",   "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 150), "color": "purple", "histtype": "step", "hatch": ""},
        # "mH_NGT_m":      {"label": r"$\tau \tau \to \tau_\mu \tau_h$",   "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 150), "color": "mediumvioletred", "histtype": "step", "hatch": ""},
        # "mH_NGT_h":      {"label": r"$\tau \tau \to \tau_h \tau_h$",   "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 150), "color": "orchid", "histtype": "step", "hatch": ""},

        # "gen_mH_NGT_e":      {"label": r"GEN $\tau \tau \to \tau_e \tau_h$",   "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 150), "color": "deepskyblue", "histtype": "fill", "hatch": ""},
        # "gen_mH_NGT_m":      {"label": r"GEN $\tau \tau \to \tau_\mu \tau_h$", "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 150), "color": "darkcyan", "histtype": "fill", "hatch": ""},
        # "gen_mH_NGT_h":      {"label": r"GEN $\tau \tau \to \tau_h \tau_h$",   "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 150), "color": "royalblue", "histtype": "fill", "hatch": ""},

        # "gen_mH_NGT":  {"label": r"GEN $H \to \tau \tau$",   "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 150), "color": "purple", "histtype": "fill", "hatch": ""},
        # "mH_NGT":      {"label": r"RECO $H \to \tau \tau$",  "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 150), "color": "mediumvioletred", "histtype": "fill", "hatch": ""},

        # "gen_mHb_NGT":  {"label": r"Gen $H \to bb$",   "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events that pass NGT Scouting", "bins": 30, "range": (0, 180), "color": "purple", "histtype": "fill", "hatch": ""},
        # "mHb_NGT":      {"label": r"Reco $H \to bb$",  "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events that pass NGT Scouting", "bins": 30, "range": (0, 180), "color": "mediumvioletred", "histtype": "fill", "hatch": ""},

        # "gen_mHH_NGT":  {"label": r"Gen. $HH \to bb \tau \tau$",   "xlabel": r"$m_{HH}$ [GeV]", "ylabel": "Events that pass NGT Scouting", "bins": 30, "range": (50, 800), "color": "purple", "histtype": "fill", "hatch": ""},
        # "mHH_NGT":      {"label": r"Reco. $HH \to bb \tau \tau$",  "xlabel": r"$m_{HH}$ [GeV]", "ylabel": "Events that pass NGT Scouting", "bins": 30, "range": (50, 800), "color": "mediumvioletred", "histtype": "fill", "hatch": ""},

        # "mHb_all":      {"label": r"All events", "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 150), "color": "lightcyan", "histtype": "fill", "hatch": ""},
        # "mHb_HLT":      {"label": r"Pass HLT",   "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 150), "color": "deepskyblue", "histtype": "step", "hatch": "//"},
        # "mHb_NGT":      {"label": r"Pass NGT",   "xlabel": r"$m_{H}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 150), "color": "royalblue", "histtype": "step", "hatch": "/"},

        # "mHH_all":      {"label": r"All events", "xlabel": r"$m_{HH}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 200), "color": "lightcyan", "histtype": "fill", "hatch": ""},
        # "mHH_HLT":      {"label": r"Pass HLT",   "xlabel": r"$m_{HH}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 200), "color": "deepskyblue", "histtype": "step", "hatch": "//"},
        # "mHH_NGT":      {"label": r"Pass NGT",   "xlabel": r"$m_{HH}$ [GeV]", "ylabel": "Events", "bins": 50, "range": (0, 200), "color": "royalblue", "histtype": "step", "hatch": "/"},
        #}

    # separate_files = False   # True: one file per histogram. False: overlay all in plot_configs onto one figure.

    # if separate_files:
    #     for col, cfg in plot_configs.items():
    #         fig, ax = plt.subplots(figsize=(10, 8))
    #         make_panel(ax, col, cfg)
    #         plt.tight_layout()
    #         plt.savefig(f"plots/mHH/muon_channel/{col}.png", dpi=300, bbox_inches="tight")
    #         plt.close()
    # else:
    #     fig, ax = plt.subplots(figsize=(10, 8))
    #     for col, cfg in plot_configs.items():
    #         make_panel(ax, col, cfg)
    #     plt.tight_layout()
    #     plt.yscale('linear')
    #     plt.savefig("plots/mHH/mH_GEN_b_vs_jet.png", dpi=300, bbox_inches="tight")
    #     plt.close()

    # separate_files = False      #True gets 1 file with all the histograms, false gets a separate file for each histogram


    # data = {      
    # }

    # plot_configs = {
    #     "e_h":  {"label": r"$\tau \tau \to \tau_e \tau_h$",   "bins": 50, "range": (100, 150), "color": "orchid", "histtype": "fill"},
    #     "mu_h": {"label": r"$\tau \tau \to \tau_{\mu} \tau_h$","bins": 50, "range": (100, 150), "color": "darkorchid", "histtype": "fill"},
    #     "h_h":  {"label": r"$\tau \tau \to \tau_h \tau_h$",  "bins": 50, "range": (100, 150), "color": "mediumvioletred", "histtype": "fill"},
        
    #     "h_h_mHH":  {"label": r"All events",  "bins": 50, "range": (100, 150), "color": "orchid", "histtype": "fill"},
    #     "h_h_mHH_HLT":  {"label": r"Events that pass trigger",  "bins": 50, "range": (100, 150), "color": "mediumvioletred", "histtype": "fill"},

    #     "hadr_pt_leading":      {"label": r"All events",            "xlabel" : r"$p_T \tau_1$", "ylabel": "Events", "bins": 30, "range": (0, 300), "color": "palevioletred", "histtype": "step", "hatch": ""},
    #     "hadr_pt_leading_HLT":  {"label": r"Events that pass HLT",  "xlabel" : r"$p_T \tau_1$", "ylabel": "Events", "bins": 30, "range": (0, 300), "color": "purple", "histtype": "fill", "hatch": "//"},
    #     "hadr_pt_leading_NGT":  {"label": r"Events that pass NGT",  "xlabel" : r"$p_T \tau_1$", "ylabel": "Events", "bins": 30, "range": (0, 300), "color": "magenta", "histtype": "step", "hatch": ""},
        
    #     "hadr_pt_subleading":      {"label": r"All events",            "xlabel" : r"$p_T \tau_2$", "ylabel": "Events", "bins": 30, "range": (0, 300), "color": "violet", "histtype": "step", "hatch": ""},
    #     "hadr_pt_subleading_HLT":  {"label": r"Events that pass HLT",  "xlabel" : r"$p_T \tau_2$", "ylabel": "Events", "bins": 30, "range": (0, 300), "color": "orchid", "histtype": "fill", "hatch": ""},
    #     "hadr_pt_subleading_NGT":  {"label": r"Events that pass NGT",  "xlabel" : r"$p_T \tau_2$", "ylabel": "Events", "bins": 30, "range": (0, 300), "color": "rebeccapurple", "histtype": "step", "hatch": "//"},

    #     "hadr_eta1":      {"label": r"All events",            "xlabel" : r"$\eta \tau_1$", "ylabel": "Events", "bins": 30, "range": (-3, 3), "color": "olivedrab", "histtype": "step", "hatch": ""},
    #     "hadr_eta1_HLT":  {"label": r"Events that pass HLT",  "xlabel" : r"$\eta \tau_1$", "ylabel": "Events", "bins": 30, "range": (-3, 3), "color": "yellowgreen", "histtype": "fill", "hatch": ""},
    #     "hadr_eta1_NGT":  {"label": r"Events that pass NGT",  "xlabel" : r"$\eta \tau_1$", "ylabel": "Events", "bins": 30, "range": (-3, 3), "color": "yellowgreen", "histtype": "step", "hatch": "//"},

    #     "hadr_pt":      {"label": r"All events",              "xlabel" : r"$p_T$ [GeV]", "ylabel": "Events", "bins": 30, "range": (0, 300), "color": "palevioletred", "histtype": "step", "hatch": ""},

    #     "hadr_eta2":      {"label": r"All events",            "xlabel" : r"$\eta \tau_2$", "ylabel": "Events", "bins": 30, "range": (-3, 3), "color": "greenyellow", "histtype": "step", "hatch": ""},
    #     "hadr_eta2_HLT":  {"label": r"Events that pass HLT",  "xlabel" : r"$\eta \tau_2$", "ylabel": "Events", "bins": 30, "range": (-3, 3), "color": "forestgreen",  "histtype": "fill", "hatch": ""},
    #     "hadr_eta2_NGT":  {"label": r"Events that pass NGT",  "xlabel" : r"$\eta \tau_2$", "ylabel": "Events", "bins": 30, "range": (-3, 3), "color": "forestgreen", "histtype": "step", "hatch": "//"},
    # }

    # def make_panel(ax, col, cfg):
    #     counts, edges = np.histogram(data[col], bins=cfg["bins"], range=cfg["range"])
    #     hep.histplot(counts, edges, ax=ax, label=cfg["label"], color=cfg["color"], hatch = cfg["hatch"], histtype=cfg["histtype"])
    #     ax.set_xlabel(cfg["xlabel"])
    #     ax.set_ylabel(cfg["ylabel"])
    #     ax.legend(loc="upper right")
    #     hep.cms.label("Preliminary", data=False, ax=ax, com=14)
    #     return counts, edges

    # eff_color = "orchid"
    # efficiency_counts = []
    # if separate_files:
    #     for col, cfg in plot_configs.items():
    #         fig, ax = plt.subplots(figsize=(10, 8))
    #         make_panel(ax, col, cfg)
    #         plt.tight_layout()
    #         plt.savefig(f"plots/{col}.png", dpi=300, bbox_inches="tight")
    #         plt.close()
    # else:
    #     fig, ax = plt.subplots(figsize=(10, 8))
    #     for col, cfg in plot_configs.items():
    #         counts, edges = make_panel(ax, col, cfg)
    #         efficiency_counts.append(counts)
            
    #     Efficiency_HLT = [num/den if den != 0 else 0 for num,den in zip(efficiency_counts[1], efficiency_counts[0])]
    #     Efficiency_NGT = [num/den if den != 0 else 0 for num,den in zip(efficiency_counts[2], efficiency_counts[0])]
    #     centers = 0.5*(edges[:-1]+edges[1:])
    #     bin_widths = np.diff(edges)
    #     y_error_HLT = [np.sqrt(eff * (1 - eff) / den) if den > 0 else 0 for eff, den in zip(Efficiency_HLT, efficiency_counts[0])] 
    #     y_error_NGT = [np.sqrt(eff * (1 - eff) / den) if den > 0 else 0 for eff, den in zip(Efficiency_NGT, efficiency_counts[0])] 
    #     ax2 = ax.twinx()

    #     ax2.errorbar(centers, Efficiency_HLT, xerr=0.5*bin_widths, yerr = y_error_HLT, color=eff_color, fmt = 'o')
    #     ax2.errorbar(centers, Efficiency_NGT, xerr=0.5*bin_widths, yerr = y_error_NGT, color='purple', fmt = '*')
    #     ax2.set_ylabel("Efficiency", color=eff_color)
    #     ax2.tick_params(axis='y', labelcolor=eff_color)
    #     ax2.set_ylim(0,1)
    #     hep.cms.label("Preliminary", data=False, ax=ax, com=14)
    #     plt.tight_layout()
    #     ax2.legend(loc='upper right')
    #     plt.savefig("plots/eta/eta2+eff.png", dpi=300, bbox_inches="tight")
    #     plt.close()


    # #─────────────────────────────────────── Define and plot histograms 2D ────────────────────────────────────────

    # separate_files_2D = True

    # data_2D = {

    # }

    # plot_configs_2D = {
    #     "pTLeading_Vs_Subleading_all": {"x": "tau_pt_lead", "y": "tau_pt_sublead", "bins": 50, "range": [[0, 300], [0, 300]], "cmap": "YlGn", "cbar_label": "All events"},
    #     "pTLeading_Vs_Subleading_HLT": {"x": "tau_pt_lead", "y": "tau_pt_sublead", "bins": 50, "range": [[0, 300], [0, 300]], "cmap": "YlGn", "cbar_label": "Events that pass HLT"},
    #     "pTLeading_Vs_Subleading_NGT": {"x": "tau_pt_lead", "y": "tau_pt_sublead", "bins": 50, "range": [[0, 300], [0, 300]], "cmap": "YlGn", "cbar_label": "Events that pass NGT"},
        
    #     "pT_Efficiency_HLT": {"x": "tau_pt_lead", "y": "tau_pt_sublead", "bins": 50, "range": [[0, 300], [0, 300]], "cmap": "YlGn", "cbar_label": "Efficiency HLT"},
    #     "pT_Efficiency_NGT": {"x": "tau_pt_lead", "y": "tau_pt_sublead", "bins": 50, "range": [[0, 300], [0, 300]], "cmap": "YlGn", "cbar_label": "Efficiency NGT"},
    # }

    # bins = plot_configs_2D["eta_Efficiency_HLT"]["bins"]
    # range_ = plot_configs_2D["eta_Efficiency_HLT"]["range"]

    # hist_all, xedges, yedges = np.histogram2d(
    #     data_2D["eta1_Vs_eta2_all"]["eta1"],
    #     data_2D["eta1_Vs_eta2_all"]["eta2"],
    #     bins=bins, range=range_
    # )
    # hist_HLT, xedges, yedges = np.histogram2d(
    #     data_2D["eta1_Vs_eta2_HLT"]["eta1"],
    #     data_2D["eta1_Vs_eta2_HLT"]["eta2"],
    #     bins=bins, range=range_
    # )
    # hist_ngt, xedges, yedges = np.histogram2d(
    #     data_2D["eta1_Vs_eta2_NGT"]["eta1"],
    #     data_2D["eta1_Vs_eta2_NGT"]["eta2"],
    #     bins=bins, range=range_
    # )
    # hist_efficiency_HLT = np.divide(hist_HLT, hist_all,
    #     out=np.zeros_like(hist_HLT, dtype=float),
    #     where=hist_all > 0
    # )
    # hist_efficiency_NGT = np.divide(hist_ngt, hist_all,
    #     out=np.zeros_like(hist_ngt, dtype=float),
    #     where=hist_all > 0
    # )
    # data_2D["eta_Efficiency_HLT"] = (hist_efficiency_HLT, xedges, yedges) #Appends Efficiency tuple in data_2D dict
    # data_2D["eta_Efficiency_NGT"] = (hist_efficiency_NGT, xedges, yedges) #Appends Efficiency NGT tuple in data_2D dict

    # def make_panel_2d(ax, col, cfg):
    #     if isinstance(data_2D[col], tuple): #For Efficiency histogram, because that one is a tuple in the dictionnary
    #         h, xedges, yedges = data_2D[col]
    #     else:                               #The other histograms are stored as dictionnaries in data, not tuples, so different treatement
    #         h, xedges, yedges = np.histogram2d(
    #             data_2D[col][cfg["x"]], data_2D[col][cfg["y"]],
    #             bins=cfg["bins"], range=cfg["range"]
    #         )

    #     mesh = ax.pcolormesh(xedges, yedges, h.T, cmap=cfg["cmap"])
    #     plt.colorbar(mesh, ax=ax, label=cfg["cbar_label"])
    #     ax.set_xlabel(r"$\eta^{\mathrm{leading}}$")
    #     ax.set_ylabel(r"$\eta^{\mathrm{subleading}}$")
    #     hep.cms.label("Preliminary", data=False, ax=ax, com=14)

    # if separate_files_2D:
    #     for col, cfg in plot_configs_2D.items():
    #         fig, ax = plt.subplots(figsize=(10, 8))
    #         make_panel_2d(ax, col, cfg)
    #         plt.tight_layout()
    #         plt.savefig(f"plots/eta/2D_{col}.png", dpi=300, bbox_inches="tight")
    #         plt.close()
    # else:
    #     fig, axes = plt.subplots(1, 4, figsize=(30, 8))
    #     for ax, (col, cfg) in zip(axes, plot_configs_2D.items()):
    #         make_panel_2d(ax, col, cfg)
    #     plt.tight_layout()
    #     plt.savefig("plots/eta/2D_eta_All+Triggered+EffHLT+EFFNGT.png", dpi=300, bbox_inches="tight")
    #     plt.close()

    # import pdb
    # pdb.set_trace()
    # all_events = df_tau_h.Histo2D(("h_all_events", "", 50, 0, 300, 50, 0, 300), "tau_pt_lead", "tau_pt_sublead")
    # all_events_HLT = df_tau_h_HLT.Histo2D(("h_all_events_HLT", "", 50, 0, 300, 50, 0, 300), "tau_pt_lead", "tau_pt_sublead")

    # h_tot = all_events.GetValue()
    # h_pas = all_events_HLT.GetValue()

    # h_eff = h_pas.Clone("h_eff")
    # h_eff.Reset("ICES")

    # h_tot.SetTitle(r"Events;$p_T\tau_{1}$;$p_T\tau_{2}$")
    # h_pas.SetTitle(r"Events;$p_T\tau_{1}$;$p_T\tau_{2}$")
    # h_eff.SetTitle(r"Efficiency;$p_T\tau_{1}$;$p_T\tau_{2}$")

    # h_eff.Divide(h_pas, h_tot, 1.0, 1.0, "B")

    # for h, name in zip([h_tot, h_pas, h_eff], ["2D_pTLeading_Vs_Subleading_total", "2D_pTLeading_Vs_Subleading_HLT", "2D_Efficiency"]):
    #     xedges = np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1, h.GetNbinsX() + 2)])
    #     yedges = np.array([h.GetYaxis().GetBinLowEdge(i) for i in range(1, h.GetNbinsY() + 2)])

    #     #Bin contents
    #     z = np.zeros((h.GetNbinsY(), h.GetNbinsX()))

    #     for ix in range(1, h.GetNbinsX() + 1):
    #         for iy in range(1, h.GetNbinsY() + 1):
    #             z[iy-1, ix-1] = h.GetBinContent(ix, iy)

    #     #Plot
    #     plt.figure(figsize=(10, 8))
    #     plt.pcolormesh(xedges, yedges, z, cmap='RdPu')
    #     plt.colorbar(label="Efficiency")
    #     hep.cms.label("Preliminary", data=False, com=14)
    #     plt.xlabel(h.GetXaxis().GetTitle())
    #     plt.ylabel(h.GetYaxis().GetTitle())
    #     plt.tight_layout()
    #     plt.savefig(f"/eos/user/s/sbenabde/CERN_Summer_student/plots/{name}")
#     #─────────────────────────────────────── Save root file output ────────────────────────────────────────

#     # HH_mass = df.Histo1D(("HH_mass", "HH_mass", 100, -10, 200), "HH_mass")

#     # out = ROOT.TFile("Signal_test.root", "RECREATE")
#     # for h in [HH_mass]:
#     #     h.Write()
#     # out.Close()
#     # print("\nSaved Signal_test.root")

