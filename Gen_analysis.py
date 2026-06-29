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
        df.Snapshot("Events", "Dataframes/df_ngt_hlt.root", snapshot_columns, snapshot_options)
    else:
        df = ROOT.RDataFrame("Events", "Dataframes/df_ngt_hlt.root")

    #─────────────────────────────────────── Define variables ────────────────────────────────────────

    df = (
        df
        .Define("H_tautau_idx",          f"Get_H_idx(15, GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)") #Index Higgs from tau
        .Define("H_bb_idx",              f"Get_H_idx(5, GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")  #Index Higgs from bb

        .Define("H_tautau_mass",         f"GenPart_mass.at(H_tautau_idx)")
        .Define("H_bb_mass",             f"GenPart_mass.at(H_bb_idx)")
        .Define("HH_mass",               f"Get_mHH(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, H_tautau_idx, H_bb_idx)")

        .Define("tau_indices",           "Get_tau_from_H_indices(GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")
        .Define("tau_channel",           "Get_tau_channel(tau_indices, GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")

        .Define("tau_pt",                "Get_tau_pt(tau_indices, GenPart_pt)")
        .Define("tau_pt_lead",           "tau_pt[0]")
        .Define("tau_pt_sublead",        "tau_pt[1]")

        .Define("tau_eta",               "Get_tau_eta(tau_indices, GenPart_pt, GenPart_eta)")
        .Define("eta1",                  "tau_eta[0]")
        .Define("eta2",                  "tau_eta[1]")
    )
    #─────────────────────────────────────── Print column names ────────────────────────────────────────

    all_cols = [str(c) for c in df.GetColumnNames()]
    L1_cols = [c for c in all_cols if "L1_p" in c]
    print("\nGenPart columns found:")
    for c in sorted(L1_cols):
        print(f"  {c}  [{df.GetColumnType(c)}]")

        df_L1 = df.Filter("tau_channel  == 2")
        df_L1_filtered = df_L1.Filter(f"{c} == true")
        
        efficiency_L1  = overall_efficiency(df_L1, df_L1_filtered)
        print(efficiency_L1)

    #─────────────────────────────────────── Filter dataframe ────────────────────────────────────────

    df_tau_e = df.Filter("tau_channel  == 0")
    df_tau_e_HLT = df_tau_e.Filter("(HLT_Ele32_WPTight_L1Seeded == true) && (HLT_Ele30_WPTight_L1Seeded_LooseDeepTauPFTauHPS30_eta2p1_CrossL1 == true)")
    df_tau_e_NGT = df_tau_e.Filter("DST_PFScouting == true")

    df_tau_m = df.Filter("tau_channel == 1")
    df_tau_m_HLT = df_tau_m.Filter("(HLT_IsoMu24_FromL1TkMuon == true) && (HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1 == true)")
    df_tau_m_NGT = df_tau_m.Filter("DST_PFScouting == true")

    df_tau_h = df.Filter("tau_channel  == 2 ")
    df_tau_h_HLT = df_tau_h.Filter("HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1 == true")
    df_tau_h_NGT = df_tau_h.Filter("DST_PFScouting == true")

    def overall_efficiency(all, trigger):
        all_c = all.Count().GetValue()
        trigger_c = trigger.Count().GetValue()
        eff = trigger_c/all_c
        return eff 

    efficiency_DiTau35 = overall_efficiency(df_tau_h, df_tau_h_HLT)
    efficiency_NGT     = overall_efficiency(df_tau_h, df_tau_h_NGT)
    print(efficiency_DiTau35)
    print(efficiency_NGT)
    #─────────────────────────────────────── Define and plot histograms 1D ────────────────────────────────────────

    separate_files = False      #True gets 1 file with all the histograms, false gets a separate file for each histogram


    data = {
        "e_h":  df_tau_e.AsNumpy(["HH_mass"])["HH_mass"],
        "mu_h": df_tau_m.AsNumpy(["HH_mass"])["HH_mass"],
        "h_h":  df_tau_h.AsNumpy(["HH_mass"])["HH_mass"],

        "h_h_mHH":      df_tau_h.AsNumpy(["HH_mass"])["HH_mass"],
        "h_h_mHH_HLT":  df_tau_h_HLT.AsNumpy(["HH_mass"])["HH_mass"],

        "hadr_pt_leading":      df_tau_h.AsNumpy(["tau_pt_lead"])["tau_pt_lead"],
        "hadr_pt_leading_HLT":  df_tau_h_HLT.AsNumpy(["tau_pt_lead"])["tau_pt_lead"],
        "hadr_pt_leading_NGT":  df_tau_h_NGT.AsNumpy(["tau_pt_lead"])["tau_pt_lead"],

        "hadr_pt_subleading":      df_tau_h.AsNumpy(["tau_pt_sublead"])["tau_pt_sublead"],
        "hadr_pt_subleading_HLT":  df_tau_h_HLT.AsNumpy(["tau_pt_sublead"])["tau_pt_sublead"],
        "hadr_pt_subleading_NGT":  df_tau_h_NGT.AsNumpy(["tau_pt_sublead"])["tau_pt_sublead"],

        "hadr_eta1":      df_tau_h.AsNumpy(["eta1"])["eta1"],
        "hadr_eta1_HLT":  df_tau_h_HLT.AsNumpy(["eta1"])["eta1"],
        "hadr_eta1_NGT":  df_tau_h_NGT.AsNumpy(["eta1"])["eta1"],

        "hadr_eta2":      df_tau_h.AsNumpy(["eta2"])["eta2"],
        "hadr_eta2_HLT":  df_tau_h_HLT.AsNumpy(["eta2"])["eta2"],
        "hadr_eta2_NGT":  df_tau_h_NGT.AsNumpy(["eta2"])["eta2"],

    }

    plot_configs = {
        # "e_h":  {"label": r"$\tau \tau \to \tau_e \tau_h$",   "bins": 50, "range": (200, 800), "color": "orchid", "histtype": "fill"},
        # "mu_h": {"label": r"$\tau \tau \to \tau_{\mu} \tau_h$","bins": 50, "range": (200, 800), "color": "darkorchid", "histtype": "fill"},
        # "h_h":  {"label": r"$\tau \tau \to \tau_h \tau_h$",  "bins": 50, "range": (200, 800), "color": "mediumvioletred", "histtype": "fill"},
        
        # "h_h_mHH":  {"label": r"All events",  "bins": 50, "range": (200, 800), "color": "orchid", "histtype": "fill"},
        # "h_h_mHH_HLT":  {"label": r"Events that pass trigger",  "bins": 50, "range": (200, 800), "color": "mediumvioletred", "histtype": "fill"},

        # "hadr_pt_leading":      {"label": r"All events",            "xlabel" : r"$p_T \tau_1$", "ylabel": "Events", "bins": 30, "range": (0, 300), "color": "palevioletred", "histtype": "step", "hatch": ""},
        # "hadr_pt_leading_HLT":  {"label": r"Events that pass HLT",  "xlabel" : r"$p_T \tau_1$", "ylabel": "Events", "bins": 30, "range": (0, 300), "color": "purple", "histtype": "fill", "hatch": "//"},
        # "hadr_pt_leading_NGT":  {"label": r"Events that pass NGT",  "xlabel" : r"$p_T \tau_1$", "ylabel": "Events", "bins": 30, "range": (0, 300), "color": "magenta", "histtype": "step", "hatch": ""},
        
        # "hadr_pt_subleading":      {"label": r"All events",            "xlabel" : r"$p_T \tau_2$", "ylabel": "Events", "bins": 30, "range": (0, 300), "color": "violet", "histtype": "step", "hatch": ""},
        # "hadr_pt_subleading_HLT":  {"label": r"Events that pass HLT",  "xlabel" : r"$p_T \tau_2$", "ylabel": "Events", "bins": 30, "range": (0, 300), "color": "orchid", "histtype": "fill", "hatch": ""},
        # "hadr_pt_subleading_NGT":  {"label": r"Events that pass NGT",  "xlabel" : r"$p_T \tau_2$", "ylabel": "Events", "bins": 30, "range": (0, 300), "color": "rebeccapurple", "histtype": "step", "hatch": "//"},

        # "hadr_eta1":      {"label": r"All events",            "xlabel" : r"$\eta \tau_1$", "ylabel": "Events", "bins": 30, "range": (-3, 3), "color": "olivedrab", "histtype": "step", "hatch": ""},
        # "hadr_eta1_HLT":  {"label": r"Events that pass HLT",  "xlabel" : r"$\eta \tau_1$", "ylabel": "Events", "bins": 30, "range": (-3, 3), "color": "yellowgreen", "histtype": "fill", "hatch": ""},
        # "hadr_eta1_NGT":  {"label": r"Events that pass NGT",  "xlabel" : r"$\eta \tau_1$", "ylabel": "Events", "bins": 30, "range": (-3, 3), "color": "yellowgreen", "histtype": "step", "hatch": "//"},

        "hadr_eta2":      {"label": r"All events",            "xlabel" : r"$\eta \tau_2$", "ylabel": "Events", "bins": 30, "range": (-3, 3), "color": "greenyellow", "histtype": "step", "hatch": ""},
        "hadr_eta2_HLT":  {"label": r"Events that pass HLT",  "xlabel" : r"$\eta \tau_2$", "ylabel": "Events", "bins": 30, "range": (-3, 3), "color": "forestgreen",  "histtype": "fill", "hatch": ""},
        "hadr_eta2_NGT":  {"label": r"Events that pass NGT",  "xlabel" : r"$\eta \tau_2$", "ylabel": "Events", "bins": 30, "range": (-3, 3), "color": "forestgreen", "histtype": "step", "hatch": "//"},
    }

    def make_panel(ax, col, cfg):
        counts, edges = np.histogram(data[col], bins=cfg["bins"], range=cfg["range"])
        hep.histplot(counts, edges, ax=ax, label=cfg["label"], color=cfg["color"], hatch = cfg["hatch"], histtype=cfg["histtype"])
        ax.set_xlabel(cfg["xlabel"])
        ax.set_ylabel(cfg["ylabel"])
        ax.legend(loc="upper right")
        hep.cms.label("Preliminary", data=False, ax=ax, com=14)
        return counts, edges

    eff_color = "orchid"
    efficiency_counts = []
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
            counts, edges = make_panel(ax, col, cfg)
            efficiency_counts.append(counts)
            
        Efficiency_HLT = [num/den if den != 0 else 0 for num,den in zip(efficiency_counts[1], efficiency_counts[0])]
        Efficiency_NGT = [num/den if den != 0 else 0 for num,den in zip(efficiency_counts[2], efficiency_counts[0])]
        centers = 0.5*(edges[:-1]+edges[1:])
        bin_widths = np.diff(edges)
        y_error_HLT = [np.sqrt(eff * (1 - eff) / den) if den > 0 else 0 for eff, den in zip(Efficiency_HLT, efficiency_counts[0])] 
        y_error_NGT = [np.sqrt(eff * (1 - eff) / den) if den > 0 else 0 for eff, den in zip(Efficiency_NGT, efficiency_counts[0])] 
        ax2 = ax.twinx()

        ax2.errorbar(centers, Efficiency_HLT, xerr=0.5*bin_widths, yerr = y_error_HLT, color=eff_color, fmt = 'o')
        ax2.errorbar(centers, Efficiency_NGT, xerr=0.5*bin_widths, yerr = y_error_NGT, color='purple', fmt = '*')
        ax2.set_ylabel("Efficiency", color=eff_color)
        ax2.tick_params(axis='y', labelcolor=eff_color)
        ax2.set_ylim(0,1)
        hep.cms.label("Preliminary", data=False, ax=ax, com=14)
        plt.tight_layout()
        ax2.legend(loc='upper right')
        plt.savefig("plots/eta/eta2+eff.png", dpi=300, bbox_inches="tight")
        plt.close()


    #─────────────────────────────────────── Define and plot histograms 2D ────────────────────────────────────────

    separate_files_2D = True

    data_2D = {
        "pTLeading_Vs_Subleading_all": df_tau_h.AsNumpy(['tau_pt_lead', 'tau_pt_sublead']),
        "pTLeading_Vs_Subleading_HLT": df_tau_h_HLT.AsNumpy(['tau_pt_lead', 'tau_pt_sublead']),
        "pTLeading_Vs_Subleading_NGT": df_tau_h_NGT.AsNumpy(['tau_pt_lead', 'tau_pt_sublead']),
        
        "eta1_Vs_eta2_all": df_tau_h.AsNumpy(['eta1', 'eta2']),
        "eta1_Vs_eta2_HLT": df_tau_h_HLT.AsNumpy(['eta1', 'eta2']),
        "eta1_Vs_eta2_NGT": df_tau_h_NGT.AsNumpy(['eta1', 'eta2']),
    }

    plot_configs_2D = {
        # "pTLeading_Vs_Subleading_all": {"x": "tau_pt_lead", "y": "tau_pt_sublead", "bins": 50, "range": [[0, 300], [0, 300]], "cmap": "YlGn", "cbar_label": "All events"},
        # "pTLeading_Vs_Subleading_HLT": {"x": "tau_pt_lead", "y": "tau_pt_sublead", "bins": 50, "range": [[0, 300], [0, 300]], "cmap": "YlGn", "cbar_label": "Events that pass HLT"},
        # "pTLeading_Vs_Subleading_NGT": {"x": "tau_pt_lead", "y": "tau_pt_sublead", "bins": 50, "range": [[0, 300], [0, 300]], "cmap": "YlGn", "cbar_label": "Events that pass NGT"},
        
        # "pT_Efficiency_HLT": {"x": "tau_pt_lead", "y": "tau_pt_sublead", "bins": 50, "range": [[0, 300], [0, 300]], "cmap": "YlGn", "cbar_label": "Efficiency HLT"},
        # "pT_Efficiency_NGT": {"x": "tau_pt_lead", "y": "tau_pt_sublead", "bins": 50, "range": [[0, 300], [0, 300]], "cmap": "YlGn", "cbar_label": "Efficiency NGT"},

        "eta1_Vs_eta2_all":  {"x": "eta1", "y": "eta2", "bins": 50, "range": [[-4, 4], [-4, 4]], "cmap": "RdPu", "cbar_label": "All events"},
        "eta1_Vs_eta2_HLT":  {"x": "eta1", "y": "eta2", "bins": 50, "range": [[-4, 4], [-4, 4]], "cmap": "RdPu", "cbar_label": "Events that pass HLT"},
        "eta1_Vs_eta2_NGT":  {"x": "eta1", "y": "eta2", "bins": 50, "range": [[-4, 4], [-4, 4]], "cmap": "RdPu", "cbar_label": "Events that pass NGT"},
        
        "eta_Efficiency_HLT":  {"x": "eta1", "y": "eta2", "bins": 50, "range": [[-4, 4], [-4, 4]], "cmap": "RdPu", "cbar_label": "Efficiency HLT"},
        "eta_Efficiency_NGT":  {"x": "eta1", "y": "eta2", "bins": 50, "range": [[-4, 4], [-4, 4]], "cmap": "RdPu", "cbar_label": "Efficiency NGT"},
    }

    bins = plot_configs_2D["eta_Efficiency_HLT"]["bins"]
    range_ = plot_configs_2D["eta_Efficiency_HLT"]["range"]

    hist_all, xedges, yedges = np.histogram2d(
        data_2D["eta1_Vs_eta2_all"]["eta1"],
        data_2D["eta1_Vs_eta2_all"]["eta2"],
        bins=bins, range=range_
    )
    hist_HLT, xedges, yedges = np.histogram2d(
        data_2D["eta1_Vs_eta2_HLT"]["eta1"],
        data_2D["eta1_Vs_eta2_HLT"]["eta2"],
        bins=bins, range=range_
    )
    hist_ngt, xedges, yedges = np.histogram2d(
        data_2D["eta1_Vs_eta2_NGT"]["eta1"],
        data_2D["eta1_Vs_eta2_NGT"]["eta2"],
        bins=bins, range=range_
    )
    hist_efficiency_HLT = np.divide(hist_HLT, hist_all,
        out=np.zeros_like(hist_HLT, dtype=float),
        where=hist_all > 0
    )
    hist_efficiency_NGT = np.divide(hist_ngt, hist_all,
        out=np.zeros_like(hist_ngt, dtype=float),
        where=hist_all > 0
    )
    data_2D["eta_Efficiency_HLT"] = (hist_efficiency_HLT, xedges, yedges) #Appends Efficiency tuple in data_2D dict
    data_2D["eta_Efficiency_NGT"] = (hist_efficiency_NGT, xedges, yedges) #Appends Efficiency NGT tuple in data_2D dict

    def make_panel_2d(ax, col, cfg):
        if isinstance(data_2D[col], tuple): #For Efficiency histogram, because that one is a tuple in the dictionnary
            h, xedges, yedges = data_2D[col]
        else:                               #The other histograms are stored as dictionnaries in data, not tuples, so different treatement
            h, xedges, yedges = np.histogram2d(
                data_2D[col][cfg["x"]], data_2D[col][cfg["y"]],
                bins=cfg["bins"], range=cfg["range"]
            )

        mesh = ax.pcolormesh(xedges, yedges, h.T, cmap=cfg["cmap"])
        plt.colorbar(mesh, ax=ax, label=cfg["cbar_label"])
        ax.set_xlabel(r"$\eta^{\mathrm{leading}}$")
        ax.set_ylabel(r"$\eta^{\mathrm{subleading}}$")
        hep.cms.label("Preliminary", data=False, ax=ax, com=14)

    if separate_files_2D:
        for col, cfg in plot_configs_2D.items():
            fig, ax = plt.subplots(figsize=(10, 8))
            make_panel_2d(ax, col, cfg)
            plt.tight_layout()
            plt.savefig(f"plots/eta/2D_{col}.png", dpi=300, bbox_inches="tight")
            plt.close()
    else:
        fig, axes = plt.subplots(1, 4, figsize=(30, 8))
        for ax, (col, cfg) in zip(axes, plot_configs_2D.items()):
            make_panel_2d(ax, col, cfg)
        plt.tight_layout()
        plt.savefig("plots/eta/2D_eta_All+Triggered+EffHLT+EFFNGT.png", dpi=300, bbox_inches="tight")
        plt.close()

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
    #─────────────────────────────────────── Save root file output ────────────────────────────────────────

    # HH_mass = df.Histo1D(("HH_mass", "HH_mass", 100, -10, 200), "HH_mass")

    # out = ROOT.TFile("Signal_test.root", "RECREATE")
    # for h in [HH_mass]:
    #     h.Write()
    # out.Close()
    # print("\nSaved Signal_test.root")

