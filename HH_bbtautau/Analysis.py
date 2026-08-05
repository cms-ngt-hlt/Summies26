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

#-------------------------Functions------------------------------------

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

#─────────────────────────────────────── Get dataframe ────────────────────────────────────────

df = ROOT.RDataFrame("Events", "/eos/user/s/sbenabde/CERN_Summer_student/Dataframes/df_ngt_hlt.root")

    #─────────────────────────────────────── Define variables ────────────────────────────────────────

df = (
    df
    #-------------- Gen -------------------

    .Define("Gen_b_idx",             f"Get_b_from_H_indices(GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")
    .Define("Gen_bJets_idx",         f"Match_b_to_GenJet(Gen_b_idx, GenPart_eta, GenPart_phi, GenJet_eta, GenJet_phi)")

    .Define("Gen_b_p4",              f"Get_b_p4(Gen_b_idx, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
    .Define("Gen_bjet_p4",           f"Get_b_p4(Gen_bJets_idx, GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass)")

    .Define("Gen_bjet_pt",           f"Get_p4_pt(Gen_bjet_p4)")
    .Define("Gen_bjet_eta",          f"Get_p4_eta(Gen_bjet_p4)")
    .Define("Gen_bjet_phi",          f"Get_p4_phi(Gen_bjet_p4)")
    .Define("Gen_bjet_mass",         f"Get_p4_mass(Gen_bjet_p4)")

    .Define("Gen_H_bb_mass",         f"Get_Gen_H_tautau_mass(Gen_b_p4)") 
    .Define("Gen_H_bbjet_mass",      f"Get_Gen_H_tautau_mass(Gen_bjet_p4)") 

    .Define("Gen_tau_idx",           f"Get_tau_from_H_indices(GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")
    .Define("tau_channel",           f"Get_tau_channel(Gen_tau_idx, GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")

    .Define("Gen_tau_p4",            f"Get_visible_tau_p4s(Gen_tau_idx, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
    .Define("Gen_tau_pt",            f"Get_p4_pt(Gen_tau_p4)")
    .Define("Gen_tau_eta",           f"Get_p4_eta(Gen_tau_p4)")
    .Define("Gen_tau_phi",           f"Get_p4_phi(Gen_tau_p4)")
    .Define("Gen_tau_mass",          f"Get_p4_mass(Gen_tau_p4)")

    .Define("Gen_Muon_idx",          f"Get_lepton_idx(13, GenPart_pdgId, Gen_tau_idx, GenPart_statusFlags, GenPart_genPartIdxMother)")
    .Define("nGen_Muon",             f"Gen_Muon_idx.size()")

    .Define("Gen_Electron_idx",      f"Get_lepton_idx(11, GenPart_pdgId, Gen_tau_idx, GenPart_statusFlags, GenPart_genPartIdxMother)")
    .Define("nGen_Electron",         f"Gen_Electron_idx.size()")


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
df_NGT = df.Filter("tau_channel >= 0 && DST_PFScouting == true")
 
# One dict keyed by channel index holds everything that used to live in
# three separate dicts (channel filter, suffix, and title), so they can
# never drift out of sync with each other.
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

ReconstructedTau =  {idx: channel[0].Filter("DST_PFScouting == true").Filter("Matched_tau_idx[0] >= 0 && Matched_tau_idx[1] >= 0") for idx, cfg in channel.items()} 
ReconstructedJet =  {idx: channel[1].Filter("DST_PFScouting == true").Filter("Matched_jet_idx[0] >= 0 && Matched_jet_idx[1] >= 0") for idx, cfg in channel.items()} 
ReconstructedBoth = {idx: channel[2].Filter("DST_PFScouting == true").Filter("Matched_tau_idx[0] >= 0 && Matched_tau_idx[1] >= 0 && Matched_jet_idx[0] >= 0 && Matched_jet_idx[1] >= 0") for idx, cfg in channel.items()} 
# ─────────────────────────────────────── Define and plot histograms 1D ────────────────────────────────────────    
 
var_configs = {
    "Gen_tau_pt":     {"branch": "Gen_tau_pt",   "xlabel": r"$p_{T}(\tau)$ [GeV]", "bins": 30, "range": (0, 300)},
    "Gen_tau_eta":    {"branch": "Gen_tau_eta",  "xlabel": r"$\eta(\tau)$",        "bins": 30, "range": (-3, 3)},
    "bscore":         {"branch": "Matched_jet_bscore", "xlabel": r"$b-score$",     "bins": 30, "range": (0, 1)},

    "nGen_Muon":       {"branch": "nGen_Muon",     "xlabel": r"n Gen $\mu$",       "bins": 4, "range": (-0.5, 3.5)},
    "nGen_Electron":   {"branch": "nGen_Electron", "xlabel": r"n Gen $e$",         "bins": 4, "range": (-0.5, 3.5)},
    
    "Gen_bjet_pt":     {"branch": "Gen_bjet_pt",   "xlabel": r"$p_{T}$ Jet [GeV]", "bins": 30, "range": (0, 300)},
    "Gen_bjet_eta":    {"branch": "Gen_bjet_eta",  "xlabel": r"$\eta$ Jet",        "bins": 30, "range": (-3, 3)},
    
    "Reco_jet_pT":     {"branch": "Reco_jet_pT",   "xlabel": r"$p_{T}$ Jet [GeV]", "bins": 30, "range": (0, 300)},
    "Reco_jet_eta":    {"branch": "Reco_jet_eta",  "xlabel": r"$\eta$ Jet",        "bins": 30, "range": (-3, 3)},
}

selection_style = {
    "all": {"label": "All events",           "color": "mediumvioletred",  "histtype": "step", "hatch": ""},
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
    recoTau_arrs  = get_leg_arrays(ReconstructedTau[ch],  branches)
    recoJet_arrs  = get_leg_arrays(ReconstructedJet[ch],  branches)
    recoBoth_arrs = get_leg_arrays(ReconstructedBoth[ch], branches)

    for var, vcfg in var_configs.items():
        data[f"{var}_all_{suf}"] = all_arrs[vcfg["branch"]]
        data[f"{var}_HLT_{suf}"] = hlt_arrs[vcfg["branch"]]
        data[f"{var}_NGT_{suf}"] = ngt_arrs[vcfg["branch"]]
        data[f"{var}_RecoTau_{suf}"]  = recoTau_arrs[vcfg["branch"]]
        data[f"{var}_RecoJet_{suf}"]  = recoJet_arrs[vcfg["branch"]]
        data[f"{var}_RecoBoth_{suf}"] = recoBoth_arrs[vcfg["branch"]]

 
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
    ax2.legend(loc="upper left", fontsize=20)
    ax.legend(loc="upper right", fontsize=20)
    plt.tight_layout()
 
    if scale == 'log':
        ax.set_yscale('log')
 
    outname = f"{path}/eff_{var}_{suf}_{scale}.{save_as}"
    fig.savefig(outname)
    plt.close(fig)
    print(f"Saved {outname}")
 

def plot_reco_efficiency(var, ch, save_as, scale, tau=False, jet=False, both=False):
    path = "/eos/user/s/sbenabde/CERN_Summer_student/HH_bbtautau/plots/Reco_Efficiency"
    os.makedirs(path, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 8))
    suf   = channels[ch]["suffix"]
    vcfg  = var_configs[var]
    bins, rng = vcfg["bins"], vcfg["range"]

    hist_data = {}

    arr_all  = data[f"{var}_all_{suf}"]
    if tau==True:
        arr_reco = data[f"{var}_RecoTau_{suf}"]
        lb= r"2 matched $\tau$s"

    if jet==True:
        arr_reco = data[f"{var}_RecoJet_{suf}"]
        lb= r"2 matched $b$-jets"
    if both==True:
        arr_reco = data[f"{var}_RecoBoth_{suf}"]
        lb= r"2 matched $\tau$'s + $b$-jets"

    counts, edges = np.histogram(arr_all, bins=bins, range=rng, density=False)
    counts2, edges2 = np.histogram(arr_reco, bins=bins, range=rng, density=False)
    hist_data['all']  = counts
    hist_data['Reco'] = counts2
    hep.histplot(counts, edges, ax=ax,   label="All", color="mediumvioletred", histtype="step")
    hep.histplot(counts2, edges2, ax=ax, label=lb,    color="purple",          histtype="fill")

    ax.set_xlabel(vcfg["xlabel"])
    ax.set_ylabel("Events")
    hep.cms.label("Preliminary", data=False, ax=ax, com=14)
    ax.text(0.125, -0.08, channels[ch]["title"], transform=ax.transAxes, fontsize=20, ha="right", va="top")

    centers     = 0.5 * (edges[:-1] + edges[1:])
    bin_widths  = np.diff(edges)
    den         = hist_data["all"]

    eff_NGT, err_NGT = eff_and_err(hist_data["Reco"], den)

    eff_color = "olivedrab"
    ax2 = ax.twinx()
    ax2.errorbar(centers, eff_NGT, xerr=0.5 * bin_widths, yerr=err_NGT, color=eff_color, fmt='*', label="Reco efficiency")
    ax2.set_ylabel("Efficiency", color=eff_color)
    ax2.tick_params(axis='y', labelcolor=eff_color)
    ax2.set_ylim(0, 1)
    ax2.legend(loc="upper left", fontsize=20)
    ax.legend(loc="upper right", fontsize=20)
    plt.tight_layout()

    if scale == 'log':
        ax.set_yscale('log')

    outname = f"{path}/eff_{var}_{suf}_{scale}.{save_as}"
    fig.savefig(outname)
    plt.close(fig)
    print(f"Saved {outname}")




for var in var_configs:
    for ch in channels:
        plot_kinematic(var, ch)

efficiency_vars = ["Gen_tau_pt", "Gen_tau_eta", "Reco_jet_pT","Reco_jet_eta"]
for var in efficiency_vars:
    for ch in channels:
        plot_efficiency(var, ch, save_as='png', scale='lin')

Reco_efficiency_vars = ["Reco_jet_pT", "Reco_jet_eta"]
# for var in Reco_efficiency_vars:
for ch in channels:
    plot_reco_efficiency("Gen_bjet_pt", ch, save_as='png', scale='lin', both=True)
    plot_reco_efficiency("Gen_bjet_eta", ch, save_as='png', scale='lin', both=True)







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

