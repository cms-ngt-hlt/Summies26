import ROOT
ROOT.gROOT.SetBatch(True)
import os
import shutil, os

output = "L1_efficiencies.txt"
SRC = "/eos/user/s/sbenabde/CERN_Summer_student/functions/functions.cc"
HDR = "/eos/user/s/sbenabde/CERN_Summer_student/headers/functions.h"
TMP = "/tmp/functions.cc"

shutil.copy(SRC, TMP)
ROOT.gROOT.ProcessLine(f'.L {TMP}+')

ROOT.ROOT.EnableImplicitMT()
df = ROOT.RDataFrame("Events", "Dataframes/df_ngt_hlt.root")

def overall_efficiency(all, trigger):
        all_c = all.Count().GetValue()
        trigger_c = trigger.Count().GetValue()
        eff = trigger_c/all_c
        return eff 

df = (
    df
    .Define("tau_indices",           "Get_tau_from_H_indices(GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")
    .Define("tau_channel",           "Get_tau_channel(tau_indices, GenPart_pdgId, GenPart_statusFlags, GenPart_genPartIdxMother)")

    .Define("tau_pt",                "Get_tau_pt(tau_indices, GenPart_pt)")
    .Define("tau_pt_lead",           "tau_pt[0]")
    .Define("tau_pt_sublead",        "tau_pt[1]")

    .Define("tau_eta",               "Get_tau_eta(tau_indices, GenPart_pt, GenPart_eta)")
    .Define("eta1",                  "tau_eta[0]")
    .Define("eta2",                  "tau_eta[1]")
)
df_L1 = df.Filter("tau_channel == 2")  

data = {}

all_cols = [str(c) for c in df.GetColumnNames()]
L1_cols = [c for c in all_cols if "L1_p" in c]

print("\nComputing the efficiencies...")
for c in sorted(L1_cols):
    column_name = f"{c}  [{df.GetColumnType(c)}]"
    df_L1_filtered = df_L1.Filter(f"{c} == true")
    efficiency_L1 = str(overall_efficiency(df_L1, df_L1_filtered))
    data[column_name] = efficiency_L1
    
sorted_data = {k: v for k, v in sorted(data.items(), key=lambda item: item[1])}

with open(output, "w") as out:
    for key, value in sorted_data.items():
        v = float(value) * 100
        out.write(f"{key:<70} {v:.2f} %\n")

print(f"Output saved in {output}")