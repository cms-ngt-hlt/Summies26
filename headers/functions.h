#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h> 

using ROOT::RVecF;
using ROOT::RVecI;
using ROOT::Math::PtEtaPhiMVector;

int Get_H_idx(int pdg, const RVecI pdgId, const RVecI status, const RVecI idx_mother);
RVecI Get_tau_from_H_indices(const RVecI& pdgId, const RVecI& status, const RVecI& idx_mother);
int tau_decay_mode(int tau_i, const RVecI& pdgId, const RVecI& status, const RVecI& idx_mother);
int Get_tau_channel(const RVecI& tau_indices, const RVecI& pdgId, const RVecI& status, const RVecI& idx_mother);
float Get_mHH(const RVecF pt, const RVecF eta, const RVecF phi, const RVecF M, int idx_tautau, int idx_bb);
RVecI Get_tau_pt(const RVecI& tau_indices, const RVecI& pt);


#endif