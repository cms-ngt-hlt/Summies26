#include <vector>
#include <cmath>
#include <algorithm>
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>         

using ROOT::RVecF;
using ROOT::RVecI;
using ROOT::Math::PtEtaPhiMVector; 

RVecI get_muon_idx(const RVecI& pdg, const RVecI& parent, const RVecI& status)
{
    RVecI results;
    for (std::size_t i = 0; i < pdg.size(); ++i)
    {
        if (!(status[i] & (1 << 13))) continue;
        
        int pdg_i = pdg.at(i);
        int parent_idx = parent.at(i);
        if (std::abs(pdg_i) == 13 and std::abs(pdg.at(parent_idx)) == 1023){
            results.push_back(i);
        }
    }
    return results;
}

RVecF get_variable(const RVecI& idx, const RVecF& variable)
{
    RVecF results;
    for (std::size_t i = 0; i < idx.size(); i++){
        int idx_i = idx.at(i);
        if (idx_i < 0) continue;
        results.push_back(variable.at(idx_i));
    }
    return results;
}

RVecF get_pair_variable(const RVecI& idx, const RVecF& variable)
{
    RVecF results;
    if (idx.size() != 2) return {-999, -999};
    for (std::size_t i = 0; i < idx.size(); i++){
        int idx_i = idx.at(i);
        if (idx_i < 0) continue;
        results.push_back(variable.at(idx_i));
    }
    return results;
}

RVecI deltaR_muon_pairing(const RVecI& idx, const RVecF& eta, const RVecF& phi, const RVecI& idx_mother)
{
    RVecI results;

    double min_dR = 999.0;
    int idx1 = -1; int idx2 = -1;
   
    for (std::size_t i = 0; i < idx.size(); i++){
        if (idx.size() < 2) continue;
        for (std::size_t j = i+1; j < idx.size(); j++){
            int idx_i = idx.at(i); int idx_j = idx.at(j);
            if (idx_i < 0 || idx_j < 0 ) continue;

            double D_eta = eta.at(idx_i) - eta.at(idx_j);
            double D_phi = phi.at(idx_i) - phi.at(idx_j);
            while (D_phi >  M_PI) D_phi -= 2.0 * M_PI; 
            while (D_phi < -M_PI) D_phi += 2.0 * M_PI;

            double current_dR = std::sqrt(D_eta*D_eta + D_phi*D_phi);
            
            if(current_dR < min_dR && idx_mother.at(idx_i)==idx_mother.at(idx_j)){
                min_dR = current_dR;
                idx1 = idx_i; idx2 = idx_j;

            }
        }
    }  
    if (idx1 != -1){
        results.push_back(idx1);
        results.push_back(idx2);
        return results;
    }  

    return {};
}

ROOT::RVec<PtEtaPhiMVector> get_pair_4V(const RVecI& pair_idx, const RVecF& pt, const RVecF& eta, const RVecF& phi, const RVecF& mass)
{
    ROOT::RVec<PtEtaPhiMVector> results;

    for (std::size_t i = 0; i < pair_idx.size(); i++){
        int idx = pair_idx[i];
        if (idx < 0) {
            results.push_back(PtEtaPhiMVector(-999.f, -999.f, -999.f, -999.f));
            continue;
        }
        results.push_back(PtEtaPhiMVector(pt.at(idx), eta.at(idx), phi.at(idx), mass.at(idx)));
    }
    return results;
}

float get_pair_pt(const RVecI& pair_idx, const RVecF& pt, const RVecF& eta, const RVecF& phi, const RVecF& mass)
{
    if (pair_idx.size() != 2) return -999.f;

    PtEtaPhiMVector mu1(pt.at(pair_idx[0]), eta.at(pair_idx[0]), phi.at(pair_idx[0]), mass.at(pair_idx[0]));
    PtEtaPhiMVector mu2(pt.at(pair_idx[1]), eta.at(pair_idx[1]), phi.at(pair_idx[1]), mass.at(pair_idx[1]));

    return (mu1 + mu2).Pt();
}


float get_pair_eta(const RVecI& pair_idx, const RVecF& pt, const RVecF& eta, const RVecF& phi, const RVecF& mass)
{
    if (pair_idx.size() != 2) return -999.f;

    PtEtaPhiMVector mu1(pt.at(pair_idx[0]), eta.at(pair_idx[0]), phi.at(pair_idx[0]), mass.at(pair_idx[0]));
    PtEtaPhiMVector mu2(pt.at(pair_idx[1]), eta.at(pair_idx[1]), phi.at(pair_idx[1]), mass.at(pair_idx[1]));

    return (mu1 + mu2).Eta();
}

float get_pair_phi(const RVecI& pair_idx, const RVecF& pt, const RVecF& eta, const RVecF& phi, const RVecF& mass)
{
    if (pair_idx.size() != 2) return -999.f;

    PtEtaPhiMVector mu1(pt.at(pair_idx[0]), eta.at(pair_idx[0]), phi.at(pair_idx[0]), mass.at(pair_idx[0]));
    PtEtaPhiMVector mu2(pt.at(pair_idx[1]), eta.at(pair_idx[1]), phi.at(pair_idx[1]), mass.at(pair_idx[1]));

    return (mu1 + mu2).Phi();
}

float get_pair_M(const RVecI& pair_idx, const RVecF& pt, const RVecF& eta, const RVecF& phi, const RVecF& mass)
{
    if (pair_idx.size() != 2) return -999.f;

    PtEtaPhiMVector mu1(pt.at(pair_idx[0]), eta.at(pair_idx[0]), phi.at(pair_idx[0]), mass.at(pair_idx[0]));
    PtEtaPhiMVector mu2(pt.at(pair_idx[1]), eta.at(pair_idx[1]), phi.at(pair_idx[1]), mass.at(pair_idx[1]));

    return (mu1 + mu2).M();
}

RVecI truth_matching(const RVecI& indices, const RVecF& gen_eta, const RVecF& gen_phi, const RVecF& reco_pt, const RVecF& reco_eta, const RVecF& reco_phi)
{
    RVecI matched_indices;
    std::vector<bool> reco_used(reco_pt.size(), false);

    for (std::size_t i = 0; i < indices.size(); i++){
        if (indices[i] < 0) { matched_indices.push_back(-1); continue; }
        
        int gen_idx = static_cast<int>(indices[i]);

        double min_dR = 999.0;
        int idx_reco_mu = -1;

        for (std::size_t j = 0; j < reco_pt.size(); j++){
            if (reco_used[j]) continue;

            double D_eta = gen_eta[gen_idx] - reco_eta[j];
            double D_phi = gen_phi[gen_idx] - reco_phi[j];
            while (D_phi >  M_PI) D_phi -= 2.0 * M_PI;
            while (D_phi < -M_PI) D_phi += 2.0 * M_PI;

            double current_dR = std::sqrt(D_eta*D_eta + D_phi*D_phi);
            if (current_dR < min_dR) { min_dR = current_dR; idx_reco_mu = j; }
        }

        if (idx_reco_mu != -1 && min_dR <= 0.3) {
            matched_indices.push_back(idx_reco_mu);
            reco_used[idx_reco_mu] = true;
        } else {
            matched_indices.push_back(-999);
        }
    }
    return matched_indices;
}

RVecF get_reco_variable(const RVecI& idx, const RVecF& variable)
{
    RVecF results;
    if (idx.size() != 2) return {-999, -999};
    for (std::size_t i = 0; i < idx.size(); ++i){
        if (idx[i] < 0) {
            results.push_back(-999);
            continue;
        }
        results.push_back(variable.at(idx.at(i)));
    }
    return results;
}

float get_pair_pt(const RVecI& pair_idx, const RVecF& pt, const RVecF& eta, const RVecF& phi, float mass)
{
    if (pair_idx.size() != 2) return -999.f;
    if (pair_idx[0] < 0 || pair_idx[1] < 0) return -999.f;

    PtEtaPhiMVector mu1(pt.at(pair_idx[0]), eta.at(pair_idx[0]), phi.at(pair_idx[0]), mass);
    PtEtaPhiMVector mu2(pt.at(pair_idx[1]), eta.at(pair_idx[1]), phi.at(pair_idx[1]), mass);

    return (mu1 + mu2).Pt();
}

float get_pair_eta(const RVecI& pair_idx, const RVecF& pt, const RVecF& eta, const RVecF& phi, float mass)
{
    if (pair_idx.size() != 2) return -999.f;
    if (pair_idx[0] < 0 || pair_idx[1] < 0) return -999.f;

    PtEtaPhiMVector mu1(pt.at(pair_idx[0]), eta.at(pair_idx[0]), phi.at(pair_idx[0]), mass);
    PtEtaPhiMVector mu2(pt.at(pair_idx[1]), eta.at(pair_idx[1]), phi.at(pair_idx[1]), mass);

    return (mu1 + mu2).Eta();
}

float get_pair_phi(const RVecI& pair_idx, const RVecF& pt, const RVecF& eta, const RVecF& phi, float mass)
{
    if (pair_idx.size() != 2) return -999.f;
    if (pair_idx[0] < 0 || pair_idx[1] < 0) return -999.f;

    PtEtaPhiMVector mu1(pt.at(pair_idx[0]), eta.at(pair_idx[0]), phi.at(pair_idx[0]), mass);
    PtEtaPhiMVector mu2(pt.at(pair_idx[1]), eta.at(pair_idx[1]), phi.at(pair_idx[1]), mass);

    return (mu1 + mu2).Phi();
}

float get_pair_M(const RVecI& pair_idx, const RVecF& pt, const RVecF& eta, const RVecF& phi, float mass)
{
    if (pair_idx.size() != 2) return -999.f;
    if (pair_idx[0] < 0 || pair_idx[1] < 0) return -999.f;

    PtEtaPhiMVector mu1(pt.at(pair_idx[0]), eta.at(pair_idx[0]), phi.at(pair_idx[0]), mass);
    PtEtaPhiMVector mu2(pt.at(pair_idx[1]), eta.at(pair_idx[1]), phi.at(pair_idx[1]), mass);

    return (mu1 + mu2).M();
}

float get_dR(const RVecI& pair_idx, const RVecF& eta, const RVecF& phi)
{
    if (pair_idx.size() != 2) return -999;
    if (pair_idx[0] < 0 || pair_idx[1] < 0) return -999;

    double D_eta = eta[pair_idx[0]] - eta[pair_idx[1]];
    double D_phi = phi[pair_idx[0]] - phi[pair_idx[1]];
    while (D_phi >  M_PI) D_phi -= 2.0 * M_PI;
    while (D_phi < -M_PI) D_phi += 2.0 * M_PI;

    double dR = std::sqrt(D_eta*D_eta + D_phi*D_phi);
    return dR;
}