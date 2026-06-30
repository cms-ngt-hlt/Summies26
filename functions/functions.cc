#include <vector>
#include <cmath>
#include <algorithm>
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>         

using ROOT::RVecF;
using ROOT::RVecI;
using ROOT::Math::PtEtaPhiMVector; 


int Get_H_idx(int pdg, const RVecI& pdgId, const RVecI& status, const RVecI& idx_mother)
{
    for (std::size_t i = 0; i < pdgId.size(); i++){
        if (!(status[i] & (1 << 13))) continue;
        
        int my_pdg = pdgId.at(i);                               // PDG of particle p at idx i 
        int my_idx_mother = idx_mother.at(i);                   // idx of mother of p
        if (my_idx_mother < 0) continue;                        // This line prevents crash if particle has no mother and returns -1
        int my_pdg_mother = pdgId.at(my_idx_mother);            // PDG of particle at idx mother

        if(std::abs(my_pdg) == pdg){                            // If particle p is a tau, we look at its parents
            int current_mother_idx = my_idx_mother;             // Define idx of current parent

            while (current_mother_idx >= 0) {                   // Mother is a tau again? Define grandmother until not a tau anymore
                int mother_pdg = pdgId.at(current_mother_idx);  // Get PDG of current parent
                if (std::abs(mother_pdg) != pdg) break;         // If not 15, go to the next step, if 15, redefine grandmother
                current_mother_idx = idx_mother.at(current_mother_idx);
            }

            if (current_mother_idx >= 0 && pdgId.at(current_mother_idx) == 25){
                // std::cout << my_pdg <<  ", "  << pdgId.at(current_mother_idx) << std::endl;       
                return current_mother_idx;                       // Mother is a Higgs? Return index Higgs
            }
        }     
    }
    return -1;
}

RVecI Get_tau_from_H_indices(const RVecI& pdgId, const RVecI& status, const RVecI& idx_mother)
{
    RVecI result;
    int n = pdgId.size();

    for (int i = 0; i < n; i++) {
        if (!(status[i] & (1 << 13))) continue;
        if (std::abs(pdgId[i]) != 15) continue;
             
        int current = idx_mother[i];
        while (current >= 0 && std::abs(pdgId[current]) == 15)
            current = idx_mother[current];

        if (current >= 0 && pdgId[current] == 25)
            result.push_back(i);
    }

    if (result.size() != 2) return {-1, -1};
    return result;
}

RVecI Get_b_from_H_indices(const RVecI& pdgId, const RVecI& status, const RVecI& idx_mother)
{
    RVecI result;
    int n = pdgId.size();

    for (int i = 0; i < n; i++) {
        if (!(status[i] & (1 << 13))) continue;
        if (std::abs(pdgId[i]) != 5) continue;
             
        int current = idx_mother[i];
        while (current >= 0 && std::abs(pdgId[current]) == 5)
            current = idx_mother[current];

        if (current >= 0 && pdgId[current] == 25)
            result.push_back(i);
    }

    if (result.size() != 2) return {-1, -1};
    return result;
}

int tau_decay_mode(int tau_i, const RVecI& pdgId, const RVecI& status, const RVecI& idx_mother)
{
    int n = pdgId.size();
    for (int i = 0; i < n; i++) {
        if (idx_mother[i] != tau_i) continue;         // must be a direct daughter
        int pdg = std::abs(pdgId[i]);
        if (pdg == 11) return 0;                       // electron: leptonic e
        if (pdg == 13) return 1;                       // muon    : leptonic μ
    }
    return 2;                                          // No lepton? hadronic
}

int Get_tau_channel(const RVecI& tau_indices, const RVecI& pdgId, const RVecI& status, const RVecI& idx_mother)
{
    if (tau_indices[0] < 0 || tau_indices[1] < 0) return -1;

    int mode0 = tau_decay_mode(tau_indices[0], pdgId, status, idx_mother);
    int mode1 = tau_decay_mode(tau_indices[1], pdgId, status, idx_mother);

    if (mode0 > mode1) std::swap(mode0, mode1);

    if (mode0 == 0 && mode1 == 2) return 0;  // e + hadronic
    if (mode0 == 1 && mode1 == 2) return 1;  // mu + hadronic
    if (mode0 == 2 && mode1 == 2) return 2;  // hadronic + hadronic
    return -1;                               // anything that's not mentioned above (e+e, e+mu...)
}

float Get_mHH(const RVecF pt, const RVecF eta, const RVecF phi, const RVecF M, int idx_tautau, int idx_bb)
{
    PtEtaPhiMVector Higgs_tautau(pt.at(idx_tautau), eta.at(idx_tautau), phi.at(idx_tautau), M.at(idx_tautau));
    PtEtaPhiMVector Higgs_bb(pt.at(idx_bb), eta.at(idx_bb), phi.at(idx_bb), M.at(idx_bb));

    PtEtaPhiMVector HH = Higgs_tautau + Higgs_bb;

    return HH.M();
}

RVecF Get_tau_pt(const RVecI& tau_indices, const RVecF& pt)
{
    if (tau_indices[0] < 0 || tau_indices[1] < 0) return {-1.f, -1.f};

    float pt1 = pt[tau_indices[0]];
    float pt2 = pt[tau_indices[1]];

    if (pt1 >= pt2) return {pt1, pt2};
    else return {pt2, pt1};
}

RVecF Get_tau_eta(const RVecI& tau_indices, const RVecF& pt, const RVecF& eta)
{
    int tau1 = -1; int tau2 = -1;
    if(pt[tau_indices[0]] > pt[tau_indices[1]]){
        tau1 = tau_indices[0]; tau2 = tau_indices[1];
    }
    else tau1 = tau_indices[1]; tau2 = tau_indices[0];
    
    float eta1 = eta[tau1];
    float eta2 = eta[tau2];

    return {eta1, eta2};
}


//////////////////////////////////// RECO FUNCTIONS

ROOT::RVec<PtEtaPhiMVector> get_4vector(const RVecF& pt, const RVecF& eta, const RVecF& phi, const RVecF& M)
{
    ROOT::RVec<PtEtaPhiMVector> four_vectors;
    four_vectors.reserve(pt.size());
    for (std::size_t i = 0; i < pt.size(); i++){
        four_vectors.emplace_back(pt.at(i), eta.at(i), phi.at(i), M.at(i));
    }
    return four_vectors;
}

float get_ID_threshold(const float pt)
{
    double t1 = 0.649, t2 = 0.441, t3 = 0.05, x1 = 35, x2 = 100, x3 = 300; 
    if (pt <= x1) return t1; 
    if (pt >= x3) return t3; 
    if (pt < x2) return (t2 - t1) / (x2 - x1) * (pt - x1) + t1; 
    return (t3 - t2) / (x3 - x2) * (pt - x2) + t2;
}

RVecI deltaR_matching(const RVecI& gen_tau_indices, const RVecF& Gen_eta, const RVecF& Gen_phi, const ROOT::RVec<PtEtaPhiMVector>& reco_tau, const RVecF& TauVSjet)
{
    RVecI matched_indices;
    std::vector<bool> reco_used(reco_tau.size(), false);
    
    for (std::size_t i = 0; i < gen_tau_indices.size(); i++){   
        int idx_gen = gen_tau_indices[i];        
        
        double min_dR = 999.0;
        int idx_reco_tau = -1;                                         // Loops through every index of gen tau i
        
        for (std::size_t j = 0; j < reco_tau.size(); j++){             // Loops through every index of reco tau j 
            if (reco_used[j]) continue;
            float id_threshold = get_ID_threshold(reco_tau[j].pt());
            if (TauVSjet[j] < id_threshold) continue;

            double D_eta = Gen_eta[idx_gen] - reco_tau[j].eta();
            double D_phi = Gen_phi[idx_gen] - reco_tau[j].phi();
            while (D_phi >  M_PI) D_phi -= 2.0 * M_PI;                  // M_PI is pi constant, this line checks the periodicity of the angle
            while (D_phi < -M_PI) D_phi += 2.0 * M_PI;                  // Making sure we remain between 0 en 180 degrees
            
            double current_dR = std::sqrt(D_eta*D_eta+ D_phi*D_phi);
            if (current_dR < min_dR) {
                min_dR = current_dR;
                idx_reco_tau = j;
            }
        }

        if (idx_reco_tau != -1 && min_dR <= 0.3) {
            matched_indices.push_back(idx_reco_tau);
            reco_used[idx_reco_tau] = true;
        } else {
            matched_indices.push_back(-1); 
        }
    }
    
    return matched_indices;
}

RVecF Get_DpT(const RVecI& gen_tau_idx, const RVecI& reco_tau_indices, const RVecF& pt_gen, const RVecF& pt_reco)
{
    RVecF DpT;

    for (std::size_t i = 0; i < gen_tau_idx.size(); i++){
        int idx_gen = gen_tau_idx[i];
        int idx_reco = reco_tau_indices[i];
        if (idx_gen < 0) continue; if (idx_reco < 0) continue;

        float gen_pt = pt_gen[idx_gen];
        float reco_pt = pt_reco[idx_reco];
        
        DpT.push_back(gen_pt - reco_pt);
    }
    return DpT;
}

RVecF Get_Deta(const RVecI& gen_tau_idx, const RVecI& reco_tau_indices, const RVecF& eta_gen, const RVecF& eta_reco)
{
    
    RVecF Deta;

    for (std::size_t i = 0; i < gen_tau_idx.size(); i++){
        int idx_gen = gen_tau_idx[i];
        int idx_reco = reco_tau_indices[i];
        if (idx_gen < 0) continue;
        if (idx_reco < 0) continue;

        float gen_eta = eta_gen[idx_gen];
        float reco_eta = eta_reco[idx_reco];
        
        Deta.push_back(gen_eta - reco_eta);
    }
    return Deta;
}

RVecF Get_Dphi(const RVecI& gen_tau_idx, const RVecI& reco_tau_indices, const RVecF& phi_gen, const RVecF& phi_reco)
{
    
    RVecF Dphi;

    for (std::size_t i = 0; i < gen_tau_idx.size(); i++){
        int idx_gen = gen_tau_idx[i];
        int idx_reco = reco_tau_indices[i];
        if (idx_gen < 0) continue;
        if (idx_reco < 0) continue;

        float gen_phi = phi_gen[idx_gen];
        float reco_phi = phi_reco[idx_reco];
        
        Dphi.push_back(gen_phi - reco_phi);
    }
    return Dphi;
}

RVecF Get_recotau_pt(const RVecI& tau_idx, const RVecF& pt)
{
    if (tau_idx[0] < 0 || tau_idx[1] < 0) return {-1.f, -1.f};

    RVecF tau_pt;
    for (std::size_t i = 0; i < tau_idx.size(); i++) {
        tau_pt.push_back(pt[tau_idx.at(i)]);
    }
    return tau_pt;
}

RVecF Get_dR(const RVecI& gen_tau_indices, const RVecI& reco_tau_indices, const RVecF& Gen_eta, const RVecF& Gen_phi, const RVecF& reco_eta, const RVecF& reco_phi)
{
    RVecF dRs;

    for (std::size_t i = 0; i < gen_tau_indices.size(); i++) {
        int idx_gen = gen_tau_indices[i];
        int idx_reco = reco_tau_indices[i];
        if (idx_gen < 0) continue;
        if (idx_reco < 0) continue;
    
        double D_eta = Gen_eta[idx_gen] - reco_eta[idx_reco];
        double D_phi = Gen_phi[idx_gen] - reco_phi[idx_reco];
        
        while (D_phi >  M_PI) D_phi -= 2.0 * M_PI;
        while (D_phi < -M_PI) D_phi += 2.0 * M_PI;
        
        dRs.push_back(std::sqrt(D_eta*D_eta + D_phi*D_phi));
        }
    return dRs;
}


RVecI deltaR_matching_jets(const RVecI& gen_b_indices, const RVecF& Gen_eta, const RVecF& Gen_phi, const ROOT::RVec<PtEtaPhiMVector>& reco_jet, const RVecF& hltAK4PuppiJet_DeepFlavour_prob_b, const RVecF& hltAK4PuppiJet_DeepFlavour_prob_bb, const RVecF& hltAK4PuppiJet_DeepFlavour_prob_c, const RVecF& hltAK4PuppiJet_DeepFlavour_prob_g, const RVecF& hltAK4PuppiJet_DeepFlavour_prob_lepb, const RVecF& hltAK4PuppiJet_DeepFlavour_prob_uds)
{
    RVecI matched_indices;
    std::vector<bool> reco_used(reco_jet.size(), false);
    
    for (std::size_t i = 0; i < gen_b_indices.size(); i++){   
        int idx_gen = gen_b_indices[i];        
        
        double min_dR = 999.0;
        int idx_reco_jet = -1;                                         // Loops through every index of gen tau i
        
        for (std::size_t j = 0; j < reco_jet.size(); j++){             // Loops through every index of reco tau j 
            if (reco_used[j]) continue;
            float total_prob = hltAK4PuppiJet_DeepFlavour_prob_b[j] + hltAK4PuppiJet_DeepFlavour_prob_bb[j] + 
                                hltAK4PuppiJet_DeepFlavour_prob_c[j] + hltAK4PuppiJet_DeepFlavour_prob_g[j] + 
                                hltAK4PuppiJet_DeepFlavour_prob_lepb[j] + hltAK4PuppiJet_DeepFlavour_prob_uds[j];
            if (total_prob == 0) continue;
            float prob_b = hltAK4PuppiJet_DeepFlavour_prob_b[j] / total_prob; 
            if (prob_b < 0.92) continue;

            double D_eta = Gen_eta[idx_gen] - reco_jet[j].eta();
            double D_phi = Gen_phi[idx_gen] - reco_jet[j].phi();
            while (D_phi >  M_PI) D_phi -= 2.0 * M_PI;                  // M_PI is pi constant, this line checks the periodicity of the angle
            while (D_phi < -M_PI) D_phi += 2.0 * M_PI;                  // Making sure we remain between 0 en 180 degrees
            
            double current_dR = std::sqrt(D_eta*D_eta+ D_phi*D_phi);
            if (current_dR < min_dR) {
                min_dR = current_dR;
                idx_reco_jet = j;
            }
        }

        if (idx_reco_jet != -1 && min_dR <= 0.3) {
            matched_indices.push_back(idx_reco_jet);
            reco_used[idx_reco_jet] = true;
        } else {
            matched_indices.push_back(-1); 
        }
    }
    return matched_indices;
}
