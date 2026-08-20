#include <vector>
#include <cmath>
#include <algorithm>
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>         

using ROOT::RVecF;
using ROOT::RVecI;
using RVecUC = ROOT::VecOps::RVec<unsigned char>;

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
    return -999;
}

RVecI Get_b_from_H_indices(const RVecI& pdgId, const RVecI& status, const RVecI& idx_mother, bool first)
{
    RVecI result;
    for (std::size_t i = 0; i < pdgId.size(); i++) {
        if (first == true) { if (!(status[i] & (1 << 12))) continue;}
        else if (first == false) { if (!(status[i] & (1 << 13))) continue;}       
        
        if (std::abs(pdgId[i]) != 5) continue;

        int mom = idx_mother[i];
        while (mom >= 0 && std::abs(pdgId[mom]) == 5) {
            mom = idx_mother[mom];
        }

        if (mom >= 0 && pdgId[mom] == 25) result.push_back(i);
    }
    return result;
}

RVecI GetPDG(const RVecI& idx, const RVecI& pdgId)
{
    RVecI result;
    for (std::size_t i = 0; i < idx.size(); i++) {
        result.push_back(pdgId[idx[i]]);
    }
    return result;
}

RVecI GetPDG_mother(const RVecI& idx, const RVecI& pdgId, const RVecI& mother_idx)
{
    RVecI result;
    for (std::size_t i = 0; i < idx.size(); i++) {
        result.push_back(pdgId[mother_idx[idx[i]]]);
    }
    return result;
}

RVecI Match_b_to_GenJet(const RVecI& gen_b_idx, const RVecF& GenPart_eta, const RVecF& GenPart_phi, const RVecF& GenJet_eta, const RVecF& GenJet_phi)
{
    RVecI matched_indices;
    std::vector<bool> genjet_used(GenJet_eta.size(), false);

    for (std::size_t i = 0; i < gen_b_idx.size(); i++) {
        int idx_b = gen_b_idx[i];
        if (idx_b < 0) { matched_indices.push_back(-1); continue; }

        double min_dR = 999.0;
        int idx_genjet = -1;

        for (std::size_t j = 0; j < GenJet_eta.size(); j++) {
            if (genjet_used[j]) continue;

            double D_eta = GenPart_eta[idx_b] - GenJet_eta[j];
            double D_phi = GenPart_phi[idx_b] - GenJet_phi[j];
            while (D_phi >  M_PI) D_phi -= 2.0 * M_PI;
            while (D_phi < -M_PI) D_phi += 2.0 * M_PI;

            double current_dR = std::sqrt(D_eta*D_eta + D_phi*D_phi);
            if (current_dR < min_dR) { min_dR = current_dR; idx_genjet = j; }
        }

        if (idx_genjet != -1 && min_dR <= 0.4) {  
            matched_indices.push_back(idx_genjet);
            genjet_used[idx_genjet] = true;
        } else {
            matched_indices.push_back(-1);
        }
    }

    // std::sort(matched_indices.begin(), matched_indices.end(),
    //           [&GenJet_pt](int a, int b) {
    //               if (a < 0) return false;      
    //               if (b < 0) return true;
    //               return GenJet_pt[a] > GenJet_pt[b];
    //           });

    return matched_indices;
}

RVecI Get_Flavour_for_Jets(const RVecI& jet_indices, const ROOT::VecOps::RVec<UChar_t>& GenJet_hadronFlavour)
{
    RVecI result;
    for (auto idx : jet_indices) {
        if (idx < 0) { result.push_back(-1); continue; } // unmatched
        result.push_back(static_cast<int>(GenJet_hadronFlavour[idx]));
    }
    return result;
}

RVecI Get_tau_from_H_indices(const RVecI& pdgId, const RVecI& status, const RVecI& idx_mother)
{
    RVecI result; // Returns a vector of size 2 that contains the two indices of the tau's coming from the H decay

    for (std::size_t i = 0; i < pdgId.size(); i++) {
        if (!(status[i] & (1 << 13))) continue;
        if (std::abs(pdgId[i]) != 15) continue;
             
        int current = idx_mother[i];
        while (current >= 0 && std::abs(pdgId[current]) == 15) // Goes up the entire tree of tau to tau to tau until mother is a Higgs
            current = idx_mother[current];

        if (current >= 0 && pdgId[current] == 25)
            result.push_back(i);
    }

    if (result.size() != 2) return {-1, -1};
    return result;
}


ROOT::RVec<PtEtaPhiMVector> Get_b_p4(const RVecI& b_indices, const RVecF& pt, const RVecF& eta, const RVecF& phi, const RVecF& mass)
{
    ROOT::RVec<PtEtaPhiMVector> result;
    for (std::size_t i = 0; i < b_indices.size(); i++) {
        int idx = b_indices[i];
        if (idx < 0) {
            PtEtaPhiMVector vector(-999.f, -999.f, -999.f, -999.f);
            result.push_back(vector);
            continue;
        }
        result.push_back(PtEtaPhiMVector(pt.at(idx), eta.at(idx), phi.at(idx), mass.at(idx)));
    }
    // std::sort(result.begin(), result.end(),
    //           [](const PtEtaPhiMVector& a, const PtEtaPhiMVector& b) {
    //               return a.Pt() > b.Pt();
    //           });

    return result;
}

int tau_decay_mode(int tau_i, const RVecI& pdgId, const RVecI& status, const RVecI& idx_mother)
{
    for (std::size_t i = 0; i < pdgId.size(); i++) {
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

RVecI Get_taus(const RVecI& tau_indices, const RVecI& pdgId, const RVecI& status, const RVecI& idx_mother, const RVecF& pt)
{
    RVecI results;

    for (std::size_t i = 0; i < pdgId.size(); i++) {
        if (std::abs(pdgId[i]) == 11 || std::abs(pdgId[i]) == 13) {
            if (idx_mother[i] == tau_indices[0]) {
                results.push_back(i);
                results.push_back(tau_indices[1]);
                break;
            }
            else if (idx_mother[i] == tau_indices[1]) {
                results.push_back(i);
                results.push_back(tau_indices[0]);
                break;
            }
        }
    }

    if (results.size() > 0) return results;

    RVecI had = tau_indices;
    if (had.size() == 2 && had[0] >= 0 && had[1] >= 0 && pt[had[0]] < pt[had[1]])
        std::swap(had[0], had[1]);

    return had;
}

RVecI Get_taus_mH(const RVecI& tau_indices, const RVecI& pdgId, const RVecI& status, const RVecI& idx_mother, const RVecF& pt)
{
    RVecI results;

    for (std::size_t i = 0; i < pdgId.size(); i++) {
        if (std::abs(pdgId[i]) == 11 || std::abs(pdgId[i]) == 13) {
            if (idx_mother[i] == tau_indices[0]) {
                results.push_back(tau_indices[0]);
                results.push_back(tau_indices[1]);
                break;
            }
            else if (idx_mother[i] == tau_indices[1]) {
                results.push_back(tau_indices[0]);
                results.push_back(tau_indices[1]);
                break;
            }
        }
    }

    if (results.size() > 0) return results;

    RVecI had = tau_indices;
    if (had.size() == 2 && had[0] >= 0 && had[1] >= 0 && pt[had[0]] < pt[had[1]])
        std::swap(had[0], had[1]);

    return had;
}

ROOT::RVec<PtEtaPhiMVector> Get_all_taus_p4(const RVecI& tau_indices, const RVecF& pt, const RVecF& eta, const RVecF& phi, const RVecF& mass)
{
    ROOT::RVec<PtEtaPhiMVector> result;
    for (std::size_t i = 0; i < tau_indices.size(); i++) {
        int idx = tau_indices[i];
        if (idx < 0) {
            result.push_back(PtEtaPhiMVector(-999.f, -999.f, -999.f, -999.f));
            continue;
        }
        result.push_back(PtEtaPhiMVector(pt.at(idx), eta.at(idx), phi.at(idx), mass.at(idx)));
    }
    // std::sort(result.begin(), result.end(),
    //           [](const PtEtaPhiMVector& a, const PtEtaPhiMVector& b) {
    //               return a.Pt() > b.Pt();
    //           });

    return result;
}

void collect_visible_daughters(int idx, const RVecI& pdgId, const RVecI& idx_mother, std::vector<int>& result)
{
    bool has_daughter = false;
    for (std::size_t j = 0; j < pdgId.size(); j++) {
        if (idx_mother[j] != idx) continue; // idx is index tau. If particle j has a mother that is the tau, then we take it as daughter
        has_daughter = true;

        int abs_pdg = std::abs(pdgId[j]);
        if (abs_pdg == 12 || abs_pdg == 14 || abs_pdg == 16) continue; // if daughter is neutrino, skip it

        collect_visible_daughters(static_cast<int>(j), pdgId, idx_mother, result); // collect all daughter that are not neutrinos
    }

    if (!has_daughter) {  // if tau has no daughters, take back tau itself 
        int abs_pdg = std::abs(pdgId[idx]);
        if (abs_pdg != 12 && abs_pdg != 14 && abs_pdg != 16)
            result.push_back(idx);
    }
}


PtEtaPhiMVector Get_visible_tau_p4(int idx, const RVecI& pdgId, const RVecI& idx_mother, const RVecF& pt, const RVecF& eta, const RVecF& phi, const RVecF& mass)
{
    if (idx < 0)
        return PtEtaPhiMVector(-999.f, -999.f, -999.f, -999.f);

    std::vector<int> visible_daughters;
    collect_visible_daughters(idx, pdgId, idx_mother, visible_daughters);

    PtEtaPhiMVector sum(0.f, 0.f, 0.f, 0.f);
    for (int dIdx : visible_daughters) {
        PtEtaPhiMVector p4(pt[dIdx], eta[dIdx], phi[dIdx], mass[dIdx]);
        sum += p4;
    }
    return sum;
}  

ROOT::RVec<PtEtaPhiMVector> Get_visible_tau_p4s(const RVecI& tau_indices, const RVecI& pdgId, const RVecI& idx_mother, const RVecF& pt, const RVecF& eta, const RVecF& phi, const RVecF& mass)
{
    ROOT::RVec<PtEtaPhiMVector> result;
    result.reserve(tau_indices.size());
    for (std::size_t i = 0; i < tau_indices.size(); i++)
        result.push_back(Get_visible_tau_p4(tau_indices[i], pdgId, idx_mother, pt, eta, phi, mass));

    return result;
}


RVecF Get_p4_pt(const ROOT::RVec<PtEtaPhiMVector>& p4s)
{
    RVecF result;
    for (const auto& p4 : p4s) result.push_back(p4.pt());
    return result;
}
RVecF Get_p4_eta(const ROOT::RVec<PtEtaPhiMVector>& p4s)
{
    RVecF result;
    for (const auto& p4 : p4s) result.push_back(p4.eta());
    return result;
}
RVecF Get_p4_phi(const ROOT::RVec<PtEtaPhiMVector>& p4s)
{
    RVecF result;
    for (const auto& p4 : p4s) result.push_back(p4.phi());
    return result;
}
RVecF Get_p4_mass(const ROOT::RVec<PtEtaPhiMVector>& p4s)
{
    RVecF result;
    for (const auto& p4 : p4s) result.push_back(p4.M());
    return result;
}

RVecF Get_variable(const RVecI& indices, const RVecF& variable)
{
    RVecF result;
    for (std::size_t i = 0; i < indices.size(); i++) {
        if (indices[i] < 0) {
            result.push_back(-999.f);
            continue;
        }
        result.push_back(variable[indices[i]]);
    }
    return result;
}


RVecI Get_lepton_idx(int pdg, const RVecI& pdgId, const RVecI& tau_i, const RVecI& status, const RVecI& idx_mother)
{
    RVecI result; // Returns a vector of size n that contains the indices of the muons coming from tau decay (we expect only one)

    for (std::size_t i = 0; i < pdgId.size(); i++) {
        if (!(status[i] & (1 << 13))) continue;
        if (std::abs(pdgId[i]) != pdg) continue;     //Particle i must be a muon
        if (idx_mother[i] != tau_i[0] && idx_mother[i] != tau_i[1]) continue;    //Muon must be a direct daughter of tau
        result.push_back(i);
    }
    return result;
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


RVecI deltaR_matching(const RVecI& idx, const RVecI& pdg,
                       const RVecF& gen_pt, const RVecF& gen_eta, const RVecF& gen_phi,
                       const RVecF& pt_tau, const RVecF& eta_tau, const RVecF& phi_tau,
                       const RVecF& pt_electron, const RVecF& eta_electron, const RVecF& phi_electron,
                       const RVecF& pt_muon, const RVecF& eta_muon, const RVecF& phi_muon,
                       const RVecF& TauVSjet, bool use_tagging)
{
    RVecI matched_indices;

    std::vector<bool> tau_used(pt_tau.size(), false);
    std::vector<bool> electron_used(pt_electron.size(), false);
    std::vector<bool> muon_used(pt_muon.size(), false);

    for (std::size_t i = 0; i < gen_pt.size(); i++){
        if (gen_pt[i] < 0) { matched_indices.push_back(-1); continue; }

        int pdg_id = std::abs(pdg[idx[i]]);

        const RVecF* pt  = nullptr;
        const RVecF* eta = nullptr;
        const RVecF* phi = nullptr;
        std::vector<bool>* used = nullptr;
        bool is_tau = false;

        if (pdg_id == 15){ pt = &pt_tau;      eta = &eta_tau;      phi = &phi_tau;      used = &tau_used;      is_tau = true; }
        else if (pdg_id == 11){ pt = &pt_electron; eta = &eta_electron; phi = &phi_electron; used = &electron_used; }
        else if (pdg_id == 13){ pt = &pt_muon;     eta = &eta_muon;     phi = &phi_muon;     used = &muon_used; }
        else { matched_indices.push_back(-1); continue; }

        double min_dR = 999.0;
        int idx_reco = -1;

        for (std::size_t j = 0; j < pt->size(); j++){
            if ((*used)[j]) continue;

            if (is_tau && use_tagging){
                float id_threshold = get_ID_threshold((*pt)[j]);
                if (TauVSjet[j] < id_threshold) continue;
            }

            double D_eta = gen_eta[i] - (*eta)[j];
            double D_phi = gen_phi[i] - (*phi)[j];
            while (D_phi >  M_PI) D_phi -= 2.0 * M_PI;
            while (D_phi < -M_PI) D_phi += 2.0 * M_PI;

            double current_dR = std::sqrt(D_eta*D_eta + D_phi*D_phi);
            if (current_dR < min_dR) { min_dR = current_dR; idx_reco = j; }
        }

        if (idx_reco != -1 && min_dR <= 0.3) {
            matched_indices.push_back(idx_reco);
            (*used)[idx_reco] = true;
        } else {
            matched_indices.push_back(-1);
        }
    }
    return matched_indices;
}

RVecI deltaR_matching_lepton(const RVecF& gen_pt, const RVecF& gen_eta, const RVecF& gen_phi, const RVecF& reco_pt, const RVecF& reco_eta, const RVecF& reco_phi)
{
    RVecI matched_indices;
    std::vector<bool> reco_used(reco_pt.size(), false);

    for (std::size_t i = 0; i < gen_pt.size(); i++){
        if (gen_pt[i] < 0) {matched_indices.push_back(-1); continue;}

        double min_dR = 999.0;
        int idx_reco_lep = -1;

        for (std::size_t j = 0; j < reco_pt.size(); j++){
            if (reco_used[j]) continue;

            double D_eta = gen_eta[i] - reco_eta[j];
            double D_phi = gen_phi[i] - reco_phi[j];
            while (D_phi >  M_PI) D_phi -= 2.0 * M_PI;
            while (D_phi < -M_PI) D_phi += 2.0 * M_PI;

            double current_dR = std::sqrt(D_eta*D_eta + D_phi*D_phi);
            if (current_dR < min_dR) {min_dR = current_dR; idx_reco_lep = j;}
        }

        if (idx_reco_lep != -1 && min_dR <= 0.3) {
            matched_indices.push_back(idx_reco_lep);
            reco_used[idx_reco_lep] = true;
        } else {
            matched_indices.push_back(-1);
        }
    }
    return matched_indices;
}


RVecI deltaR_matching_jets(const RVecI& gen_b_idx, const RVecF& Gen_eta, const RVecF& Gen_phi, const RVecF& reco_eta, const RVecF& reco_phi,
                            const RVecF& prob_b, const RVecF& prob_bb, const RVecF& prob_c, const RVecF& prob_g, const RVecF& prob_lepb, const RVecF& prob_uds, bool use_tagging)
{
    RVecI matched_indices;
    std::vector<bool> reco_used(reco_eta.size(), false);

    for (std::size_t i = 0; i < gen_b_idx.size(); i++){
        int idx_gen = gen_b_idx[i];
        if (idx_gen < 0) { matched_indices.push_back(-1); continue; }

        double min_dR = 999.0;
        int idx_reco_jet = -1;

        for (std::size_t j = 0; j < reco_eta.size(); j++){
            if (reco_used[j]) continue;
            float total_prob = prob_b[j] + prob_bb[j] + prob_c[j] + prob_g[j] + prob_lepb[j] + prob_uds[j];
            if (total_prob == 0) continue;
            float b_disc = (prob_b[j] + prob_bb[j] + prob_lepb[j]) / total_prob;
            if (use_tagging){
                if (b_disc < 0.92) continue;
            }

            double D_eta = Gen_eta[idx_gen] - reco_eta[j];
            double D_phi = Gen_phi[idx_gen] - reco_phi[j];
            while (D_phi >  M_PI) D_phi -= 2.0 * M_PI;
            while (D_phi < -M_PI) D_phi += 2.0 * M_PI;

            double current_dR = std::sqrt(D_eta*D_eta + D_phi*D_phi);
            if (current_dR < min_dR) { min_dR = current_dR; idx_reco_jet = j; }
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

RVecF Get_bscore(const RVecI& jet_idx, const RVecF& prob_b, const RVecF& prob_bb, const RVecF& prob_c, const RVecF& prob_g, const RVecF& prob_lepb, const RVecF& prob_uds){
    
    RVecF score;

    for (std::size_t i=0; i < jet_idx.size() ; i++){
        int idx = jet_idx[i];
        float total_prob = prob_b[idx] + prob_bb[idx] + prob_c[idx] + prob_g[idx] + prob_lepb[idx] + prob_uds[idx];
        if (total_prob == 0) continue;
        float bscore = (prob_b[idx] + prob_bb[idx] + prob_lepb[idx]) / total_prob;
        score.push_back(bscore); 
    }
    return score;
}

RVecF Get_dR_tau(const ROOT::RVec<PtEtaPhiMVector>& gen_vis_tau_p4, const RVecI& matched_tau_idx,
                  const RVecF& reco_eta, const RVecF& reco_phi)
{
    RVecF dRs;
    for (std::size_t i = 0; i < gen_vis_tau_p4.size(); i++) {
        int idx_reco = matched_tau_idx[i];
        if (gen_vis_tau_p4[i].pt() < 0 || idx_reco < 0) {
            dRs.push_back(-999.f);
            continue;
        }

        double D_eta = gen_vis_tau_p4[i].eta() - reco_eta[idx_reco];
        double D_phi = gen_vis_tau_p4[i].phi() - reco_phi[idx_reco];
        while (D_phi >  M_PI) D_phi -= 2.0 * M_PI;
        while (D_phi < -M_PI) D_phi += 2.0 * M_PI;

        dRs.push_back(std::sqrt(D_eta*D_eta + D_phi*D_phi));
    }
    return dRs;
}

RVecF Get_dR_jet(const RVecI& gen_b_idx, const RVecI& matched_jet_idx,
                  const RVecF& gen_eta, const RVecF& gen_phi,
                  const RVecF& reco_eta, const RVecF& reco_phi)
{
    RVecF dRs;
    for (std::size_t i = 0; i < gen_b_idx.size(); i++) {
        int idx_gen  = gen_b_idx[i];
        int idx_reco = matched_jet_idx[i];

        if (idx_gen < 0 || idx_reco < 0) {
            dRs.push_back(-999.f);
            continue;
        }

        double D_eta = gen_eta[idx_gen] - reco_eta[idx_reco];
        double D_phi = gen_phi[idx_gen] - reco_phi[idx_reco];
        while (D_phi >  M_PI) D_phi -= 2.0 * M_PI;
        while (D_phi < -M_PI) D_phi += 2.0 * M_PI;

        dRs.push_back(std::sqrt(D_eta*D_eta + D_phi*D_phi));
    }
    return dRs;
}

RVecF Get_DpT(const RVecI& reco_tau_indices, const RVecF& gen_pt, const RVecF& reco_pt)
{
    RVecF result;
    for (std::size_t i = 0; i < reco_tau_indices.size(); i++) {
        int idx_reco = reco_tau_indices[i];
        if (gen_pt[i] < 0 || idx_reco < 0) {
            result.push_back(-999.f);
            continue;
        }
        result.push_back(reco_pt[idx_reco] - gen_pt[i]);
    }
    return result;
}

RVecF Get_Deta(const RVecI& reco_tau_indices, const RVecF& gen_eta, const RVecF& reco_eta)
{
    RVecF result;
    for (std::size_t i = 0; i < reco_tau_indices.size(); i++) {
        int idx_reco = reco_tau_indices[i];
        if (gen_eta[i] < 0 || idx_reco < 0) {
            result.push_back(-999.f);
            continue;
        }
        result.push_back(reco_eta[idx_reco] - gen_eta[i]);
    }
    return result;
}

RVecF Get_Dphi(const RVecI& reco_tau_indices, const RVecF& gen_phi, const RVecF& reco_phi)
{
    RVecF result;
    for (std::size_t i = 0; i < reco_tau_indices.size(); i++) {
        int idx_reco = reco_tau_indices[i];
        if (gen_phi[i] < 0 || idx_reco < 0) {
            result.push_back(-999.f);
            continue;
        }
        double dphi =  reco_phi[idx_reco] - gen_phi[i];
        while (dphi >  M_PI) dphi -= 2.0 * M_PI;
        while (dphi < -M_PI) dphi += 2.0 * M_PI;
        result.push_back(dphi);
    }
    return result;
}

PtEtaPhiMVector Build_Higgs_p4(const RVecI& idx_tau, const RVecF& pt, const RVecF& eta, const RVecF& phi, const RVecF& mass)
{
    if (idx_tau.size() != 2 || idx_tau[0] < 0 || idx_tau[1] < 0)
        return PtEtaPhiMVector(-999.f, 0.f, 0.f, 0.f);

    PtEtaPhiMVector tau1(pt.at(idx_tau[0]), eta.at(idx_tau[0]), phi.at(idx_tau[0]), mass.at(idx_tau[0]));
    PtEtaPhiMVector tau2(pt.at(idx_tau[1]), eta.at(idx_tau[1]), phi.at(idx_tau[1]), mass.at(idx_tau[1]));
    
    return tau1+tau2;
}

float Get_Higgs_mass(const PtEtaPhiMVector& H)
{
    if (H.pt() < 0) return -999.f;
    return H.M();
}

PtEtaPhiMVector Get_HH(const PtEtaPhiMVector& H1, const PtEtaPhiMVector& H2)
{

    if (H1.Pt() < 0 || H2.Pt() < 0) return PtEtaPhiMVector(-999.f, -999.f, -999.f, -999.f);
    PtEtaPhiMVector HH = H1 + H2;
    return HH;
}

// float Get_mHH(const PtEtaPhiMVector& H1, const PtEtaPhiMVector& H2)
// {
//     if (H1.pt() < 0 || H2.pt() < 0) return -999.f;
//     return (H1 + H2).M();
// }

PtEtaPhiMVector Get_H_fromdecay(const ROOT::RVec<PtEtaPhiMVector>& p4)
{
    if (p4.size() != 2) return PtEtaPhiMVector(-999.f, -999.f, -999.f, -999.f);
    if (p4[0].pt() < 0 || p4[1].pt() < 0) return PtEtaPhiMVector(-999.f, -999.f, -999.f, -999.f);

    PtEtaPhiMVector H = p4[0] + p4[1];
    return H;
}

PtEtaPhiMVector Get_H_fromdecay_idx(const RVecI& pair_idx, const RVecF& pt, const RVecF& eta, const RVecF& phi, const RVecF& mass)
{
    PtEtaPhiMVector H;
    for (std::size_t i = 0; i < pair_idx.size(); i++){
        int idx = pair_idx[i];
        if (idx < 0) { return (PtEtaPhiMVector(-999.f, -999.f, -999.f, -999.f));
            continue;
        }
        H += PtEtaPhiMVector(pt.at(idx), eta.at(idx), phi.at(idx), mass.at(idx));
    }
    return H;
}