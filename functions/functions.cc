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
    for (int i = 0; i < pdgId.size(); i++){
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

RVecF Get_tau_pt(const RVecI& tau_indices, const RVecI& pt)
{
    if (tau_indices[0] < 0 || tau_indices[1] < 0) return {-1.f, -1.f};

    float pt1 = pt[tau_indices[0]];
    float pt2 = pt[tau_indices[1]];

    if (pt1 >= pt2) return {pt1, pt2};
    else return {pt2, pt1};
}