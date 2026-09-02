#include "Tools/Tools/interface/SVfitInterface.h"
SVfitInterface::SVfitInterface () {};
SVfitInterface::~SVfitInterface() {}
std::vector<float> SVfitInterface::FitAndGetResultWithInputs(
    int verbosity, int pairType, int DM1, int DM2,
    double tau1_pt, double tau1_eta, double tau1_phi, double tau1_mass,
    double tau2_pt, double tau2_eta, double tau2_phi, double tau2_mass,
    double met_pt, double met_phi,
    double met_covXX, double met_covXY, double met_covYY){

  // Start timing here
  auto start_time = std::chrono::high_resolution_clock::now();

  // Measured quantities
  std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons_;
  measuredTauLeptons_.clear();
  double METx_;
  double METy_;
  TMatrixD covMET_(2, 2);
  double kappa_;

  // MET
  covMET_(0,0) = met_covXX;
  covMET_(0,1) = met_covXY;
  covMET_(1,0) = met_covXY;
  covMET_(1,1) = met_covYY;
  METx_ = met_pt * cos(met_phi);
  METy_ = met_pt * sin(met_phi);

  // Leptons
  double mass1 = -1, mass2 = -1;
  int decay1 =  -1, decay2 = -1;
  classic_svFit::MeasuredTauLepton::kDecayType l1Type, l2Type;

  if (pairType == 0) // EleTau
   {
    l1Type = classic_svFit::MeasuredTauLepton::kTauToElecDecay;
    mass1  = 0.51100e-3;
    decay1 = -1;
    l2Type = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
    mass2  = tau2_mass;
    decay2 = DM2;
    kappa_ = 4.;
  }
  else if (pairType == 1) // MuTau
  {
    l1Type = classic_svFit::MeasuredTauLepton::kTauToMuDecay;
    mass1  = 105.658e-3;
    decay1 = -1;
    l2Type = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
    mass2  = tau2_mass;
    decay2 = DM2;
    kappa_ = 4.;
  }
  else if (pairType == 2)// TauTau
  {
    l1Type = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
    mass1  = tau1_mass;
    decay1 = DM1;
    l2Type = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
    mass2  = tau2_mass;
    decay2 = DM2;
    kappa_ = 5.;
  }

  // Fill the measuredTauLeptons
  if (pairType == 2 ) {
    if (tau1_pt > tau2_pt) {
      measuredTauLeptons_.push_back(classic_svFit::MeasuredTauLepton(l1Type, tau1_pt, tau1_eta, tau1_phi, mass1, decay1, verbosity));
      measuredTauLeptons_.push_back(classic_svFit::MeasuredTauLepton(l2Type, tau2_pt, tau2_eta, tau2_phi, mass2, decay2, verbosity));
    } else {
      measuredTauLeptons_.push_back(classic_svFit::MeasuredTauLepton(l2Type, tau2_pt, tau2_eta, tau2_phi, mass2, decay2, verbosity));
      measuredTauLeptons_.push_back(classic_svFit::MeasuredTauLepton(l1Type, tau1_pt, tau1_eta, tau1_phi, mass1, decay1, verbosity));
    }
  } else { // pairType 0 (etau) or 1 (mutau): lepton first, hadronic tau second
    measuredTauLeptons_.push_back(classic_svFit::MeasuredTauLepton(l1Type, tau1_pt, tau1_eta, tau1_phi, mass1, decay1, verbosity));
    measuredTauLeptons_.push_back(classic_svFit::MeasuredTauLepton(l2Type, tau2_pt, tau2_eta, tau2_phi, mass2, decay2, verbosity));
 }

  // Declare result: vector of SVfit (Pt,Eta,Phi,Mass,timing)
  std::vector<float> result(5,-999.);

  // Declare algo
  ClassicSVfit algo(verbosity);

  // Configure more options
  algo.addLogM_fixed(false, kappa_);
  algo.addLogM_dynamic(false);

  // Actually integrate
  algo.integrate(measuredTauLeptons_, METx_, METy_, covMET_);

  // Return SVfit quantities if the integration succeeded
  // otherwise vector of -999 if the integration failed
  if (algo.isValidSolution())
  {
    result.at(0) = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getPt();
    result.at(1) = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getEta();
    result.at(2) = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getPhi();
    result.at(3) = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getMass();
  }

  // Add timing information
  auto end_time = std::chrono::high_resolution_clock::now();
  result.at(4) = (end_time - start_time).count() * 1.e-9;

  return result;
}
