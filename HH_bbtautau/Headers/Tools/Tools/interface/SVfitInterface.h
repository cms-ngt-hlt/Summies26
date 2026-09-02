#ifndef SVfitInterface_h
#define SVfitInterface_h

// -------------------------------------------------------------------------------------------------------------- //
//                                                                                                                //
//   class SVfitInterface                                                                                         //
//                                                                                                                //
//   Class to compute ClassicSVfit. Modified from                                                                 //
//    - https://github.com/LLRCMS/KLUBAnalysis/blob/VBF_legacy/src/SVfitKLUBinterface.cc                          //
//                                                                                                                //
//   Based on:                                                                                                    //
//    - https://github.com/LLRCMS/LLRHiggsTauTau/blob/102X_HH/NtupleProducer/plugins/ClassicSVfitInterface.cc     //
//                                                                                                                //
//   Original SVfit package from:                                                                                 //
//    - https://github.com/SVfit/ClassicSVfit                                                                     //
//                                                                                                                //
//   Author: Francesco Brivio (Milano-Bicocca)                                                                    //
//   Date  : April 2020                                                                                           //
//                                                                                                                //
// -------------------------------------------------------------------------------------------------------------- //

// Standard libraries
#include <vector>
#include <string>
#include <cmath>
#include <chrono>

// ROOT libraries
#include <TLorentzVector.h>
#include "TMatrixD.h"

// ClassicSVfit libraries
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfitIntegrand.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

using namespace classic_svFit;

// SVfitInterface class
class SVfitInterface {

  public:
  SVfitInterface ();
  ~SVfitInterface ();
  std::vector<float> FitAndGetResultWithInputs(
    int verbosity, int pairType, int DM1, int DM2,
    double tau1_pt, double tau1_eta, double tau1_phi, double tau1_mass,
    double tau2_pt, double tau2_eta, double tau2_phi, double tau2_mass,
    double met_pt, double met_phi,
    double met_covXX, double met_covXY, double met_covYY);

};

#endif // SVfitInterface_h