// The contents of this file are in the public domain. See LICENSE_FOR_EXAMPLE_PROGRAMS.txt
/*

    This is an example illustrating the use the general purpose non-linear 
    least squares optimization routines from the dlib C++ Library.

    This example program will demonstrate how these routines can be used for data fitting.
    In particular, we will generate a set of data and then use the least squares  
    routines to infer the parameters of the model which generated the data.
*/

#include <dlib/optimization.h>
#include <dlib/config_reader.h>
#include <dlib/error.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <math.h>
#include "datatypes.h"

using namespace std;
using namespace dlib;

// ----------------------------------------------------------------------------------------
double LNGIVIM_residual (const std::pair<double, double>& data, const vectorParam2& params);
double LNGIVIM_residual_weights (const std::pair<double, double>& data, const vectorParam2& params);
double LNGIVIM_residual_MSE (const matrix<std::pair<double, double>, 0, 1>& data, const matrix<double,0,1>&);
double NGIVIM_residual (const std::pair<double, double>& data, const vectorParam5& params);
double NGIVIM_residual_weights (const std::pair<double, double>& data, const vectorParam5& params);
double NGIVIM_residual_MSE (const matrix<std::pair<double, double>, 0, 1>& data, const matrix<double,0,1>&);
double MonoExponential_residual (const std::pair<double, double>& data, const vectorParam2& params);
double MonoExponential_residual_weights (const std::pair<double, double>& data, const vectorParam2& params);
double MonoExponential_residual_MSE (const matrix<std::pair<double, double>, 0, 1>& data, const matrix<double,0,1>&);
double MonoExponentialN_residual_MSE (const matrix<std::pair<double, double>, 0, 1>& data, const matrix<double,0,1>&);
double MonoExponentialN_residual_MSE_mtx (const matrix<std::pair<double, double>, 0, 1>& data, const matrix<double,0,1>&);
double MonoN_residual_MSE (const std::pair<double, double>& data, const matrix<double,1,1>& params);
double MonoT1rho_residual (const std::pair<double, double>& data, const vectorParam2& params);
double MonoT1rho_residual_weights (const std::pair<double, double>& data, const vectorParam2& params);
double MonoT1rho_residual_MSE (const matrix<std::pair<double, double>, 0, 1>& data, const matrix<double,0,1>&);
double Kurt_residual (const std::pair<double, double>& data, const vectorParam3& params);
double Kurt_residual_weights (const std::pair<double, double>& data, const vectorParam3& params);
double Kurt_const(const vectorParam3& params);
double Kurt_residual_MSE (const matrix<std::pair<double, double>, 0, 1>& data, const matrix<double,0,1>&);
double KurtN_residual_MSE (const matrix<std::pair<double, double>, 0, 1>& data, const matrix<double,0,1>&);
double Stretched_residual (const std::pair<double, double>& data, const vectorParam3& params);
double Stretched_residual_weights (const std::pair<double, double>& data, const vectorParam3& params);
double Stretched_residual_MSE (const matrix<std::pair<double, double>, 0, 1>& data, const matrix<double,0,1>&);
double StretchedN_residual_MSE (const matrix<std::pair<double, double>, 0, 1>& data, const matrix<double,0,1>&);
double Biexp_model(const double& xvalue, const vectorParam4& params);
double Biexp_residual (const std::pair<double, double>& data, const vectorParam4& params);
double Biexp_residual_weights (const std::pair<double, double>& data, const vectorParam4& params);
double Biexp_residual_fixedDs (const std::pair<double, double>& data, const vectorParam3& params);
double Biexp_residual_MSE (const matrix<std::pair<double, double>, 0, 1>& data, const matrix<double,0,1>&);
double Biexp_residual_MSE_fixedDs (const matrix<std::pair<double, double>, 0, 1>& data, const matrix<double,0,1>& params);
double Biexp_residual_MSE_fixedDsf (const matrix<std::pair<double, double>, 0, 1>& data, const matrix<double,0,1>& params);
double BiexpN_residual_fixedDs_mtx (const std::pair<double, double>& data, const vectorParam2& params);
double BiexpN_residual_MSE (const matrix<std::pair<double, double>, 0, 1>& data, const matrix<double,0,1>&);
double BiexpN_residual_MSE_fixedDs_mtx (const matrix<std::pair<double, double>, 0, 1>&, const matrix<double,0,1>&);
double BiexpN_residual_MSE_fixedDs (const std::pair<double, double>&, const matrix<double,2,1>&);

// ----------------------------------------------------------------------------------------
SIMU_IMPLEMENTATION(Simulate_4params, 4)
// ----------------------------------------------------------------------------------------
//define constrained non-normalized
FIT_CONST_MI_IMPLEMENTATION(FitConst_2params_2_bvalues_MultipleInitializations, 2, 2)
FIT_CONST_MI_IMPLEMENTATION(FitConst_2params_3_bvalues_MultipleInitializations, 2, 3)
FIT_CONST_MI_IMPLEMENTATION(FitConst_2params_4_bvalues_MultipleInitializations, 2, 4)
FIT_CONST_MI_IMPLEMENTATION(FitConst_2params_5_bvalues_MultipleInitializations, 2, 5)
FIT_CONST_MI_IMPLEMENTATION(FitConst_2params_6_bvalues_MultipleInitializations, 2, 6)
FIT_CONST_MI_IMPLEMENTATION(FitConst_2params_7_bvalues_MultipleInitializations, 2, 7)
FIT_CONST_MI_IMPLEMENTATION(FitConst_2params_8_bvalues_MultipleInitializations, 2, 8)
FIT_CONST_MI_IMPLEMENTATION(FitConst_2params_9_bvalues_MultipleInitializations, 2, 9)
FIT_CONST_MI_IMPLEMENTATION(FitConst_2params_10_bvalues_MultipleInitializations, 2, 10)
FIT_CONST_MI_IMPLEMENTATION(FitConst_2params_11_bvalues_MultipleInitializations, 2, 11)
FIT_CONST_MI_IMPLEMENTATION(FitConst_2params_12_bvalues_MultipleInitializations, 2, 12)
FIT_CONST_MI_IMPLEMENTATION(FitConst_2params_13_bvalues_MultipleInitializations, 2, 13)
FIT_CONST_MI_IMPLEMENTATION(FitConst_2params_14_bvalues_MultipleInitializations, 2, 14)
FIT_CONST_MI_IMPLEMENTATION(FitConst_2params_15_bvalues_MultipleInitializations, 2, 15)

FIT_CONST_MI_IMPLEMENTATION(FitConst_3params_2_bvalues_MultipleInitializations, 3, 2)
FIT_CONST_MI_IMPLEMENTATION(FitConst_3params_3_bvalues_MultipleInitializations, 3, 3)
FIT_CONST_MI_IMPLEMENTATION(FitConst_3params_4_bvalues_MultipleInitializations, 3, 4)
FIT_CONST_MI_IMPLEMENTATION(FitConst_3params_5_bvalues_MultipleInitializations, 3, 5)
FIT_CONST_MI_IMPLEMENTATION(FitConst_3params_6_bvalues_MultipleInitializations, 3, 6)
FIT_CONST_MI_IMPLEMENTATION(FitConst_3params_7_bvalues_MultipleInitializations, 3, 7)
FIT_CONST_MI_IMPLEMENTATION(FitConst_3params_8_bvalues_MultipleInitializations, 3, 8)
FIT_CONST_MI_IMPLEMENTATION(FitConst_3params_9_bvalues_MultipleInitializations, 3, 9)
FIT_CONST_MI_IMPLEMENTATION(FitConst_3params_10_bvalues_MultipleInitializations, 3, 10)
FIT_CONST_MI_IMPLEMENTATION(FitConst_3params_11_bvalues_MultipleInitializations, 3, 11)
FIT_CONST_MI_IMPLEMENTATION(FitConst_3params_12_bvalues_MultipleInitializations, 3, 12)
FIT_CONST_MI_IMPLEMENTATION(FitConst_3params_14_bvalues_MultipleInitializations, 3, 14)
FIT_CONST_MI_IMPLEMENTATION(FitConst_3params_15_bvalues_MultipleInitializations, 3, 15)
FIT_CONST_MI_IMPLEMENTATION(FitConst_5params_12_bvalues_MultipleInitializations, 5, 12)
FIT_CONST_MI_IMPLEMENTATION(FitConst_5params_14_bvalues_MultipleInitializations, 5, 14)
FIT_CONST_MI_IMPLEMENTATION(FitConst_5params_15_bvalues_MultipleInitializations, 5, 15)
FIT_CONST_MI_IMPLEMENTATION_BIEXP(FitBiexpConst_4params_2_bvalues_MultipleInitializations, 4, 2)
FIT_CONST_MI_IMPLEMENTATION_BIEXP(FitBiexpConst_4params_3_bvalues_MultipleInitializations, 4, 3)
FIT_CONST_MI_IMPLEMENTATION_BIEXP(FitBiexpConst_4params_4_bvalues_MultipleInitializations, 4, 4)
FIT_CONST_MI_IMPLEMENTATION_BIEXP(FitBiexpConst_4params_5_bvalues_MultipleInitializations, 4, 5)
FIT_CONST_MI_IMPLEMENTATION_BIEXP(FitBiexpConst_4params_6_bvalues_MultipleInitializations, 4, 6)
FIT_CONST_MI_IMPLEMENTATION_BIEXP(FitBiexpConst_4params_7_bvalues_MultipleInitializations, 4, 7)
FIT_CONST_MI_IMPLEMENTATION_BIEXP(FitBiexpConst_4params_8_bvalues_MultipleInitializations, 4, 8)
FIT_CONST_MI_IMPLEMENTATION_BIEXP(FitBiexpConst_4params_9_bvalues_MultipleInitializations, 4, 9)
FIT_CONST_MI_IMPLEMENTATION_BIEXP(FitBiexpConst_4params_10_bvalues_MultipleInitializations, 4, 10)
FIT_CONST_MI_IMPLEMENTATION_BIEXP(FitBiexpConst_4params_11_bvalues_MultipleInitializations, 4, 11)
FIT_CONST_MI_IMPLEMENTATION_BIEXP(FitBiexpConst_4params_12_bvalues_MultipleInitializations, 4, 12)
FIT_CONST_MI_IMPLEMENTATION_BIEXP(FitBiexpConst_4params_13_bvalues_MultipleInitializations, 4, 13)
FIT_CONST_MI_IMPLEMENTATION_BIEXP(FitBiexpConst_4params_14_bvalues_MultipleInitializations, 4, 14)
FIT_CONST_MI_IMPLEMENTATION_BIEXP(FitBiexpConst_4params_15_bvalues_MultipleInitializations, 4, 15)
FIT_CONST_MI_IMPLEMENTATION_MARZISEGMENTEDBIEXP(FitMarziSegmentedBiexpConst_params_12_bvalues_3tail_MultipleInitializations, 12, 3)
FIT_CONST_MI_IMPLEMENTATION_MARZISEGMENTEDBIEXP(FitMarziSegmentedBiexpConst_params_14_bvalues_3tail_MultipleInitializations, 14, 3)
FIT_CONST_MI_IMPLEMENTATION_MARZISEGMENTEDBIEXP(FitMarziSegmentedBiexpConst_params_15_bvalues_3tail_MultipleInitializations, 15, 3)
FIT_CONST_MI_IMPLEMENTATION_MARZISEGMENTEDBIEXP(FitMarziSegmentedBiexpConst_params_12_bvalues_4tail_MultipleInitializations, 12, 4)
FIT_CONST_MI_IMPLEMENTATION_MARZISEGMENTEDBIEXP(FitMarziSegmentedBiexpConst_params_14_bvalues_4tail_MultipleInitializations, 14, 4)
FIT_CONST_MI_IMPLEMENTATION_MARZISEGMENTEDBIEXP(FitMarziSegmentedBiexpConst_params_15_bvalues_4tail_MultipleInitializations, 15, 4)
FIT_CONST_MI_IMPLEMENTATION_MARZISEGMENTEDBIEXP(FitMarziSegmentedBiexpConst_params_16_bvalues_7tail_MultipleInitializations, 16, 7)
FIT_CONST_MI_IMPLEMENTATION_MARZISEGMENTEDBIEXP(FitMarziSegmentedBiexpConst_params_16_bvalues_8tail_MultipleInitializations, 16, 8)
FIT_CONST_MI_IMPLEMENTATION_MARZISEGMENTEDBIEXP(FitMarziSegmentedBiexpConst_params_16_bvalues_9tail_MultipleInitializations, 16, 9)
FIT_CONST_MI_IMPLEMENTATION_CHOSEGMENTEDBIEXP(FitChoSegmentedBiexpConst_params_12_bvalues_3tail_MultipleInitializations, 12, 3)
FIT_CONST_MI_IMPLEMENTATION_CHOSEGMENTEDBIEXP(FitChoSegmentedBiexpConst_params_14_bvalues_3tail_MultipleInitializations, 14, 3)
FIT_CONST_MI_IMPLEMENTATION_CHOSEGMENTEDBIEXP(FitChoSegmentedBiexpConst_params_15_bvalues_3tail_MultipleInitializations, 15, 3)
FIT_CONST_MI_IMPLEMENTATION_CHOSEGMENTEDBIEXP(FitChoSegmentedBiexpConst_params_12_bvalues_4tail_MultipleInitializations, 12, 4)
FIT_CONST_MI_IMPLEMENTATION_CHOSEGMENTEDBIEXP(FitChoSegmentedBiexpConst_params_14_bvalues_4tail_MultipleInitializations, 14, 4)
FIT_CONST_MI_IMPLEMENTATION_CHOSEGMENTEDBIEXP(FitChoSegmentedBiexpConst_params_15_bvalues_4tail_MultipleInitializations, 15, 4)
FIT_CONST_MI_IMPLEMENTATION_CHOSEGMENTEDBIEXP(FitChoSegmentedBiexpConst_params_16_bvalues_7tail_MultipleInitializations, 16, 7)
FIT_CONST_MI_IMPLEMENTATION_CHOSEGMENTEDBIEXP(FitChoSegmentedBiexpConst_params_16_bvalues_8tail_MultipleInitializations, 16, 8)
FIT_CONST_MI_IMPLEMENTATION_CHOSEGMENTEDBIEXP(FitChoSegmentedBiexpConst_params_16_bvalues_9tail_MultipleInitializations, 16, 9)

//define normalized versions for constrained fits
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_1params_2_bvalues_MultipleInitializations, 1, 2)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_1params_3_bvalues_MultipleInitializations, 1, 3)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_1params_4_bvalues_MultipleInitializations, 1, 4)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_1params_5_bvalues_MultipleInitializations, 1, 5)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_1params_6_bvalues_MultipleInitializations, 1, 6)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_1params_7_bvalues_MultipleInitializations, 1, 7)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_1params_8_bvalues_MultipleInitializations, 1, 8)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_1params_9_bvalues_MultipleInitializations, 1, 9)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_1params_10_bvalues_MultipleInitializations, 1, 10)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_1params_11_bvalues_MultipleInitializations, 1, 11)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_1params_12_bvalues_MultipleInitializations, 1, 12)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_1params_14_bvalues_MultipleInitializations, 1, 14)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_1params_15_bvalues_MultipleInitializations, 1, 15)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_2params_2_bvalues_MultipleInitializations, 2, 2)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_2params_3_bvalues_MultipleInitializations, 2, 3)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_2params_4_bvalues_MultipleInitializations, 2, 4)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_2params_5_bvalues_MultipleInitializations, 2, 5)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_2params_6_bvalues_MultipleInitializations, 2, 6)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_2params_7_bvalues_MultipleInitializations, 2, 7)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_2params_8_bvalues_MultipleInitializations, 2, 8)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_2params_9_bvalues_MultipleInitializations, 2, 9)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_2params_10_bvalues_MultipleInitializations, 2, 10)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_2params_11_bvalues_MultipleInitializations, 2, 11)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_2params_12_bvalues_MultipleInitializations, 2, 12)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_2params_14_bvalues_MultipleInitializations, 2, 14)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(FitConstN_2params_15_bvalues_MultipleInitializations, 2, 15)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION_BIEXP(FitBiexpConstN_3params_12_bvalues_MultipleInitializations, 3, 12)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION_BIEXP(FitBiexpConstN_3params_14_bvalues_MultipleInitializations, 3, 14)
FIT_CONST_MI_NORMALIZED_IMPLEMENTATION_BIEXP(FitBiexpConstN_3params_15_bvalues_MultipleInitializations, 3, 15)

//define linearized constrained BFGS fit for different number of parameters in macros, see datatypes.h for implementation
FIT_CONST_MI_LIN_IMPLEMENTATION(Fit_2params_12_bvalues_3tail_MultipleInitializations,2,12)
FIT_CONST_MI_LIN_IMPLEMENTATION(Fit_2params_14_bvalues_3tail_MultipleInitializations,2,14)
FIT_CONST_MI_LIN_IMPLEMENTATION(Fit_2params_15_bvalues_3tail_MultipleInitializations,2,15)
FIT_CONST_MI_SIMPLE_LIN_IMPLEMENTATION(SFit_2params_12_bvalues_3tail_MultipleInitializations,12,3)
FIT_CONST_MI_SIMPLE_LIN_IMPLEMENTATION(SFit_2params_14_bvalues_3tail_MultipleInitializations,14,3)
FIT_CONST_MI_SIMPLE_LIN_IMPLEMENTATION(SFit_2params_15_bvalues_3tail_MultipleInitializations,15,3)
//define asymptotic bi-exponential fits for tail (Pekar, Pang)
FIT_MI_IMPLEMENTATION_ASYMPTOTICBIEXP(AsympFit_2params_12_bvalues_3tail_MultipleInitializations, 12, 3)
FIT_MI_IMPLEMENTATION_ASYMPTOTICBIEXP(AsympFit_2params_14_bvalues_3tail_MultipleInitializations, 14, 3)
FIT_MI_IMPLEMENTATION_ASYMPTOTICBIEXP(AsympFit_2params_15_bvalues_3tail_MultipleInitializations, 15, 3)
//define log-linear bi-exponential fits for tail
FIT_MI_IMPLEMENTATION_LOGLINASYMPTOTICBIEXP(LogLinAsympFit_2params_12_bvalues_2tail_MultipleInitializations, 12, 2)
FIT_MI_IMPLEMENTATION_LOGLINASYMPTOTICBIEXP(LogLinAsympFit_2params_14_bvalues_2tail_MultipleInitializations, 14, 2)
FIT_MI_IMPLEMENTATION_LOGLINASYMPTOTICBIEXP(LogLinAsympFit_2params_15_bvalues_2tail_MultipleInitializations, 15, 2)
FIT_MI_IMPLEMENTATION_LOGLINASYMPTOTICBIEXP(LogLinAsympFit_2params_12_bvalues_3tail_MultipleInitializations, 12, 3)
FIT_MI_IMPLEMENTATION_LOGLINASYMPTOTICBIEXP(LogLinAsympFit_2params_14_bvalues_3tail_MultipleInitializations, 14, 3)
FIT_MI_IMPLEMENTATION_LOGLINASYMPTOTICBIEXP(LogLinAsympFit_2params_15_bvalues_3tail_MultipleInitializations, 15, 3)
// ----------------------------------------------------------------------------------------
//define LM fit for different number of parameters in macros, see datatypes.h for implementation
FIT_LM_MI_IMPLEMENTATION(Fit_2params_2_bvalues_MultipleInitializations,2,2)
FIT_LM_MI_IMPLEMENTATION(Fit_2params_3_bvalues_MultipleInitializations,2,3)
FIT_LM_MI_IMPLEMENTATION(Fit_2params_4_bvalues_MultipleInitializations,2,4)
FIT_LM_MI_IMPLEMENTATION(Fit_2params_5_bvalues_MultipleInitializations,2,5)
FIT_LM_MI_IMPLEMENTATION(Fit_2params_6_bvalues_MultipleInitializations,2,6)
FIT_LM_MI_IMPLEMENTATION(Fit_2params_7_bvalues_MultipleInitializations,2,7)
FIT_LM_MI_IMPLEMENTATION(Fit_2params_8_bvalues_MultipleInitializations,2,8)
FIT_LM_MI_IMPLEMENTATION(Fit_2params_12_bvalues_MultipleInitializations,2,12)
FIT_LM_MI_IMPLEMENTATION(Fit_2params_14_bvalues_MultipleInitializations,2,14)
FIT_LM_MI_IMPLEMENTATION(Fit_2params_15_bvalues_MultipleInitializations,2,15)
FIT_LM_MI_IMPLEMENTATION(Fit_3params_12_bvalues_MultipleInitializations,3,12)
FIT_LM_MI_IMPLEMENTATION(Fit_3params_14_bvalues_MultipleInitializations,3,14)
FIT_LM_MI_IMPLEMENTATION(Fit_3params_15_bvalues_MultipleInitializations,3,15)
FIT_LM_MI_IMPLEMENTATION_BIEXP(FitBiexp_4params_12_bvalues_MultipleInitializations,4,12)
FIT_LM_MI_IMPLEMENTATION_BIEXP(FitBiexp_4params_14_bvalues_MultipleInitializations,4,14)
FIT_LM_MI_IMPLEMENTATION_BIEXP(FitBiexp_4params_15_bvalues_MultipleInitializations,4,15)
FIT_LM_MI_IMPLEMENTATION_MARZISEGMENTEDBIEXP(FitMarziSegmentedBiexp_params_12_bvalues_3tail_MultipleInitializations, 12, 3)
FIT_LM_MI_IMPLEMENTATION_MARZISEGMENTEDBIEXP(FitMarziSegmentedBiexp_params_14_bvalues_3tail_MultipleInitializations, 14, 3)
FIT_LM_MI_IMPLEMENTATION_MARZISEGMENTEDBIEXP(FitMarziSegmentedBiexp_params_15_bvalues_3tail_MultipleInitializations, 15, 3)
FIT_LM_MI_IMPLEMENTATION_MARZISEGMENTEDBIEXP(FitMarziSegmentedBiexp_params_12_bvalues_4tail_MultipleInitializations, 12, 4)
FIT_LM_MI_IMPLEMENTATION_MARZISEGMENTEDBIEXP(FitMarziSegmentedBiexp_params_14_bvalues_4tail_MultipleInitializations, 14, 4)
FIT_LM_MI_IMPLEMENTATION_MARZISEGMENTEDBIEXP(FitMarziSegmentedBiexp_params_15_bvalues_4tail_MultipleInitializations, 15, 4)
// ----------------------------------------------------------------------------------------
// Fast fits for IVIM
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_16_bvalues_15tail_MultipleInitializations_BFGS, 16, 15)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_16_bvalues_14tail_MultipleInitializations_BFGS, 16, 14)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_16_bvalues_13tail_MultipleInitializations_BFGS, 16, 13)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_16_bvalues_12tail_MultipleInitializations_BFGS, 16, 12)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_16_bvalues_11tail_MultipleInitializations_BFGS, 16, 11)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_16_bvalues_10tail_MultipleInitializations_BFGS, 16, 10)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_16_bvalues_9tail_MultipleInitializations_BFGS, 16, 9)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_16_bvalues_8tail_MultipleInitializations_BFGS, 16, 8)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_16_bvalues_7tail_MultipleInitializations_BFGS, 16, 7)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_16_bvalues_6tail_MultipleInitializations_BFGS, 16, 6)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_16_bvalues_5tail_MultipleInitializations_BFGS, 16, 5)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_16_bvalues_4tail_MultipleInitializations_BFGS, 16, 4)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_16_bvalues_3tail_MultipleInitializations_BFGS, 16, 3)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_16_bvalues_2tail_MultipleInitializations_BFGS, 16, 2)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_15_bvalues_14tail_MultipleInitializations_BFGS, 15, 14)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_15_bvalues_13tail_MultipleInitializations_BFGS, 15, 13)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_15_bvalues_12tail_MultipleInitializations_BFGS, 15, 12)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_15_bvalues_11tail_MultipleInitializations_BFGS, 15, 11)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_15_bvalues_10tail_MultipleInitializations_BFGS, 15, 10)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_15_bvalues_9tail_MultipleInitializations_BFGS, 15, 9)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_15_bvalues_8tail_MultipleInitializations_BFGS, 15, 8)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_15_bvalues_7tail_MultipleInitializations_BFGS, 15, 7)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_15_bvalues_6tail_MultipleInitializations_BFGS, 15, 6)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_15_bvalues_5tail_MultipleInitializations_BFGS, 15, 5)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_15_bvalues_4tail_MultipleInitializations_BFGS, 15, 4)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_15_bvalues_3tail_MultipleInitializations_BFGS, 15, 3)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_15_bvalues_2tail_MultipleInitializations_BFGS, 15, 2)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_14_bvalues_13tail_MultipleInitializations_BFGS, 14, 13)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_14_bvalues_12tail_MultipleInitializations_BFGS, 14, 12)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_14_bvalues_11tail_MultipleInitializations_BFGS, 14, 11)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_14_bvalues_10tail_MultipleInitializations_BFGS, 14, 10)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_14_bvalues_9tail_MultipleInitializations_BFGS, 14, 9)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_14_bvalues_8tail_MultipleInitializations_BFGS, 14, 8)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_14_bvalues_7tail_MultipleInitializations_BFGS, 14, 7)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_14_bvalues_6tail_MultipleInitializations_BFGS, 14, 6)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_14_bvalues_5tail_MultipleInitializations_BFGS, 14, 5)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_14_bvalues_4tail_MultipleInitializations_BFGS, 14, 4)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_14_bvalues_3tail_MultipleInitializations_BFGS, 14, 3)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_14_bvalues_2tail_MultipleInitializations_BFGS, 14, 2)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_13_bvalues_12tail_MultipleInitializations_BFGS, 13, 12)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_13_bvalues_11tail_MultipleInitializations_BFGS, 13, 11)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_13_bvalues_10tail_MultipleInitializations_BFGS, 13, 10)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_13_bvalues_9tail_MultipleInitializations_BFGS, 13, 9)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_13_bvalues_8tail_MultipleInitializations_BFGS, 13, 8)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_13_bvalues_7tail_MultipleInitializations_BFGS, 13, 7)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_13_bvalues_6tail_MultipleInitializations_BFGS, 13, 6)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_13_bvalues_5tail_MultipleInitializations_BFGS, 13, 5)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_13_bvalues_4tail_MultipleInitializations_BFGS, 13, 4)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_13_bvalues_3tail_MultipleInitializations_BFGS, 13, 3)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_13_bvalues_2tail_MultipleInitializations_BFGS, 13, 2)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_12_bvalues_11tail_MultipleInitializations_BFGS, 12, 11)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_12_bvalues_10tail_MultipleInitializations_BFGS, 12, 10)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_12_bvalues_9tail_MultipleInitializations_BFGS, 12, 9)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_12_bvalues_8tail_MultipleInitializations_BFGS, 12, 8)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_12_bvalues_7tail_MultipleInitializations_BFGS, 12, 7)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_12_bvalues_6tail_MultipleInitializations_BFGS, 12, 6)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_12_bvalues_5tail_MultipleInitializations_BFGS, 12, 5)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_12_bvalues_4tail_MultipleInitializations_BFGS, 12, 4)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_12_bvalues_3tail_MultipleInitializations_BFGS, 12, 3)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_12_bvalues_2tail_MultipleInitializations_BFGS, 12, 2)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_11_bvalues_10tail_MultipleInitializations_BFGS, 11, 10)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_11_bvalues_9tail_MultipleInitializations_BFGS, 11, 9)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_11_bvalues_8tail_MultipleInitializations_BFGS, 11, 8)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_11_bvalues_7tail_MultipleInitializations_BFGS, 11, 7)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_11_bvalues_6tail_MultipleInitializations_BFGS, 11, 6)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_11_bvalues_5tail_MultipleInitializations_BFGS, 11, 5)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_11_bvalues_4tail_MultipleInitializations_BFGS, 11, 4)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_11_bvalues_3tail_MultipleInitializations_BFGS, 11, 3)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_11_bvalues_2tail_MultipleInitializations_BFGS, 11, 2)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_10_bvalues_9tail_MultipleInitializations_BFGS, 10, 9)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_10_bvalues_8tail_MultipleInitializations_BFGS, 10, 8)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_10_bvalues_7tail_MultipleInitializations_BFGS, 10, 7)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_10_bvalues_6tail_MultipleInitializations_BFGS, 10, 6)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_10_bvalues_5tail_MultipleInitializations_BFGS, 10, 5)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_10_bvalues_4tail_MultipleInitializations_BFGS, 10, 4)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_10_bvalues_3tail_MultipleInitializations_BFGS, 10, 3)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_10_bvalues_2tail_MultipleInitializations_BFGS, 10, 2)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_9_bvalues_8tail_MultipleInitializations_BFGS, 9, 8)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_9_bvalues_7tail_MultipleInitializations_BFGS, 9, 7)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_9_bvalues_6tail_MultipleInitializations_BFGS, 9, 6)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_9_bvalues_5tail_MultipleInitializations_BFGS, 9, 5)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_9_bvalues_4tail_MultipleInitializations_BFGS, 9, 4)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_9_bvalues_3tail_MultipleInitializations_BFGS, 9, 3)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_9_bvalues_2tail_MultipleInitializations_BFGS, 9, 2)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_8_bvalues_7tail_MultipleInitializations_BFGS, 8, 7)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_8_bvalues_6tail_MultipleInitializations_BFGS, 8, 6)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_8_bvalues_5tail_MultipleInitializations_BFGS, 8, 5)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_8_bvalues_4tail_MultipleInitializations_BFGS, 8, 4)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_8_bvalues_3tail_MultipleInitializations_BFGS, 8, 3)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_8_bvalues_2tail_MultipleInitializations_BFGS, 8, 2)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_7_bvalues_6tail_MultipleInitializations_BFGS, 7, 6)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_7_bvalues_5tail_MultipleInitializations_BFGS, 7, 5)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_7_bvalues_4tail_MultipleInitializations_BFGS, 7, 4)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_7_bvalues_3tail_MultipleInitializations_BFGS, 7, 3)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_7_bvalues_2tail_MultipleInitializations_BFGS, 7, 2)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_6_bvalues_5tail_MultipleInitializations_BFGS, 6, 5)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_6_bvalues_4tail_MultipleInitializations_BFGS, 6, 4)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_6_bvalues_3tail_MultipleInitializations_BFGS, 6, 3)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_6_bvalues_2tail_MultipleInitializations_BFGS, 6, 2)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_5_bvalues_4tail_MultipleInitializations_BFGS, 5, 4)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_5_bvalues_3tail_MultipleInitializations_BFGS, 5, 3)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_5_bvalues_2tail_MultipleInitializations_BFGS, 5, 2)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_4_bvalues_3tail_MultipleInitializations_BFGS, 4, 3)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_4_bvalues_2tail_MultipleInitializations_BFGS, 4, 2)
FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(FitLogInitBiexp_3_bvalues_2tail_MultipleInitializations_BFGS, 3, 2)
// ----------------------------------------------------------------------------------------
string resolveMultipleInitializations(const std::vector<std::vector<double> >& limits, std::vector<std::vector<double> >& pinit, const executionInfo& info, int info_i);
template <typename T> std::vector<T> readKey(config_reader& cr, string keyname);
std::vector<std::vector<double> > resolveLimits(const config_reader& cr, executionInfo& info);
// ----------------------------------------------------------------------------------------
std::vector<double>  bvalues;
std::vector<double>  weights;
// ----------------------------------------------------------------------------------------

/**
 * Tokenizes string  
 *
 * @param str string to be tokenize
 * @param tokens (output) vector of tokens
 * @param delimeter delimeter string
 */
void tokenize(const string& str, std::vector<string>& tokens, const string& delimiter = " ")
{
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiter, 0);
  // Find first non-delimiter.
  string::size_type pos = str.find_first_of(delimiter, lastPos);

  while (string::npos != pos || string::npos != lastPos) {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));

    // Skip delimiters.
    lastPos = str.find_first_not_of(delimiter, pos);

    // Find next non-delimiter.
    pos = str.find_first_of(delimiter, lastPos);
  }
}

/**
 * Read Signal Intensity (SI) curves from an ASCII file
 *
 * @see source code for file formatting
 * @param filename ASCII file to be read, with proper formatting
 * @param info (output) file header information, for setting b-value set
 * @return list of SI curves, where each has list of pairs {<b-value(empty)>,<intensity>}
 */
std::vector<std::vector<double> > readSIfile(string filename, executionInfo& info)  {

	std::vector<std::vector<double> > SIs;
	string line;
	string dump;
    int row = 0;
	int col = 0;

	//open file
	ifstream pFile;

	pFile.open(filename.c_str(), std::ios_base::in);
	std::vector<string> tokens;
	if (pFile.is_open())
	{
		bool header_finished = false;
		while(!pFile.eof())
		{
			getline(pFile, line);
			row++;
			//remove [ ] from the string  
			line.erase(std::remove(line.begin(), line.end(), '['), line.end());
			line.erase(std::remove(line.begin(), line.end(), ']'), line.end());
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
            line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
			//tokenize string
			tokens.clear();
			tokenize(line,  tokens, " ");
			if((int)tokens.size() > 0) {
				if(tokens[0].compare("subwindow:") == 0 && (int)tokens.size() > 1) {
					for(int i = 1; i < (int)tokens.size(); i++) {
						info.subwindow.push_back(atoi(tokens[i].c_str()));
					}
					continue;
				}
				if(tokens[0].compare("number:") == 0 && (int)tokens.size() > 1) {
					info.number = atoi(tokens[1].c_str());
					continue;
				}
				if(tokens[0].compare("bset:") == 0 && (int)tokens.size() > 1) {
					for(int i = 1; i < (int)tokens.size(); i++) {
						info.bset.push_back(atof(tokens[i].c_str()));
					}
					continue;
				}
				if(tokens[0].compare("ROIslice:") == 0 && (int)tokens.size() > 1) {
					for(int i = 1; i < (int)tokens.size(); i++) {
						info.ROIslice.push_back(atoi(tokens[i].c_str()));
					}
					continue;
				}
				if(tokens[0].compare("name:") == 0 && (int)tokens.size() > 1) {
					info.ROIName = tokens[1];
					continue;
				}
				if(tokens[0].compare("SIs:") == 0) {
					header_finished = true;
					continue;
				}
				if(!header_finished)
					continue;
				std::vector<double> SI;
				double output;
				for(int i = 0; i < (int)tokens.size(); i++) {
					output = atof(tokens[i].c_str());
					SI.push_back(output);
				}
				SIs.push_back(SI);
			}
		}
		pFile.close();
	}
	else {
		throw string("Unable to open file");
	}

	return SIs;
}

/**
 * Write Signal Intensity (SI) curves to an ASCII file
 *
 * @see source code for file formatting
 * @param filename ASCII file to be written, with proper formatting
 * @param info (output) file header information, for setting b-value set
 * @param SIs curves, where each has list of pairs {<b-value(empty)>,<intensity>}
 */
void writeSIfile(string filename, executionInfo& info, const 	std::vector<std::vector<double> > SIs)  {

	//open file
	ofstream pFile;
	pFile.open(filename.c_str(), std::ios_base::out);
	pFile << std::fixed;
	std::vector<string> tokens;
    if (pFile.is_open())
    {
		//write header lines 
		pFile << "subwindow: [";
		for(int i = 0; i < (int)info.subwindow.size(); i++)
			pFile << info.subwindow[i] << " ";
		pFile << "]" << endl;
		pFile << "number: " << info.number << endl;
		pFile << "bset: [";
		for(int i = 0; i < (int)info.bset.size(); i++)
			pFile << info.bset[i] << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed;
		pFile << "]" << endl;
		pFile << "ROIslice: [";
		for(int i = 0; i < (int)info.ROIslice.size(); i++)
			pFile << info.ROIslice[i] << " ";
		pFile << "]" << endl;
		pFile << "name: " << info.ROIName << endl;
		pFile << "executiontime: " << info.executiontime << " seconds"<< endl;
		pFile << "description:" << info.desc << endl;
		//write parameter values, including RMSE
		pFile << "parameters: ";
		for(int i = 0; i < (int)info.parameterNames.size(); i++)
			pFile << info.parameterNames[i] << " ";
		pFile << endl;
		pFile << "SIs:" << endl;
		for(int i = 0; i < (int)SIs.size(); i++) {
			for(int j = 0; j < (int)SIs[i].size(); j++) {
				pFile << SIs[i][j] << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed;
			}						
			pFile << endl;
		}
        pFile.close();
    }
    else
        cout << "Unable to open file " << filename << " for writing"; 
}

/**
 * Writes curves to a file
 *
 * @see ModelName for model code numbers 
 * @param filename file where to write
 * @param p_vals list of pairs {<parameter vector>,<RMSE value>}
 * @param info header information to be written to the file
 */
void writefile(string filename, std::vector<std::vector<double> > p_vals, executionInfo& info) {

	//open file
	ofstream pFile;
	pFile.open(filename.c_str(), std::ios_base::out);
	pFile << std::fixed;
	std::vector<string> tokens;
    if (pFile.is_open())
    {
		//write header lines 
		pFile << "subwindow: [";
		for(int i = 0; i < (int)info.subwindow.size(); i++)
			pFile << info.subwindow[i] << " ";
		pFile << "]" << endl;
		pFile << "number: " << info.number << endl;
		pFile << "bset: [";
		for(int i = 0; i < (int)info.bset.size(); i++)
			pFile << info.bset[i] << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed;
		pFile << "]" << endl;
		pFile << "ROIslice: [";
		for(int i = 0; i < (int)info.ROIslice.size(); i++)
			pFile << info.ROIslice[i] << " ";
		pFile << "]" << endl;
		pFile << "name: " << info.ROIName << endl;
		pFile << "executiontime: " << info.executiontime << " seconds"<< endl;
		pFile << "description:" << info.desc << endl;

		//write parameter values, including RMSE
		pFile << "parameters: ";
		for(int i = 0; i < (int)info.parameterNames.size(); i++)
			pFile << info.parameterNames[i] << " ";
		pFile << "RMSE" << endl;
		for(int i = 0; i < (int)p_vals.size(); i++) {
			for(int j = 0; j < (int)p_vals[i].size(); j++) {
				pFile << p_vals[i][j] << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed;
			}						
			pFile << endl;
		}
        pFile.close();
    }
    else
        cout << "Unable to open file " << filename << " for writing"; 
}

/**
 * Prints usage information to stdout
 *
 * @param programname
 */
void printUsage(string programname) {
	cout << "DWI B-value fitting " << programname << endl;
	cout << "Compilation infostring:" << getCompilationInfoString() << endl;
	cout << endl;
	cout << "usage: <executable>.exe <ASCII file> <model name> <[Weighted]|BFGS|LM | SIMULATION> [settings filename] [<SD> <iterations> [<param1>..<paramN>]] [No tail values]" << endl;
	cout << endl;
	cout << "<ASCII file>:\tfile from where SI curves are located" << endl;
	cout << endl;
	cout << "<model name>: either of:\t" << endl;
	cout << "\tMono - Monoexponential" << endl;
	cout << "\tT1rho - T1rho monoexponential" << endl;
	cout << "\tKurt - Kurtosis model" << endl;
	cout << "\tStretched - Stretched exponential fitting function" << endl;
	cout << "\tBiexp - Bi-exponential model" << endl;
	cout << "\tMarziSegmentedBiexp - Fitting tail first, then f, Df together" << endl;
	cout << "\tChoSegmentedBiexp - Fitting tail first, then asymptotic f, and Df fitted alone" << endl;
	cout << "\tAsymptoticBiexp" << endl;
	cout << "\tLogLinAsymptoticBiexp" << endl;
	cout << "\tNonGaussianIVIM" << endl;
	cout << "\tLinNonGaussianIVIM" << endl;
	cout << "\tSLinNonGaussianIVIM" << endl;
	cout << "\tLogInitBiexp - Asymptotic and log-linear slope for initialization of f and Df" << endl;
	cout << endl;
	cout << "Weighted: use constant weightings for b-values' LSQ value in fittings. Use as prefix, e. g. 'WeightedLM'" << endl;
	cout << endl;
	cout << "BFGS|LM|SIMULATION: optimization methods" << endl;
	cout << "\tLM: Levenberg-Marquardt algorithm" << endl;
	cout << "\tBFGS: Broyden Fletcher Goldfarb Shanno algorithm" << endl;
	cout << "\tSIMULATION: SI curve simulation. Reads additional parameters from 'Simulation' in settings file" << endl;
	cout << endl;
	cout << "settings filename: filename of settings in ASCII format, default is 'settings.ini'" << endl;
	cout << endl;
	cout << "SD:\tstandard deviation of Rician noise" << endl;
	cout << endl;
	cout << "iterations:\tnumber of noise realizations" << endl;
	cout << endl;
	cout << "param1..N:\tsimulation parameter values that are used instead of ones found in the settings file,\n\t all must be given or none for reading values from settings" << endl;
	cout << endl;
	cout << "SD:\tNumber of tail values to use in segmented fits, mandatory for:" << endl;
	cout << "\tMarziSegmentedBiexp" << endl;
	cout << "\tChoSegmentedBiexp" << endl;
	cout << "\tAsymptoticBiexp" << endl;
	cout << "\tLogLinAsymptoticBiexp" << endl;
	cout << endl;
}

// ----------------------------------------------------------------------------------------
/**
 * Main function: 
 * 1. Reads SI data from ASCII file
 * 2. Makes model fitting to all SI curves
 * 3. Writes fitting results to ASCII file
 *
 * @param argc number of arguments
 * @param argv arguments as character arrays
 */
int main(int argc,char *argv[])
{
	//true for weighted version of fitting
	bool weighted = false;
	//number of fitting method, (default) 1 == Traditional Levenberg-Marquardt, 2 == Trust-Region Reflective
	int fitting_method = OP_LM;
	//SI data
	std::vector<std::vector<double> > SIdata;
    try
    {
		if(argc < 3) {
			printUsage(argv[0]);
			throw "not enough input parameters was given";
		}
		string mfilename = argv[0];
		mfilename = mfilename + getCompilationInfoString();
		string filename = argv[1];
		executionInfo info;
		executionInfo simulation_info;
		info.modelName = argv[2];
		cout << "Modelname:" << info.modelName << endl;

		//resolve settings
		string sfilename;
		if(argc > 4) {
			sfilename = argv[4];
		} else {
			sfilename = "settings.ini";
		}
		cout << "reading " << sfilename << " for settings"<< endl;
		config_reader cr(sfilename);

		//set b-values to all curves that we are going to fit
		cout << "resolving b-values" << endl;
		bvalues = readKey<double>(cr, "bvalues");
		int no_bvalues = (int)bvalues.size();
		cout << no_bvalues << " b-values found" << endl;
		std::vector<double> indexes;
		indexes = bvalues;	

		//resolve fitting method
		#pragma region resolve fitting method
		if(argc > 3) {
			string fitmethodstr = argv[3];
			//resolve weightings
			if(fitmethodstr.find("Weighted") != string::npos) {
				weighted = true;
				if(cr.is_key_defined("weights")) {
					weights = readKey<double>(cr, "weights");
					cout << weights.size() << " weights found" << endl;
					indexes = std::vector<double>(no_bvalues);			
					cout << "indexes=[";
					for(unsigned int i = 0; i < indexes.size(); i++) {
						indexes[i] = i;
						cout << i << " ";
					}
					cout << "]" << endl;
				} else {
					weighted = false;
					std::fill(weights.begin(), weights.end(), 1.0);
				}
			}
			//resolve fitting method
            //Leveberg-Marquard, normalized data
			if(fitmethodstr.find("LM_NORMALIZED") != string::npos) {
				fitting_method = OP_LM_NORMALIZED;
			}
			//Broyden-Fletcher.Goldfarb-Shanno algorithm
			if(fitmethodstr.find("BFGS") != string::npos) {
				fitting_method = OP_BFGS;
			}
			//Broyden-Fletcher.Goldfarb-Shanno algorithm, normalized
			if(fitmethodstr.find("BFGS_NORMALIZED") != string::npos) {
				fitting_method = OP_BFGS_NORMALIZED;
			}
			if(fitmethodstr.find("SIMULATION") != string::npos) {
				fitting_method = OP_SIMULATION;
			}
		}
		if(weighted && (int)weights.size() != no_bvalues) {
			cout << "weights array length " << (int)SIdata.size() << " in settings";
			cout << " is not consistent with b-value array length " << no_bvalues << " in settings" << endl;
			return 1;
		}
		#pragma endregion

		//read signal intensity data
		#pragma region read signal intensity data
		if(fitting_method != OP_SIMULATION) {
			try {
				cout << "Opening " << filename << " .. ";
				SIdata = readSIfile(filename, info);
			} catch(string str) {
				throw str;
			} catch (...) { 
				cout << "unknown exception while reading SI file";
				return 1;
			}
			cout << (int)SIdata.size() << " SI curves found" << endl;
			if((int)SIdata.size() == 0) {
				cout << (int)SIdata.size() << " data was not found for processing" << endl;			
				return 1;
			}
			if(weighted && (int)SIdata[0].size() != (int)weights.size()) {
				cout << "SI curve length " << (int)SIdata[0].size();
				cout << " is not consistent with weights array length " << (int)weights.size() << " in settings" << endl;			
				return 1;
			}
			if((int)SIdata[0].size() != no_bvalues) {
				cout << "SI curve length " << (int)SIdata[0].size();
				cout << " is not consistent with b-value array length " << no_bvalues << " in settings" << endl;			
				return 1;
			}
			for(int b_i = 0; b_i < no_bvalues; b_i++) {
				indexes[b_i] = info.bset[b_i];
				bvalues[b_i] = info.bset[b_i];
			}
		//place default information to info
		} else {
			info.subwindow.push_back(0);
			info.subwindow.push_back(0);
			info.subwindow.push_back(0);
			info.subwindow.push_back(0);
			info.number = 0;
			simulation_info.subwindow.push_back(0);
			simulation_info.subwindow.push_back(0);
			simulation_info.subwindow.push_back(0);
			simulation_info.subwindow.push_back(0);
			simulation_info.number = 0;
			for(int b_i = 0; b_i < no_bvalues; b_i++) {
				info.bset.push_back(bvalues[b_i]);
				simulation_info.bset.push_back(bvalues[b_i]);
			}
			info.ROIslice.push_back(0);
			info.ROIName = "simulation";
			simulation_info.ROIslice.push_back(0);
			simulation_info.ROIName = "simulation";
			simulation_info.modelName = "Simulation";
		}
		#pragma endregion
		
		std::vector<std::vector<double> > limits;
		std::vector<std::vector<double> > limits_1st;
		std::vector<std::vector<double> > limits_2nd;
		std::vector<double> C_values;
		//simulation limits only used in simulations
		std::vector<std::vector<double> > simulation_limits;
		//number of tail values used only in specific models
		int No_tail_values = 3;
		//string containing additional information about the fit 
		std::ostringstream os;
		//resolve parameter limits for C to be percentages of maximum intensity
		#pragma region resolve parameter limits for C to be percentages of maximum intensity
		try {
			limits = resolveLimits(cr, info);
			if((int)info.parameterNames.size() == 0) {
				THROW_FATAL("No parameter names found for model")
			}
            VERBOSE2("last parameter " << info.parameterNames[(int)info.parameterNames.size()-1])
			if(info.parameterNames[(int)info.parameterNames.size()-1].compare("C") == 0) {
				//convert percentage values to coefficients for C, use 1.0 if none declared
                VERBOSE2("limits.size() " << limits.size())
                VERBOSE2("limits[limits.size()-1][0] " << limits[limits.size()-1][0])
                VERBOSE2("limits[limits.size()-1][1] " << limits[limits.size()-1][1])
                VERBOSE2("limits[limits.size()-1][2] " << limits[limits.size()-1][2])
				for(double C = limits[limits.size()-1][0]; C <= limits[limits.size()-1][2]; C += limits[limits.size()-1][1]) {
                    VERBOSE2("C_values.push_back(" << C/100.0 << ")")
					C_values.push_back(C/100.0);
				}
                if(C_values.size() == 0) {
                    C_values.push_back(limits[limits.size()-1][0]/100.0);
                }
			} else {
				C_values.push_back(1.0);
			}
			os << "C=[" << C_values[0] << ".." << C_values[C_values.size()-1] << "]";
			if(fitting_method == OP_SIMULATION) {
				//Read 'Simulation' section from settings file
				if(argc > 6) {
					//SD
					std::vector<double> limSD(3);
					limSD[0] = atof(argv[5]);
					limSD[1] = atof(argv[5]);
					limSD[2] = atof(argv[5]);
					simulation_limits.push_back(limSD);
					simulation_info.parameterNames.push_back("SD");
					simulation_info.parameterLimit_explicit.push_back(false);
					cout << "simulation noise SD " << limSD[0] << endl;
					//iterations
					std::vector<double> limIT(3);
					limIT[0] = atof(argv[6]);
					limIT[1] = atof(argv[6]);
					limIT[2] = atof(argv[6]);
					simulation_limits.push_back(limIT);
					simulation_info.parameterNames.push_back("iterations");
					simulation_info.parameterLimit_explicit.push_back(false);
					cout << "simulation iterations for each noise level " << limIT[0] << endl;
					//read simulation parameters from command line
					if(argc > 7) {
						cout << "looking for " << ((int)limits.size()) << " parameters for " << info.modelName << endl;
						if(argc > 6+(int)limits.size()) {
							for(unsigned long l_i = 0; l_i < limits.size(); l_i++) {
								limits[l_i][0] = atof(argv[7+l_i]);
								limits[l_i][1] = atof(argv[7+l_i]);
								limits[l_i][2] = atof(argv[7+l_i]);
								cout << "simulating with " << info.parameterNames[l_i] << "=" << limits[l_i][0] << endl;
							}
						} else {
							THROW_FATAL("simution parameters were not properly given in command line")
						}
					}
				} else {
					THROW_FATAL("SD and number of iterations must be given as parameters to simulations")
				}
			} else {
                VERBOSE2("resolving tail values")
				if((info.modelName.compare(MARZISEGMENTEDBIEXP) == 0) || 
				   (info.modelName.compare(CHOSEGMENTEDBIEXP) == 0) || 
				   (info.modelName.compare(ASYMPTOTICBIEXP) == 0) ||
				   (info.modelName.compare(LOGINITIVIM) == 0) ||
				   (info.modelName.compare(LOGLINASYMPTOTICBIEXP) == 0)) {
					if(argc > 5) {
						No_tail_values = (int)atof(argv[5]);
					} else {
						cout << "Number of tail values was not given for model:" << info.modelName << " using default 3" << endl;
					}
					os << " [" << No_tail_values << " tail values used in fitting]";
				}
                VERBOSE2("resolving tail values done")
			}

		} catch(string s) {
			cout << "exception while resolving parameter limits:" << s << endl;
			return 1;
		} catch(const char* s) {
			cout << "exception while resolving parameter limits:" << s << endl;
			return 1;
		} catch(...) {
			cout << "unknown exception while resolving parameter limits" << endl;
			return 1;
		}
		#pragma endregion

		//resolve base output name
		std::vector<string> tokens;
		tokenize(filename,  tokens, ".");
		string output_filename;
		for(int i = 0; i < (int)tokens.size()-1; i++)
			output_filename = output_filename + tokens[i];

		//run fittings 
		std::vector<std::vector<double> > pinit;
		std::vector<std::vector<double> > pinit2;
		std::vector<std::vector<double> > simulation_pinit;
		std::vector<std::vector<double> > p_vals;
		if(info.modelName.compare(MONO) == 0) {
			#pragma region MONO
			matrix<double,2,1> x_lo;
			matrix<double,2,1> x_hi;
			//ADCm
			x_lo(0) = 0;
			x_hi(0) = 1;
			//C (in normalized intensity scale)
			x_lo(1) = 0;
			x_hi(1) = 10;
			switch(fitting_method) {
				case OP_BFGS:
					info.desc = mfilename + " dlib constrained BFGS technique MI " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();
					cout << info.desc << endl;
					switch(no_bvalues) {
						case 2:
							FitConst_2params_2_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponential_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 3:
							FitConst_2params_3_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponential_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 4:
							FitConst_2params_4_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponential_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 5:
							FitConst_2params_5_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponential_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 6:
							FitConst_2params_6_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponential_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 7:
							FitConst_2params_7_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponential_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 8:
							FitConst_2params_8_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponential_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 9:
							FitConst_2params_9_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponential_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 10:
							FitConst_2params_10_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponential_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 11:
							FitConst_2params_11_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponential_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 12:
							FitConst_2params_12_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponential_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 13:
							FitConst_2params_12_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponential_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 14:
							FitConst_2params_14_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponential_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 15:
							FitConst_2params_15_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponential_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						default:
							THROW_FATAL("Unsupported number of b-values")
					}
					break;
				case OP_BFGS_NORMALIZED:
					info.desc = mfilename + " dlib constrained BFGS technique MI, normalized SI data " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();
					cout << info.desc << endl;
					switch(no_bvalues) {
						case 2:
							FitConstN_1params_2_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponentialN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 3:
							FitConstN_1params_3_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponentialN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 4:
							FitConstN_1params_4_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponentialN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 5:
							FitConstN_1params_5_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponentialN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 6:
							FitConstN_1params_6_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponentialN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 7:
							FitConstN_1params_7_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponentialN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 8:
							FitConstN_1params_8_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponentialN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 9:
							FitConstN_1params_9_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponentialN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 10:
							FitConstN_1params_10_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponentialN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 11:
							FitConstN_1params_11_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponentialN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 12:
							FitConstN_1params_12_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponentialN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 14:
							FitConstN_1params_14_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponentialN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 15:
							FitConstN_1params_15_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoExponentialN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						default:
							THROW_FATAL("Unsupported number of b-values")
					}
					break;
				case OP_LM:
					double (*residual_fun)(const std::pair<double, double>&, const matrix<double,2,1>&);
					if(weighted) {
						info.desc = mfilename + " dlib Traditional Levenberg-Marquardt technique, weighted " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();
						residual_fun = MonoExponential_residual_weights;
					} else {
						info.desc = mfilename + " dlib Traditional Levenberg-Marquardt technique " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();			
						residual_fun = MonoExponential_residual;
					}
					cout << info.desc << endl;
					switch(no_bvalues) {
						case 12:
							Fit_2params_12_bvalues_MultipleInitializations(indexes, SIdata, pinit, residual_fun, p_vals, info, C_values);
							break;
						case 14:
							Fit_2params_14_bvalues_MultipleInitializations(indexes, SIdata, pinit, residual_fun, p_vals, info, C_values);
							break;
						case 15:
							Fit_2params_15_bvalues_MultipleInitializations(indexes, SIdata, pinit, residual_fun, p_vals, info, C_values);
							break;
						default:
							THROW_FATAL("Unsupported number of b-values")
					}
					break;
				default:
					THROW_FATAL("Unsupported fitting method")
			}
			#pragma endregion
		} else if(info.modelName.compare(T1RHO) == 0) {
			#pragma region T1RHO
			matrix<double,2,1> x_lo;
			matrix<double,2,1> x_hi;
			//rho
			x_lo(0) = 0;
			x_hi(0) = 10e8;
			//C (in normalized intensity scale)
			x_lo(1) = 0;
			x_hi(1) = 10;
			switch(fitting_method) {
			case OP_BFGS:
					info.desc = mfilename + " dlib constrained BFGS technique MI " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();
					cout << info.desc << endl;
					switch(no_bvalues) {
						case 2:
							FitConst_2params_2_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoT1rho_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 3:
							FitConst_2params_3_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoT1rho_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 4:
							FitConst_2params_4_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoT1rho_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 5:
							FitConst_2params_5_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoT1rho_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 6:
							FitConst_2params_6_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoT1rho_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 7:
							FitConst_2params_7_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoT1rho_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 8:
							FitConst_2params_8_bvalues_MultipleInitializations(indexes, SIdata, pinit, MonoT1rho_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						default:
							THROW_FATAL("Unsupported number of b-values")
					}
					break;
				case OP_LM:
					double (*residual_fun)(const std::pair<double, double>&, const matrix<double,2,1>&);
					if(weighted) {
						info.desc = mfilename + " dlib Traditional Levenberg-Marquardt technique, weighted " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();
						residual_fun = MonoT1rho_residual_weights;
					} else {
						info.desc = mfilename + " dlib Traditional Levenberg-Marquardt technique " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();			
						residual_fun = MonoT1rho_residual;
					}
					cout << info.desc << endl;
					switch(no_bvalues) {
						case 2:
							Fit_2params_2_bvalues_MultipleInitializations(indexes, SIdata, pinit, residual_fun, p_vals, info, C_values);
							break;
						case 3:
							Fit_2params_3_bvalues_MultipleInitializations(indexes, SIdata, pinit, residual_fun, p_vals, info, C_values);
							break;
						case 4:
							Fit_2params_4_bvalues_MultipleInitializations(indexes, SIdata, pinit, residual_fun, p_vals, info, C_values);
							break;
						case 5:
							Fit_2params_5_bvalues_MultipleInitializations(indexes, SIdata, pinit, residual_fun, p_vals, info, C_values);
							break;
						case 6:
							Fit_2params_6_bvalues_MultipleInitializations(indexes, SIdata, pinit, residual_fun, p_vals, info, C_values);
							break;
						case 7:
							Fit_2params_7_bvalues_MultipleInitializations(indexes, SIdata, pinit, residual_fun, p_vals, info, C_values);
							break;
						case 8:
							Fit_2params_8_bvalues_MultipleInitializations(indexes, SIdata, pinit, residual_fun, p_vals, info, C_values);
							break;
						default:
							THROW_FATAL("Unsupported number of b-values")
					}
					break;
				default:
					THROW_FATAL("Unsupported fitting method")
			}
			#pragma endregion
		} else if(info.modelName.compare(KURT) == 0) {
			#pragma region KURT
			//un-linearized version of Kurtosis model
			matrix<double,3,1> x_lo;
			matrix<double,3,1> x_hi;
			//ADCk
			x_lo(0) = 0;
			x_hi(0) = 1;
			//K
			x_lo(1) = 0;
			x_hi(1) = 10;
			//C (in normalized intensity scale)
			x_lo(2) = 0;
			x_hi(2) = 10;
			//calibration with largets b-value in order to avoid hitting maximum of double precision
			//largest allowed exponent is ~700 for keeping resulting numbers below 1.7977e+308 
			//C is assumed to be in [0..2] in normalized data
			double largest_bvalue = indexes[no_bvalues-1];
			double largest_ADCk = x_hi(0)*1.01;
			double largest_K = x_hi(1)*1.01;
			double largest_C = 2;
			double double_max = std::numeric_limits<double>::max();
			double exp_max = std::log(double_max);
			//resolve largest x in ax^2 + bx + c == 0
			double a_hi = (largest_ADCk*largest_ADCk/6.0)*largest_K;
			double b_hi = -largest_ADCk;
			double c_hi = -exp_max;
			double largest_allowed_bvalue = (-b_hi+std::sqrt(b_hi*b_hi-4.0*a_hi*c_hi))/(2.0*a_hi);
			//divide by 2 to account for calculating squared error with maximum output
			largest_allowed_bvalue /= 2;
			//subtraction to accound for calculating squared error with maximum C, and for maximum value range in double
			largest_allowed_bvalue -= std::log(largest_C*largest_C)+std::log(0.5);
			double coeff = 1.0;
			if(largest_bvalue > largest_allowed_bvalue) {
				cout << "Values are re-scaled in order to keep inside numerical limits of the implementation" << endl;
				coeff = largest_allowed_bvalue/largest_bvalue;
				for(int b_i = 0; b_i < no_bvalues; b_i++) {
					indexes[b_i] *= coeff;
					for(unsigned long s_i = 0; s_i < SIdata.size(); s_i++) {
						SIdata[s_i][b_i] *= coeff;
					}
				}
			}
			switch(fitting_method) {
				case OP_BFGS:
					info.desc = mfilename + " dlib constrained BFGS technique MI " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();
					cout << info.desc << endl;
					switch(no_bvalues) {
						case 2:
							FitConst_3params_2_bvalues_MultipleInitializations(indexes, SIdata, pinit, Kurt_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 3:
							FitConst_3params_3_bvalues_MultipleInitializations(indexes, SIdata, pinit, Kurt_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 4:
							FitConst_3params_4_bvalues_MultipleInitializations(indexes, SIdata, pinit, Kurt_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 5:
							FitConst_3params_5_bvalues_MultipleInitializations(indexes, SIdata, pinit, Kurt_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 6:
							FitConst_3params_6_bvalues_MultipleInitializations(indexes, SIdata, pinit, Kurt_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 7:
							FitConst_3params_7_bvalues_MultipleInitializations(indexes, SIdata, pinit, Kurt_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 8:
							FitConst_3params_8_bvalues_MultipleInitializations(indexes, SIdata, pinit, Kurt_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 9:
							FitConst_3params_9_bvalues_MultipleInitializations(indexes, SIdata, pinit, Kurt_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 10:
							FitConst_3params_10_bvalues_MultipleInitializations(indexes, SIdata, pinit, Kurt_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 11:
							FitConst_3params_11_bvalues_MultipleInitializations(indexes, SIdata, pinit, Kurt_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 12:
							FitConst_3params_12_bvalues_MultipleInitializations(indexes, SIdata, pinit, Kurt_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 14:
							FitConst_3params_14_bvalues_MultipleInitializations(indexes, SIdata, pinit, Kurt_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 15:
							FitConst_3params_15_bvalues_MultipleInitializations(indexes, SIdata, pinit, Kurt_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						default:
							THROW_FATAL("Unsupported number of b-values")
					}
					break;
				case OP_BFGS_NORMALIZED:
					info.desc = mfilename + " dlib constrained BFGS technique MI, normalized SI data " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();
					cout << info.desc << endl;
					switch(no_bvalues) {
						case 2:
							FitConstN_2params_2_bvalues_MultipleInitializations(indexes, SIdata, pinit, KurtN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 3:
							FitConstN_2params_3_bvalues_MultipleInitializations(indexes, SIdata, pinit, KurtN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 4:
							FitConstN_2params_4_bvalues_MultipleInitializations(indexes, SIdata, pinit, KurtN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 5:
							FitConstN_2params_5_bvalues_MultipleInitializations(indexes, SIdata, pinit, KurtN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 6:
							FitConstN_2params_6_bvalues_MultipleInitializations(indexes, SIdata, pinit, KurtN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 7:
							FitConstN_2params_7_bvalues_MultipleInitializations(indexes, SIdata, pinit, KurtN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 8:
							FitConstN_2params_8_bvalues_MultipleInitializations(indexes, SIdata, pinit, KurtN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 9:
							FitConstN_2params_9_bvalues_MultipleInitializations(indexes, SIdata, pinit, KurtN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 10:
							FitConstN_2params_10_bvalues_MultipleInitializations(indexes, SIdata, pinit, KurtN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 11:
							FitConstN_2params_11_bvalues_MultipleInitializations(indexes, SIdata, pinit, KurtN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 12:
							FitConstN_2params_12_bvalues_MultipleInitializations(indexes, SIdata, pinit, KurtN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 14:
							FitConstN_2params_14_bvalues_MultipleInitializations(indexes, SIdata, pinit, KurtN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 15:
							FitConstN_2params_15_bvalues_MultipleInitializations(indexes, SIdata, pinit, KurtN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						default:
							THROW_FATAL("Unsupported number of b-values")
					}
					break;
				case OP_LM:
					double (*residual_fun)(const std::pair<double, double>&, const matrix<double,3,1>&);
					if(weighted) {
						info.desc = mfilename + " dlib Traditional Levenberg-Marquardt technique MI, weighted " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();
						residual_fun = Kurt_residual_weights;
					} else {
						info.desc = mfilename + " dlib Traditional Levenberg-Marquardt technique MI " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();			
						residual_fun = Kurt_residual;
					}
					cout << info.desc << endl;
					switch(no_bvalues) {
						case 12:
							Fit_3params_12_bvalues_MultipleInitializations(indexes, SIdata, pinit, residual_fun, p_vals, info, C_values);
							break;
						case 14:
							Fit_3params_14_bvalues_MultipleInitializations(indexes, SIdata, pinit, residual_fun, p_vals, info, C_values);
							break;
						case 15:
							Fit_3params_15_bvalues_MultipleInitializations(indexes, SIdata, pinit, residual_fun, p_vals, info, C_values);
							break;
						default:
							THROW_FATAL("Unsupported number of b-values")
					}
					break;
				default:
					THROW_FATAL("Unsupported fitting method")
			}
			//calibration of parameter values after fittings with b-value normalized data
			for(unsigned long p_i = 0; p_i < p_vals.size(); p_i++) {
				if(fitting_method == OP_BFGS_NORMALIZED) {
					//ADCk
					p_vals[p_i][0] = p_vals[p_i][0]*coeff;
					//RMSE
					p_vals[p_i][2] = p_vals[p_i][2]/coeff;
				} else {
					//ADCk
					p_vals[p_i][0] = p_vals[p_i][0]*coeff;
					//C
					p_vals[p_i][2] = p_vals[p_i][2]/coeff;
					//RMSE
					p_vals[p_i][3] = p_vals[p_i][3]/coeff;
				}
			}
			#pragma endregion
		} else if(info.modelName.compare(STRETCHED) == 0) {
			#pragma region STRETCHED
			matrix<double,3,1> x_lo;
			matrix<double,3,1> x_hi;
			//ADCs
			x_lo(0) = 0;
			x_hi(0) = 1;
			//Alpha
			x_lo(1) = 0;
			x_hi(1) = 1;
			//C (in normalized intensity scale)
			x_lo(2) = 0;
			x_hi(2) = 10;
			switch(fitting_method) {
				case OP_BFGS:
						info.desc = mfilename + " dlib constrained BFGS technique MI " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();
						cout << info.desc << endl;
						switch(no_bvalues) {
							case 2:
								FitConst_3params_2_bvalues_MultipleInitializations(indexes, SIdata, pinit, Stretched_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
								break;
							case 3:
								FitConst_3params_3_bvalues_MultipleInitializations(indexes, SIdata, pinit, Stretched_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
								break;
							case 4:
								FitConst_3params_4_bvalues_MultipleInitializations(indexes, SIdata, pinit, Stretched_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
								break;
							case 5:
								FitConst_3params_5_bvalues_MultipleInitializations(indexes, SIdata, pinit, Stretched_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
								break;
							case 6:
								FitConst_3params_6_bvalues_MultipleInitializations(indexes, SIdata, pinit, Stretched_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
								break;
							case 7:
								FitConst_3params_7_bvalues_MultipleInitializations(indexes, SIdata, pinit, Stretched_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
								break;
							case 8:
								FitConst_3params_8_bvalues_MultipleInitializations(indexes, SIdata, pinit, Stretched_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
								break;
							case 9:
								FitConst_3params_9_bvalues_MultipleInitializations(indexes, SIdata, pinit, Stretched_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
								break;
							case 10:
								FitConst_3params_10_bvalues_MultipleInitializations(indexes, SIdata, pinit, Stretched_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
								break;
							case 11:
								FitConst_3params_11_bvalues_MultipleInitializations(indexes, SIdata, pinit, Stretched_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
								break;
							case 12:
								FitConst_3params_12_bvalues_MultipleInitializations(indexes, SIdata, pinit, Stretched_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
								break;
							case 14:
								FitConst_3params_14_bvalues_MultipleInitializations(indexes, SIdata, pinit, Stretched_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
								break;
							case 15:
								FitConst_3params_15_bvalues_MultipleInitializations(indexes, SIdata, pinit, Stretched_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
								break;
							default:
								THROW_FATAL("Unsupported number of b-values")
						}
						break;
				case OP_BFGS_NORMALIZED:
					info.desc = mfilename + " dlib constrained BFGS technique MI, normalized SI data " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();
					cout << info.desc << endl;
					switch(no_bvalues) {
						case 2:
							FitConstN_2params_2_bvalues_MultipleInitializations(indexes, SIdata, pinit, StretchedN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 3:
							FitConstN_2params_3_bvalues_MultipleInitializations(indexes, SIdata, pinit, StretchedN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 4:
							FitConstN_2params_4_bvalues_MultipleInitializations(indexes, SIdata, pinit, StretchedN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 5:
							FitConstN_2params_5_bvalues_MultipleInitializations(indexes, SIdata, pinit, StretchedN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 6:
							FitConstN_2params_6_bvalues_MultipleInitializations(indexes, SIdata, pinit, StretchedN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 7:
							FitConstN_2params_7_bvalues_MultipleInitializations(indexes, SIdata, pinit, StretchedN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 8:
							FitConstN_2params_8_bvalues_MultipleInitializations(indexes, SIdata, pinit, StretchedN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 9:
							FitConstN_2params_9_bvalues_MultipleInitializations(indexes, SIdata, pinit, StretchedN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 10:
							FitConstN_2params_10_bvalues_MultipleInitializations(indexes, SIdata, pinit, StretchedN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 11:
							FitConstN_2params_11_bvalues_MultipleInitializations(indexes, SIdata, pinit, StretchedN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 12:
							FitConstN_2params_12_bvalues_MultipleInitializations(indexes, SIdata, pinit, StretchedN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 14:
							FitConstN_2params_14_bvalues_MultipleInitializations(indexes, SIdata, pinit, StretchedN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 15:
							FitConstN_2params_15_bvalues_MultipleInitializations(indexes, SIdata, pinit, StretchedN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						default:
							THROW_FATAL("Unsupported number of b-values")
					}
					break;
				case OP_LM:
					double (*residual_fun)(const std::pair<double, double>&, const matrix<double,3,1>&);
					if(weighted) {
						info.desc = mfilename + " dlib Traditional Levenberg-Marquardt technique MI, weighted " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();
						residual_fun = Stretched_residual_weights;
					} else {
						info.desc = mfilename + " dlib Traditional Levenberg-Marquardt technique MI " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();			
						residual_fun = Stretched_residual;
					}
					cout << info.desc << endl;
					switch(no_bvalues) {
						case 12:
							Fit_3params_12_bvalues_MultipleInitializations(indexes, SIdata, pinit, residual_fun, p_vals, info, C_values);
							break;
						case 14:
							Fit_3params_14_bvalues_MultipleInitializations(indexes, SIdata, pinit, residual_fun, p_vals, info, C_values);
							break;
						case 15:
							Fit_3params_15_bvalues_MultipleInitializations(indexes, SIdata, pinit, residual_fun, p_vals, info, C_values);
							break;
						default:
							THROW_FATAL("Unsupported number of b-values")
					}
					break;
				default:
					THROW_FATAL("Unsupported fitting method")
			}
			#pragma endregion
		} else if(info.modelName.compare(BIEXP) == 0) {
			#pragma region BIEXP
			matrix<double,4,1> x_lo;
			matrix<double,4,1> x_hi;
			//f
			x_lo(0) = 0;
			x_hi(0) = 1;
			//Df
			x_lo(1) = 0;
			x_hi(1) = 1;
			//Ds
			x_lo(2) = 0;
			x_hi(2) = 1;
			//C (in normalized intensity scale)
			x_lo(3) = 0;
			x_hi(3) = 10;
			switch(fitting_method) {
				case OP_BFGS:
					info.desc = mfilename + " dlib constrained BFGS technique MI " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();
					cout << info.desc << endl;
					switch(no_bvalues) {
						case 2:
							FitBiexpConst_4params_2_bvalues_MultipleInitializations(indexes, SIdata, pinit, Biexp_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 3:
							FitBiexpConst_4params_3_bvalues_MultipleInitializations(indexes, SIdata, pinit, Biexp_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 4:
							FitBiexpConst_4params_4_bvalues_MultipleInitializations(indexes, SIdata, pinit, Biexp_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 5:
							FitBiexpConst_4params_5_bvalues_MultipleInitializations(indexes, SIdata, pinit, Biexp_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 6:
							FitBiexpConst_4params_6_bvalues_MultipleInitializations(indexes, SIdata, pinit, Biexp_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 7:
							FitBiexpConst_4params_7_bvalues_MultipleInitializations(indexes, SIdata, pinit, Biexp_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 8:
							FitBiexpConst_4params_8_bvalues_MultipleInitializations(indexes, SIdata, pinit, Biexp_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 9:
							FitBiexpConst_4params_9_bvalues_MultipleInitializations(indexes, SIdata, pinit, Biexp_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 10:
							FitBiexpConst_4params_10_bvalues_MultipleInitializations(indexes, SIdata, pinit, Biexp_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 11:
							FitBiexpConst_4params_11_bvalues_MultipleInitializations(indexes, SIdata, pinit, Biexp_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 12:
							FitBiexpConst_4params_12_bvalues_MultipleInitializations(indexes, SIdata, pinit, Biexp_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 13:
							FitBiexpConst_4params_13_bvalues_MultipleInitializations(indexes, SIdata, pinit, Biexp_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 14:
							FitBiexpConst_4params_14_bvalues_MultipleInitializations(indexes, SIdata, pinit, Biexp_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 15:
							FitBiexpConst_4params_15_bvalues_MultipleInitializations(indexes, SIdata, pinit, Biexp_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						default:
                            THROW_FATAL("Unsupported number of b-values")
					}
					break;
				case OP_BFGS_NORMALIZED:
					info.desc = mfilename + " dlib constrained BFGS technique MI, normalized SI data " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();
					cout << info.desc << endl;
					switch(no_bvalues) {
						case 12:
							FitBiexpConstN_3params_12_bvalues_MultipleInitializations(indexes, SIdata, pinit, BiexpN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 14:
							FitBiexpConstN_3params_14_bvalues_MultipleInitializations(indexes, SIdata, pinit, BiexpN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						case 15:
							FitBiexpConstN_3params_15_bvalues_MultipleInitializations(indexes, SIdata, pinit, BiexpN_residual_MSE, p_vals, info, x_lo, x_hi);
							break;
						default:
                            THROW_FATAL("Unsupported number of b-values")
					}
					break;
				case OP_LM:
					double (*residual_fun)(const std::pair<double, double>&, const matrix<double,4,1>&);
					if(weighted) {
						info.desc = mfilename + " dlib Traditional Levenberg-Marquardt technique MI, weighted " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();
						residual_fun = Biexp_residual_weights;
					} else {
						info.desc = mfilename + " dlib Traditional Levenberg-Marquardt technique MI " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();			
						residual_fun = Biexp_residual;
					}
					cout << info.desc << endl;
					switch(no_bvalues) {
						case 12:
							FitBiexp_4params_12_bvalues_MultipleInitializations(indexes, SIdata, pinit, residual_fun, p_vals, info, C_values);
							break;
						case 14:
							FitBiexp_4params_14_bvalues_MultipleInitializations(indexes, SIdata, pinit, residual_fun, p_vals, info, C_values);
							break;
						case 15:
							FitBiexp_4params_15_bvalues_MultipleInitializations(indexes, SIdata, pinit, residual_fun, p_vals, info, C_values);
							break;
						default:
                            THROW_FATAL("Unsupported number of b-values")
					}
					break;
				//simulations
				case OP_SIMULATION:
					resolveMultipleInitializations(limits, pinit, info, 0);
					//go through parameter set, writing file for each set
					unsigned int p_i;
					for(p_i = 0; p_i < (int)pinit.size(); p_i++) {
						std::stringstream output_suffix;
						std::vector<std::vector<double> > limits_instance;
						for(int param_i = 0; param_i < 4; param_i++) {
							output_suffix << "_" << info.parameterNames[param_i] <<  pinit[p_i][param_i];
							std::vector<double> lim(3);
							lim[0] = pinit[p_i][param_i];
							lim[1] = pinit[p_i][param_i];
							lim[2] = pinit[p_i][param_i];
							limits_instance.push_back(lim);
						}
						std::vector<std::vector<double> > pinit_instance;
						std::vector<std::vector<double> > simulation_pinit_instance;
						info.desc = mfilename + " creation of simulated curves " + resolveMultipleInitializations(limits_instance, pinit_instance, info, 0) + 
							" " + resolveMultipleInitializations(simulation_limits, simulation_pinit_instance, simulation_info, 0) + os.str();
						cout << info.desc << endl;

						output_suffix << "_SD" << simulation_pinit_instance[0][0];
						output_suffix << "_iterations" << simulation_pinit_instance[0][1];
						std::vector<std::vector<double> > SIdata_instance;
						Simulate_4params(indexes, SIdata_instance, pinit_instance, simulation_pinit_instance, Biexp_model, info);
						/*
						switch(no_bvalues) {
							case 12:
								Simulate_4params(indexes, SIdata_instance, pinit_instance, simulation_pinit_instance, Biexp_model, info);
								break;
							case 14:
								Simulate_4params(indexes, SIdata_instance, pinit_instance, simulation_pinit_instance, Biexp_model, info);
								break;
							case 15:
								Simulate_4params(indexes, SIdata_instance, pinit_instance, simulation_pinit_instance, Biexp_model, info);
								break;
							case 16:
								Simulate_4params(indexes, SIdata_instance, pinit_instance, simulation_pinit_instance, Biexp_model, info);
								break;
							default:
								THROW_FATAL("Unsupported number of b-values")
						}
						*/
						//write results
						if(fitting_method == OP_SIMULATION) {
							output_filename = filename + "_" + info.modelName + "_" + output_suffix.str() + "_simulation_ASCII.txt";
							cout << "writing :" << output_filename << endl;		
							writeSIfile(output_filename, info, SIdata_instance);
						}
					}
					break;
				default:
					THROW_FATAL("Unsupported fitting method")
			}
			#pragma endregion
		} else if(info.modelName.compare(MARZISEGMENTEDBIEXP) == 0) {
			#pragma region MARZISEGMENTEDBIEXP
			matrix<double,2,1> x_lo1;
			matrix<double,2,1> x_hi1;
			//Ds
			x_lo1(0) = 0;
			x_hi1(0) = 1;
			//C (in normalized intensity scale)
			x_lo1(1) = 0;
			x_hi1(1) = 10;
			matrix<double,3,1> x_lo2;
			matrix<double,3,1> x_hi2;
			//f
			x_lo2(0) = 0;
			x_hi2(0) = 1;
			//Df (Dp)
			x_lo2(1) = 0;
			x_hi2(1) = 1;
			//C (in normalized intensity scale)
			x_lo2(2) = 0;
			x_hi2(2) = 10;
			switch(fitting_method) {
				case OP_LM:
					double (*residual_fun_1st)(const std::pair<double, double>&, const matrix<double,2,1>&);
					double (*residual_fun_2nd)(const std::pair<double, double>&, const matrix<double,3,1>&);
					if(weighted) {
						THROW_FATAL("Weighted fit is not supported")
					} else {
                        VERBOSE2("resolving multiple initializations")
						info.desc = mfilename + " dlib Traditional Levenberg-Marquardt technique " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();
						residual_fun_1st = MonoExponential_residual;
						residual_fun_2nd = Biexp_residual_fixedDs;
                        VERBOSE2("resolving multiple initializations done")
					}
					cout << info.desc << endl;
					switch(No_tail_values) {
						case 3:
							switch(no_bvalues) {
								case 12:
									FitMarziSegmentedBiexp_params_12_bvalues_3tail_MultipleInitializations(indexes, SIdata, pinit, residual_fun_1st, residual_fun_2nd, p_vals, info, C_values);
									break;
								case 14:
									FitMarziSegmentedBiexp_params_14_bvalues_3tail_MultipleInitializations(indexes, SIdata, pinit, residual_fun_1st, residual_fun_2nd, p_vals, info, C_values);
									break;
								case 15:
									FitMarziSegmentedBiexp_params_15_bvalues_3tail_MultipleInitializations(indexes, SIdata, pinit, residual_fun_1st, residual_fun_2nd, p_vals, info, C_values);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						case 4:
							switch(no_bvalues) {
								case 12:
									FitMarziSegmentedBiexp_params_12_bvalues_4tail_MultipleInitializations(indexes, SIdata, pinit, MonoExponential_residual, Biexp_residual_fixedDs, p_vals, info, C_values);
									break;
								case 14:
									FitMarziSegmentedBiexp_params_14_bvalues_4tail_MultipleInitializations(indexes, SIdata, pinit, MonoExponential_residual, Biexp_residual_fixedDs, p_vals, info, C_values);
									break;
								case 15:
									FitMarziSegmentedBiexp_params_15_bvalues_4tail_MultipleInitializations(indexes, SIdata, pinit, MonoExponential_residual, Biexp_residual_fixedDs, p_vals, info, C_values);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						default:
							THROW_FATAL("Unsupported number of tail values")
						break;
					}
					break;
				case OP_BFGS:
					double (*residual_MSE_fun_1st)(const matrix<std::pair<double, double>, 0, 1>&, const matrix<double,0,1>&);
					double (*residual_MSE_fun_2nd)(const matrix<std::pair<double, double>, 0, 1>&, const matrix<double,0,1>&);
					if(weighted) {
						THROW_FATAL("Weighted fit is not supported")
					} else {
						limits_1st.push_back(limits[2]);
						limits_2nd.push_back(limits[0]);
						limits_2nd.push_back(limits[1]);
						info.desc = mfilename + " dlib constrained BFGS technique MI " + resolveMultipleInitializations(limits_1st, pinit, info, 2) + " " + resolveMultipleInitializations(limits_2nd, pinit2, info, 0) + os.str();
						residual_MSE_fun_1st = MonoExponential_residual_MSE;
						residual_MSE_fun_2nd = Biexp_residual_MSE_fixedDs;
					}
					cout << info.desc << endl;
					switch(No_tail_values) {
						case 3:
							switch(no_bvalues) {
								case 12:
									FitMarziSegmentedBiexpConst_params_12_bvalues_3tail_MultipleInitializations(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 14:
									FitMarziSegmentedBiexpConst_params_14_bvalues_3tail_MultipleInitializations(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 15:
									FitMarziSegmentedBiexpConst_params_15_bvalues_3tail_MultipleInitializations(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						case 4:
							switch(no_bvalues) {
								case 12:
									FitMarziSegmentedBiexpConst_params_12_bvalues_4tail_MultipleInitializations(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 14:
									FitMarziSegmentedBiexpConst_params_14_bvalues_4tail_MultipleInitializations(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 15:
									FitMarziSegmentedBiexpConst_params_15_bvalues_4tail_MultipleInitializations(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						case 7:
							switch(no_bvalues) {
								case 16:
									FitMarziSegmentedBiexpConst_params_16_bvalues_7tail_MultipleInitializations(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						case 8:
							switch(no_bvalues) {
								case 16:
									FitMarziSegmentedBiexpConst_params_16_bvalues_8tail_MultipleInitializations(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						case 9:
							switch(no_bvalues) {
								case 16:
									FitMarziSegmentedBiexpConst_params_16_bvalues_9tail_MultipleInitializations(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						default:
							THROW_FATAL("Unsupported number of tail values")
					}
					break;
				default:
					THROW_FATAL("Unsupported fitting method")
			}
			#pragma endregion
		} else if(info.modelName.compare(CHOSEGMENTEDBIEXP) == 0) {
			#pragma region CHOSEGMENTEDBIEXP
			matrix<double,2,1> x_lo1;
			matrix<double,2,1> x_hi1;
			//Ds
			x_lo1(0) = 0;
			x_hi1(0) = 1;
			//C (in normalized intensity scale)
			x_lo1(1) = 0;
			x_hi1(1) = 10;
			matrix<double,2,1> x_lo2;
			matrix<double,2,1> x_hi2;
			//Df (Dp)
			x_lo2(0) = 0;
			x_hi2(0) = 1;
			//C (in normalized intensity scale)
			x_lo2(1) = 0;
			x_hi2(1) = 10;
			switch(fitting_method) {
				case OP_BFGS:
					double (*residual_MSE_fun_1st)(const matrix<std::pair<double, double>, 0, 1>&, const matrix<double,0,1>&);
					double (*residual_MSE_fun_2nd)(const matrix<std::pair<double, double>, 0, 1>&, const matrix<double,0,1>&);
					if(weighted) {
						THROW_FATAL("Weighted fit is not supported")
					} else {
						limits_1st.push_back(limits[1]);
						limits_2nd.push_back(limits[0]);
						info.desc = mfilename + " dlib constrained BFGS technique MI " + resolveMultipleInitializations(limits_1st, pinit, info, 1) + " " + resolveMultipleInitializations(limits_2nd, pinit2, info, 0) + os.str();
						residual_MSE_fun_1st = MonoExponential_residual_MSE;
						residual_MSE_fun_2nd = Biexp_residual_MSE_fixedDsf;
					}
					cout << info.desc << endl;
					switch(No_tail_values) {
						case 3:
							switch(no_bvalues) {
								case 12:
									FitChoSegmentedBiexpConst_params_12_bvalues_3tail_MultipleInitializations(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 14:
									FitChoSegmentedBiexpConst_params_14_bvalues_3tail_MultipleInitializations(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 15:
									FitChoSegmentedBiexpConst_params_15_bvalues_3tail_MultipleInitializations(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						case 4:
							switch(no_bvalues) {
								case 12:
									FitChoSegmentedBiexpConst_params_12_bvalues_4tail_MultipleInitializations(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 14:
									FitChoSegmentedBiexpConst_params_14_bvalues_4tail_MultipleInitializations(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 15:
									FitChoSegmentedBiexpConst_params_15_bvalues_4tail_MultipleInitializations(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
						case 7:
							switch(no_bvalues) {
								case 16:
									FitChoSegmentedBiexpConst_params_16_bvalues_7tail_MultipleInitializations(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
						case 8:
							switch(no_bvalues) {
								case 16:
									FitChoSegmentedBiexpConst_params_16_bvalues_8tail_MultipleInitializations(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
						case 9:
							switch(no_bvalues) {
								case 16:
									FitChoSegmentedBiexpConst_params_16_bvalues_9tail_MultipleInitializations(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						default:
							THROW_FATAL("Unsupported number of tail values")
					}
					break;
				default:
					THROW_FATAL("Unsupported fitting method")
			}
			//Insert C as 3rd parameter name, because it was not used in 
			//fittings, but still reported for sake of completeness
			info.parameterNames.insert(info.parameterNames.begin(), "f");
			#pragma endregion
		} else if(info.modelName.compare(LOGINITIVIM) == 0) {
                        #pragma region LOGINITIVIM
			matrix<double,2,1> x_lo1;
			matrix<double,2,1> x_hi1;
			//Ds
			x_lo1(0) = 0;
			x_hi1(0) = 1;
			//C (in normalized intensity scale)
			x_lo1(1) = 0;
			x_hi1(1) = 10;
			matrix<double,3,1> x_lo2;
			matrix<double,3,1> x_hi2;
			//f
			x_lo2(0) = 0;
			x_hi2(0) = 1;
			//Df (Dp)
			x_lo2(1) = 0;
			x_hi2(1) = 1;
			//C (in normalized intensity scale)
			x_lo2(2) = 0;
			x_hi2(2) = 10;
			switch(fitting_method) {
				case OP_BFGS_NORMALIZED:
					double (*residual_MSE_fun_1st)(const matrix<std::pair<double, double>, 0, 1>&, const matrix<double,0,1>&);
					double (*residual_MSE_fun_2nd)(const matrix<std::pair<double, double>, 0, 1>&, const matrix<double,0,1>&);
					if(weighted) {
						THROW_FATAL("Weighted fit is not supported")
					} else {
                        VERBOSE2("number of limits " << limits.size())
						limits_1st.push_back(limits[2]);
						limits_2nd.push_back(limits[0]);
						limits_2nd.push_back(limits[1]);
                        VERBOSE2("resolving multiple initializations")
						info.desc = mfilename + " dlib constrained BFGS technique MI " + resolveMultipleInitializations(limits_1st, pinit, info, 2) + " " + resolveMultipleInitializations(limits_2nd, pinit2, info, 0) + os.str();
                        VERBOSE2("resolving multiple initializations done")
						residual_MSE_fun_1st = MonoExponential_residual_MSE;
						residual_MSE_fun_2nd = Biexp_residual_MSE_fixedDs;
					}
					cout << info.desc << endl;
					switch(No_tail_values) {
						case 2:
							switch(no_bvalues) {
								case 3:
									FitLogInitBiexp_3_bvalues_2tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 4:
									FitLogInitBiexp_4_bvalues_2tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 5:
									FitLogInitBiexp_5_bvalues_2tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 6:
									FitLogInitBiexp_6_bvalues_2tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 7:
									FitLogInitBiexp_7_bvalues_2tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 8:
									FitLogInitBiexp_8_bvalues_2tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 9:
									FitLogInitBiexp_9_bvalues_2tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 10:
									FitLogInitBiexp_10_bvalues_2tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 11:
									FitLogInitBiexp_11_bvalues_2tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 12:
									FitLogInitBiexp_12_bvalues_2tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 13:
									FitLogInitBiexp_13_bvalues_2tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 14:
									FitLogInitBiexp_14_bvalues_2tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 15:
									FitLogInitBiexp_15_bvalues_2tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 16:
									FitLogInitBiexp_16_bvalues_2tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						case 3:
							switch(no_bvalues) {
								case 4:
									FitLogInitBiexp_4_bvalues_3tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 5:
									FitLogInitBiexp_5_bvalues_3tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 6:
									FitLogInitBiexp_6_bvalues_3tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 7:
									FitLogInitBiexp_7_bvalues_3tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 8:
									FitLogInitBiexp_8_bvalues_3tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 9:
									FitLogInitBiexp_9_bvalues_3tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 10:
									FitLogInitBiexp_10_bvalues_3tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 11:
									FitLogInitBiexp_11_bvalues_3tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 12:
									FitLogInitBiexp_12_bvalues_3tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 13:
									FitLogInitBiexp_13_bvalues_3tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 14:
									FitLogInitBiexp_14_bvalues_3tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 15:
									FitLogInitBiexp_15_bvalues_3tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 16:
									FitLogInitBiexp_16_bvalues_3tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						case 4:
							switch(no_bvalues) {
								case 5:
									FitLogInitBiexp_5_bvalues_4tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 6:
									FitLogInitBiexp_6_bvalues_4tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 7:
									FitLogInitBiexp_7_bvalues_4tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 8:
									FitLogInitBiexp_8_bvalues_4tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 9:
									FitLogInitBiexp_9_bvalues_4tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 10:
									FitLogInitBiexp_10_bvalues_4tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 11:
									FitLogInitBiexp_11_bvalues_4tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 12:
									FitLogInitBiexp_12_bvalues_4tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 13:
									FitLogInitBiexp_13_bvalues_4tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 14:
									FitLogInitBiexp_14_bvalues_4tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 15:
									FitLogInitBiexp_15_bvalues_4tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 16:
									FitLogInitBiexp_16_bvalues_4tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						case 5:
							switch(no_bvalues) {
								case 6:
									FitLogInitBiexp_6_bvalues_5tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 7:
									FitLogInitBiexp_7_bvalues_5tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 8:
									FitLogInitBiexp_8_bvalues_5tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 9:
									FitLogInitBiexp_9_bvalues_5tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 10:
									FitLogInitBiexp_10_bvalues_5tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 11:
									FitLogInitBiexp_11_bvalues_5tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 12:
									FitLogInitBiexp_12_bvalues_5tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 13:
									FitLogInitBiexp_13_bvalues_5tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 14:
									FitLogInitBiexp_14_bvalues_5tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 15:
									FitLogInitBiexp_15_bvalues_5tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 16:
									FitLogInitBiexp_16_bvalues_5tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						case 6:
							switch(no_bvalues) {
								case 7:
									FitLogInitBiexp_7_bvalues_6tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 8:
									FitLogInitBiexp_8_bvalues_6tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 9:
									FitLogInitBiexp_9_bvalues_6tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 10:
									FitLogInitBiexp_10_bvalues_6tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 11:
									FitLogInitBiexp_11_bvalues_6tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 12:
									FitLogInitBiexp_12_bvalues_6tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 13:
									FitLogInitBiexp_13_bvalues_6tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 14:
									FitLogInitBiexp_14_bvalues_6tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 15:
									FitLogInitBiexp_15_bvalues_6tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 16:
									FitLogInitBiexp_16_bvalues_6tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						case 7:
							switch(no_bvalues) {
								case 8:
									FitLogInitBiexp_8_bvalues_7tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 9:
									FitLogInitBiexp_9_bvalues_7tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 10:
									FitLogInitBiexp_10_bvalues_7tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 11:
									FitLogInitBiexp_11_bvalues_7tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 12:
									FitLogInitBiexp_12_bvalues_7tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 13:
									FitLogInitBiexp_13_bvalues_7tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 14:
									FitLogInitBiexp_14_bvalues_7tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 15:
									FitLogInitBiexp_15_bvalues_7tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 16:
									FitLogInitBiexp_16_bvalues_7tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						case 8:
							switch(no_bvalues) {
								case 9:
									FitLogInitBiexp_9_bvalues_8tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 10:
									FitLogInitBiexp_10_bvalues_8tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 11:
									FitLogInitBiexp_11_bvalues_8tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 12:
									FitLogInitBiexp_12_bvalues_8tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 13:
									FitLogInitBiexp_13_bvalues_8tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 14:
									FitLogInitBiexp_14_bvalues_8tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 15:
									FitLogInitBiexp_15_bvalues_8tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 16:
									FitLogInitBiexp_16_bvalues_8tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						case 9:
							switch(no_bvalues) {
								case 10:
									FitLogInitBiexp_10_bvalues_9tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 11:
									FitLogInitBiexp_11_bvalues_9tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 12:
									FitLogInitBiexp_12_bvalues_9tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 13:
									FitLogInitBiexp_13_bvalues_9tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 14:
									FitLogInitBiexp_14_bvalues_9tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 15:
									FitLogInitBiexp_15_bvalues_9tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 16:
									FitLogInitBiexp_16_bvalues_9tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						case 10:
							switch(no_bvalues) {
								case 11:
									FitLogInitBiexp_11_bvalues_10tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 12:
									FitLogInitBiexp_12_bvalues_10tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 13:
									FitLogInitBiexp_13_bvalues_10tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 14:
									FitLogInitBiexp_14_bvalues_10tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 15:
									FitLogInitBiexp_15_bvalues_10tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 16:
									FitLogInitBiexp_16_bvalues_10tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						case 11:
							switch(no_bvalues) {
								case 12:
									FitLogInitBiexp_12_bvalues_11tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 13:
									FitLogInitBiexp_13_bvalues_11tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 14:
									FitLogInitBiexp_14_bvalues_11tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 15:
									FitLogInitBiexp_15_bvalues_11tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 16:
									FitLogInitBiexp_16_bvalues_11tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						case 12:
							switch(no_bvalues) {
								case 13:
									FitLogInitBiexp_13_bvalues_12tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 14:
									FitLogInitBiexp_14_bvalues_12tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 15:
									FitLogInitBiexp_15_bvalues_12tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 16:
									FitLogInitBiexp_16_bvalues_12tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						case 13:
							switch(no_bvalues) {
								case 14:
									FitLogInitBiexp_14_bvalues_13tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 15:
									FitLogInitBiexp_15_bvalues_13tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 16:
									FitLogInitBiexp_16_bvalues_13tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						case 14:
							switch(no_bvalues) {
								case 15:
									FitLogInitBiexp_15_bvalues_14tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								case 16:
									FitLogInitBiexp_16_bvalues_14tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						case 15:
							switch(no_bvalues) {
								case 16:
									FitLogInitBiexp_16_bvalues_15tail_MultipleInitializations_BFGS(indexes, SIdata, pinit, pinit2, residual_MSE_fun_1st, residual_MSE_fun_2nd, p_vals, info, C_values, x_lo1, x_hi1, x_lo2, x_hi2);
									break;
								default:
									THROW_FATAL("Unsupported number of b-values")
							}
							break;
						default:
							THROW_FATAL("Unsupported number of tail values")
					}
					break;
				default:
					THROW_FATAL("Unsupported fitting method")
			}
            #pragma endregion
		// Non-Gaussian IntraVoxel Incoherent Motion model (NG-IVIM)
		} else if(info.modelName.compare(NGIVIM) == 0) {
			#pragma region NGIVIM
			matrix<double,5,1> x_lo;
			matrix<double,5,1> x_hi;
			//f weighting for pseudo-diffusion component
			x_lo(0) = 0;
			x_hi(0) = 1;
			//Dp Pseudo-Diffusion in IVIM
			x_lo(1) = 0;
			x_hi(1) = 1;
			//D Diffusion
			x_lo(2) = 0;
			x_hi(2) = 1;
			//K deviation from Gaussian diffusion
			x_lo(3) = 0;
			x_hi(3) = 10;
			//C intensity
			x_lo(4) = 0;
			x_hi(4) = 10;
			info.desc = mfilename + " dlib constrained BFGS technique MI " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();
			cout << info.desc << endl;
			switch(no_bvalues) {
				case 12:
					FitConst_5params_12_bvalues_MultipleInitializations(indexes, SIdata, pinit, NGIVIM_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
					break;
				case 14:
					FitConst_5params_14_bvalues_MultipleInitializations(indexes, SIdata, pinit, NGIVIM_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
					break;
				case 15:
					FitConst_5params_15_bvalues_MultipleInitializations(indexes, SIdata, pinit, NGIVIM_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
					break;
				default:
					THROW_FATAL("Unsupported number of b-values")
			}
			#pragma endregion
		// Non-Gaussian IntraVoxel Incoherent Motion model (NG-IVIM), linearized simplification
		} else if(info.modelName.compare(LNGIVIM) == 0) {
			#pragma region LNGIVIM
			matrix<double,2,1> x_lo;
			matrix<double,2,1> x_hi;
			//f weighting for pseudo-diffusion component
			x_lo(0) = 0;
			x_hi(0) = 1;
			//D Diffusion
			x_lo(1) = 0;
			x_hi(1) = 1;
			info.desc = mfilename + " dlib constrained BFGS technique MI " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();
			cout << info.desc << endl;
			switch(no_bvalues) {
				case 12:
					Fit_2params_12_bvalues_3tail_MultipleInitializations(indexes, SIdata, pinit, LNGIVIM_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
					break;
				case 14:
					Fit_2params_14_bvalues_3tail_MultipleInitializations(indexes, SIdata, pinit, LNGIVIM_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
					break;
				case 15:
					Fit_2params_15_bvalues_3tail_MultipleInitializations(indexes, SIdata, pinit, LNGIVIM_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
					break;
				default:
					THROW_FATAL("Unsupported number of b-values")
			}
			#pragma endregion
		// Non-Gaussian IntraVoxel Incoherent Motion model (NG-IVIM), linearized simplification
		// [1] J. Pekar, C.T. Moonen, P.C. van Zijl, 
		// On the precision of diffusion/perfusion imaging by gradient sensitization., 
		// Magn. Reson. Med. 23 (1992) 122–9.
		} else if(info.modelName.compare(ASYMPTOTICBIEXP) == 0) {
			#pragma region ASYMPTOTICBIEXP
			matrix<double,2,1> x_lo;
			matrix<double,2,1> x_hi;
			//D Diffusion
			x_lo(0) = 0;
			x_hi(0) = 1;
			//C (in normalized intensity scale)
			x_lo(1) = 0;
			x_hi(1) = 10;
			info.desc = mfilename + " dlib constrained BFGS technique MI " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();
			cout << info.desc << endl;
			switch(No_tail_values) {
				case 3:
					switch(no_bvalues) {
						case 12:
							AsympFit_2params_12_bvalues_3tail_MultipleInitializations(indexes, SIdata, pinit, MonoExponential_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 14:
							AsympFit_2params_14_bvalues_3tail_MultipleInitializations(indexes, SIdata, pinit, MonoExponential_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						case 15:
							AsympFit_2params_15_bvalues_3tail_MultipleInitializations(indexes, SIdata, pinit, MonoExponential_residual_MSE, p_vals, info, C_values, x_lo, x_hi);
							break;
						default:
							THROW_FATAL("Unsupported number of b-values")
					}
					break;
				default:
					THROW_FATAL("Unsupported number of tail values")
			}
			//Insert C as 3rd parameter name, because it was not used in 
			//fittings, but still reported for sake of completeness
			info.parameterNames.insert(info.parameterNames.begin()+2, "C");
			#pragma endregion
		// Non-Gaussian IntraVoxel Incoherent Motion model (NG-IVIM), log-linearized simplification
		// [1] Y. Pang, B. Turkbey, M. Bernardo, J. Kruecker, S. Kadoury, M.J. Merino, et al., 
		// Intravoxel incoherent motion MR imaging for prostate cancer: an evaluation of perfusion 
		// fraction and diffusion coefficient derived from different b-value combinations., Magn. Reson. Med. 69 (2013) 553–62. 
		// doi:10.1002/mrm.24277. 
		} else if(info.modelName.compare(LOGLINASYMPTOTICBIEXP) == 0) {
			#pragma region LOGLINASYMPTOTICBIEXP
			matrix<double,2,1> x_lo;
			matrix<double,2,1> x_hi;
			//D Diffusion
			x_lo(0) = 0;
			x_hi(0) = 1;
			//C (in normalized intensity scale)
			x_lo(1) = 0;
			x_hi(1) = 10;
			info.desc = mfilename + " dlib constrained BFGS technique (linear regression fit with no initializations) " + os.str();
			cout << info.desc << endl;
			switch(No_tail_values) {
				case 2:
					switch(no_bvalues) {
						case 12:
							LogLinAsympFit_2params_12_bvalues_2tail_MultipleInitializations(indexes, SIdata, MonoExponential_residual_MSE, p_vals, info);
							break;
						case 14:
							LogLinAsympFit_2params_14_bvalues_2tail_MultipleInitializations(indexes, SIdata, MonoExponential_residual_MSE, p_vals, info);
							break;
						case 15:
							LogLinAsympFit_2params_15_bvalues_2tail_MultipleInitializations(indexes, SIdata, MonoExponential_residual_MSE, p_vals, info);
							break;
						default:
							THROW_FATAL("Unsupported number of b-values")
					}
					break;
				case 3:
					switch(no_bvalues) {
						case 12:
							LogLinAsympFit_2params_12_bvalues_3tail_MultipleInitializations(indexes, SIdata, MonoExponential_residual_MSE, p_vals, info);
							break;
						case 14:
							LogLinAsympFit_2params_14_bvalues_3tail_MultipleInitializations(indexes, SIdata, MonoExponential_residual_MSE, p_vals, info);
							break;
						case 15:
							LogLinAsympFit_2params_15_bvalues_3tail_MultipleInitializations(indexes, SIdata, MonoExponential_residual_MSE, p_vals, info);
							break;
						default:
							THROW_FATAL("Unsupported number of b-values")
					}
					break;
				default:
					THROW_FATAL("Unsupported number of tail values")
			}
			#pragma endregion
		// Non-Gaussian IntraVoxel Incoherent Motion model (SNG-IVIM), linearized simplification with simple linear regression
		} else if(info.modelName.compare(SLNGIVIM) == 0) {
			#pragma region SLNGIVIM
			matrix<double,2,1> x_lo;
			matrix<double,2,1> x_hi;
			//f weighting for pseudo-diffusion component
			x_lo(0) = 0;
			x_hi(0) = 1;
			//D Diffusion
			x_lo(1) = 0;
			x_hi(1) = 1;
			info.desc = mfilename + " dlib simple linear regression line " + resolveMultipleInitializations(limits, pinit, info, 0) + os.str();
			cout << info.desc << endl;
			switch(no_bvalues) {
				case 12:
					SFit_2params_12_bvalues_3tail_MultipleInitializations(indexes, SIdata, pinit, LNGIVIM_residual_MSE, p_vals, info, x_lo, x_hi);
					break;
				case 14:
					SFit_2params_14_bvalues_3tail_MultipleInitializations(indexes, SIdata, pinit, LNGIVIM_residual_MSE, p_vals, info, x_lo, x_hi);
					break;
				case 15:
					SFit_2params_15_bvalues_3tail_MultipleInitializations(indexes, SIdata, pinit, LNGIVIM_residual_MSE, p_vals, info, x_lo, x_hi);
					break;
				default:
					THROW_FATAL("Unsupported number of b-values")
			}
			#pragma endregion
		} else {
			THROW_FATAL("Model was not recognized")
		}
		//write results for other than simulation
		if(fitting_method != OP_SIMULATION) {
			output_filename = output_filename + "_" + info.modelName + (weighted ? "_Weighted" : "") + (fitting_method == OP_BFGS_NORMALIZED ? "N" : "") + "_results.txt";
			cout << "writing :" << output_filename << endl;		
			writefile(output_filename, p_vals, info);
		}
    }
	catch (std::exception& e)
    {
        cout << "EXCEPTION:" << e.what() << endl;
		return 1;
    }
    catch (std::string s)
    {
        cout << "ERROR MESSAGE:" << s << endl;
		return 1;
    }
	catch (...) { 
		cout << "unknown exception";
	}
	return 0;
}
// ----------------------------------------------------------------------------------------
#pragma region MonoExponential
double MonoN_model(const double& xvalue, const vectorParam1& params)
{
    double ADCm = params(0);
	double b = xvalue;
	double r = std::exp(-b*ADCm);
    //cout << "r=" << r << "ADCm=" << ADCm << "b=" << b << std::setprecision(OUTPUT_PRECISION) << std::fixed << endl; 
	return r;
}
// This function is the "residual" for a least squares problem.   It takes an input/output
// pair and compares it to the output of our model and returns the amount of error.  The idea
// is to find the set of parameters which makes the residual small on all the data pairs.
double MonoN_residual_MSE (
    const std::pair<double, double>& data,
    const matrix<double,1,1>& params
)
{
	return MonoN_model(data.first, params) - data.second;
}

double MonoExponential_model(const double& xvalue, const vectorParam2& params)
{
    double ADCm = params(0);
    double C = params(1);
	double b = xvalue;

	double r = C*std::exp(-b*ADCm);
	return r;
}
// This function is the "residual" for a least squares problem.   It takes an input/output
// pair and compares it to the output of our model and returns the amount of error.  The idea
// is to find the set of parameters which makes the residual small on all the data pairs.
double MonoExponential_residual (
    const std::pair<double, double>& data,
    const matrix<double,2,1>& params
)
{
	return MonoExponential_model(data.first, params) - data.second;
}
double MonoExponential_residual_weights(
    const std::pair<double, double>& data,
    const matrix<double,2,1>& params
)
{
	double r = (MonoExponential_model(bvalues[(int)data.first], params) - data.second)*weights[(int)data.first];
	return r;
}
double MonoExponential_residual_MSE (const matrix<std::pair<double, double>, 0, 1>& data, const matrix<double, 0, 1>& params)
{
	double MSE = 0.0;
	double err;
	std::pair<double, double> pair;
	if(VERBOSE > 1) {
		for(long i = 0; i < params.size(); i++) {
			cout << " " << params(i);
		}
		cout << endl;
	}
	for(long i = 0; i < data.size(); i++) {
		pair = data(i);
		err = MonoExponential_residual(pair, params);
		if(VERBOSE > 1)
			cout << " " << err;
		MSE = MSE + err*err;
	}
	if(VERBOSE > 1)
		cout << endl;
	if((int)data.size() == 0) return -1;
	MSE = MSE/((double)data.size());
#if VERBOSE==2
	cout << "mono " << params(0) << " " << params(1) << " > " << MSE << endl;
#endif
	return MSE;
}
#pragma endregion

// ----------------------------------------------------------------------------------------
#pragma region MonoT1rho

double MonoT1rho_model(const double& xvalue, const vectorParam2& params)
{
    double rho = params(0);
    double C = params(1);
	double b = xvalue;
	if(rho == 0)
		rho = 1e-7;
	double r = C*std::exp(-b/rho);
	return r;
}
// This function is the "residual" for a least squares problem.   It takes an input/output
// pair and compares it to the output of our model and returns the amount of error.  The idea
// is to find the set of parameters which makes the residual small on all the data pairs.
double MonoT1rho_residual (
    const std::pair<double, double>& data,
    const vectorParam2& params
)
{
	return MonoT1rho_model(data.first, params) - data.second;
}
double MonoT1rho_residual_weights(
    const std::pair<double, double>& data,
    const vectorParam2& params
)
{
	double r = (MonoT1rho_model(bvalues[(int)data.first], params) - data.second)*weights[(int)data.first];
	return r;
}
double MonoT1rho_residual_MSE (const matrix<std::pair<double, double>, 0, 1>& data, const matrix<double,0,1>& params)
{
	double MSE = 0.0;
	double err;
	std::pair<double, double> pair;
	for(long i = 0; i < data.size(); i++) {
		pair = data(i);
		err = MonoT1rho_residual(pair, params);
		MSE = MSE + err*err;
	}
	if((int)data.size() == 0) return -1;
	MSE = MSE/((double)data.size());
#if VERBOSE==2
	cout << "monoT1rho " << params(0) << " " << params(1) << " > " << MSE << endl;
#endif
	return MSE;
}

#pragma endregion

// ----------------------------------------------------------------------------------------
#pragma region Log-Linearized IVIM

double LNGIVIM_model(const double& xvalue, const vectorParam2& params)
{
    double f = params(0);
	double D = params(1);
	double b = xvalue;

	double r;
	r = -f-b*D;
	return r;
}
// This function is the "residual" for a least squares problem.   It takes an input/output
// pair and compares it to the output of our model and returns the amount of error.  The idea
// is to find the set of parameters which makes the residual small on all the data pairs.
double LNGIVIM_residual (
    const std::pair<double, double>& data,
    const vectorParam2& params
)
{
    return LNGIVIM_model(data.first, params) - data.second;
}
//Mean Squared Error of data to the fitted curve
double LNGIVIM_residual_MSE (
    const matrix<std::pair<double, double>, 0, 1>& data,
    const matrix<double,0,1>& params
)
{
	double MSE = 0.0;
	double err;
	std::pair<double, double> pair;
	for(long i = 0; i < data.size(); i++) {
		pair = data(i);
		err = LNGIVIM_residual(pair, params);
		MSE = MSE + err*err;
	}
	if((int)data.size() == 0) return -1;
	MSE = MSE/((double)data.size());
	return MSE;
}
double LNGIVIM_residual_weights (
    const std::pair<double, double>& data,
    const vectorParam2& params
)
{
	return (LNGIVIM_model(bvalues[(int)data.first], params) - data.second)*weights[(int)data.first];
}

#pragma endregion

// ----------------------------------------------------------------------------------------
#pragma region Non-Gaussian IVIM, i. e. with kurtosity

double NGIVIM_model(const double& xvalue, const vectorParam5& params)
{
    double f = params(0);
    double Dp = params(1);
	double D = params(2);
	double K = params(3);
    double C = params(4);
	double b = xvalue;

	double alpha = b*D-b*b*D*D*K/6.0;
	double beta = b*Dp;

	double r;
	r = C*((1-f)*std::exp(-alpha)+f*std::exp(-alpha-beta));
	return r;
}
// This function is the "residual" for a least squares problem.   It takes an input/output
// pair and compares it to the output of our model and returns the amount of error.  The idea
// is to find the set of parameters which makes the residual small on all the data pairs.
double NGIVIM_residual (
    const std::pair<double, double>& data,
    const vectorParam5& params
)
{
    return NGIVIM_model(data.first, params) - data.second;
}
//Mean Squared Error of data to the fitted curve
double NGIVIM_residual_MSE (
    const matrix<std::pair<double, double>, 0, 1>& data,
    const matrix<double,0,1>& params
)
{
	double MSE = 0.0;
	double err;
	std::pair<double, double> pair;
	for(long i = 0; i < data.size(); i++) {
		pair = data(i);
		err = NGIVIM_residual(pair, params);
		MSE = MSE + err*err;
	}
	if((int)data.size() == 0) return -1;
	MSE = MSE/((double)data.size());
	return MSE;
}
double NGIVIM_residual_weights (
    const std::pair<double, double>& data,
    const vectorParam5& params
)
{
	return (NGIVIM_model(bvalues[(int)data.first], params) - data.second)*weights[(int)data.first];
}

#pragma endregion

// ----------------------------------------------------------------------------------------
#pragma region Kurt

double Kurt_model(const double& xvalue, const vectorParam3& params)
{
    double ADCk = params(0);
	double K = params(1);
    double C = params(2);
	double b = xvalue;
	double r;
	r = C*std::exp(-b*ADCk+1.0/6.0*b*b*ADCk*ADCk*K);
	if(!is_finite(r)) {
		if(r > 0)
			r = std::numeric_limits<double>::max();
		else
			r = std::numeric_limits<double>::min();
	}
	return r;
}
// This function is the "residual" for a least squares problem.   It takes an input/output
// pair and compares it to the output of our model and returns the amount of error.  The idea
// is to find the set of parameters which makes the residual small on all the data pairs.
double Kurt_residual (
    const std::pair<double, double>& data,
    const vectorParam3& params
)
{
    return Kurt_model(data.first, params) - data.second;
}
//Mean Squared Error of data to the fitted curve
double Kurt_residual_MSE (
    const matrix<std::pair<double, double>, 0, 1>& data,
    const matrix<double,0,1>& params
)
{
	double MSE = 0.0;
	double err;
	std::pair<double, double> pair;
	if(VERBOSE > 1) {
		for(long i = 0; i < params.size(); i++) {
			cout << " " << params(i);
		}
		cout << endl;
	}
	for(long i = 0; i < data.size(); i++) {
		pair = data(i);
		err = Kurt_residual(pair, params);
		if(VERBOSE > 1)
			cout << " " << err;
		MSE = MSE + err*err;
	}
	if(VERBOSE > 1)
		cout << endl;
	if((int)data.size() == 0) return -1;
	MSE = MSE/((double)data.size());
	return MSE;
}
double Kurt_residual_weights (
    const std::pair<double, double>& data,
    const vectorParam3& params
)
{
	return (Kurt_model(bvalues[(int)data.first], params) - data.second)*weights[(int)data.first];
}
//constrains for multiple initilization fits
double Kurt_const(const vectorParam3& params)
{
    double ADCk = params(0);
	double K = params(1);
    double C = params(2);
	double r = 0;
	if(K < 0) r = -K;
	return r;
}
#pragma endregion

// ----------------------------------------------------------------------------------------
#pragma region Stretched

double Stretched_model(const double& xvalue, const vectorParam3& params)
{
    double ADCs = params(0);
	double Alpha = params(1);
    double C = params(2);
	double b = xvalue;
	double ADCss;
	//explicit assertion that the base of exponential
	//is allways non-negative, because negative base are not defined for this model
	if(b*ADCs < 0)
		ADCss = 0;
	else
		ADCss = std::pow(b*ADCs, Alpha);
	double r = C*std::exp(-ADCss);
	return r;
}
// This function is the "residual" for a least squares problem.   It takes an input/output
// pair and compares it to the output of our model and returns the amount of error.  The idea
// is to find the set of parameters which makes the residual small on all the data pairs.
double Stretched_residual (
    const std::pair<double, double>& data,
    const vectorParam3& params
)
{
    return Stretched_model(data.first, params) - data.second;
}
double Stretched_residual_weights (
    const std::pair<double, double>& data,
    const vectorParam3& params
)
{
	return (Stretched_model(bvalues[(int)data.first], params) - data.second)*weights[(int)data.first];
}
//Mean Squared Error of data to the fitted curve
double Stretched_residual_MSE (
    const matrix<std::pair<double, double>, 0, 1>& data,
    const matrix<double,0,1>& params
)
{
	double MSE = 0.0;
	double err;
	std::pair<double, double> pair;
//	cout << params(0) << " " << params(1) << " " << params(2) << endl;
	for(long i = 0; i < data.size(); i++) {
		pair = data(i);
		err = Stretched_residual(pair, params);
	//	cout << err << " ";
		MSE = MSE + err*err;
	}
//	cout << endl;
	if((int)data.size() == 0) return -1;
	MSE = MSE/((double)data.size());
	return MSE;
}

#pragma endregion

// ----------------------------------------------------------------------------------------
#pragma region Biexp

double Biexp_model(const double& xvalue, const vectorParam4& params)
{
    double f = params(0);
    double Df = params(1);
	double Ds = params(2);
    double C = params(3);
	double b = xvalue;

	double r;
	r = C*(f*std::exp(-b*Df)+(1-f)*std::exp(-b*Ds));
	return r;
}
// This function is the "residual" for a least squares problem.   It takes an input/output
// pair and compares it to the output of our model and returns the amount of error.  The idea
// is to find the set of parameters which makes the residual small on all the data pairs.
double Biexp_residual (
    const std::pair<double, double>& data,
    const vectorParam4& params
)
{
    return Biexp_model(data.first, params) - data.second;
}
double Biexp_residual_weights (
    const std::pair<double, double>& data,
    const vectorParam4& params
)
{
	return (Biexp_model(bvalues[(int)data.first], params) - data.second)*weights[(int)data.first];
}
// This function is the "residual" for a least squares problem.   It takes an input/output
// pair and compares it to the output of our model and returns the amount of error.  The idea
// is to find the set of parameters which makes the residual small on all the data pairs.
double Biexp_residual_fixedDsf (
    const std::pair<double, double>& data,
    const vectorParam2& params
)
{
    double f = fixed_parameters[0];
	double Df = params(0);
	double Ds = fixed_parameters[1];
    double C = params(1);
	double b = data.first;
	double r = C*(f*std::exp(-b*Df)+(1-f)*std::exp(-b*Ds));

    return r - data.second;
}

// This function is the "residual" for a least squares problem.   It takes an input/output
// pair and compares it to the output of our model and returns the amount of error.  The idea
// is to find the set of parameters which makes the residual small on all the data pairs.
double Biexp_residual_fixedDs (const std::pair<double, double>& data, const vectorParam3& params)
{
    double f = params(0);
	double Df = params(1);
	double Ds = fixed_parameters[0];
    double C = params(2);
	double b = data.first;
	double r = C*(f*std::exp(-b*Df)+(1-f)*std::exp(-b*Ds));

    return r - data.second;
}
double Biexp_residual_MSE_fixedDs (const matrix<std::pair<double, double>, 0, 1>& data, const matrix<double,0,1>& params)
{
	double MSE = 0.0;
	double err;
	std::pair<double, double> pair;
	for(long i = 0; i < data.size(); i++) {
		pair = data(i);
		err = Biexp_residual_fixedDs(pair, params);
		MSE = MSE + err*err;
	}
	if((int)data.size() == 0) return -1;
	MSE = MSE/((double)data.size());
	return MSE;
}

double BiexpN_residual_fixedDs_mtx (const std::pair<double, double>& data, const vectorParam2& params)
{
    double f = params(0);
	double Df = params(1);
	double Ds = fixed_parameters[0];
	double b = data.first;
	double r = (f*std::exp(-b*Df)+(1-f)*std::exp(-b*Ds));

    return r - data.second;
}
double BiexpN_residual_MSE_fixedDs_mtx (const matrix<std::pair<double, double>, 0, 1>& data, const matrix<double,0,1>& params)
{
	double MSE = 0.0;
	double err;
	std::pair<double, double> pair;
	for(long i = 0; i < data.size(); i++) {
		pair = data(i);
		err = BiexpN_residual_fixedDs_mtx(pair, params);
		MSE = MSE + err*err;
	}
	if((int)data.size() == 0) return -1;
	MSE = MSE/((double)data.size());
	return MSE;
}

//Mean Squared Error of data to the fitted curve
double Biexp_residual_MSE (const matrix<std::pair<double, double>, 0, 1>& data, const matrix<double,0,1>& params)
{
	double MSE = 0.0;
	double err;
	std::pair<double, double> pair;
	for(long i = 0; i < data.size(); i++) {
		pair = data(i);
		err = Biexp_residual(pair, params);
		MSE = MSE + err*err;
	}
	if((int)data.size() == 0) return -1;
	MSE = MSE/((double)data.size());
	return MSE;
}
double BiexpN_residual_fixedDs (const double& xvalue, const vectorParam2& params)
{
    double f = params(0);
	double Df = params(1);
	double Ds = fixed_parameters[0];
	double b = xvalue;
	double r = (f*std::exp(-b*Df)+(1-f)*std::exp(-b*Ds));
    return r;
}
double BiexpN_residual_MSE_fixedDs (const std::pair<double, double>& data, const vectorParam2& params)
{
    return BiexpN_residual_fixedDs(data.first, params) - data.second;
}
double Biexp_residual_MSE_fixedDsf (const matrix<std::pair<double, double>, 0, 1>& data, const matrix<double,0,1>& params)
{
	double MSE = 0.0;
	double err;
	std::pair<double, double> pair;
	for(long i = 0; i < data.size(); i++) {
		pair = data(i);
		err = Biexp_residual_fixedDsf(pair, params);
		MSE = MSE + err*err;
	}
	if((int)data.size() == 0) return -1;
	MSE = MSE/((double)data.size());
	return MSE;
}

#pragma endregion

// ----------------------------------------------------------------------------------------
#pragma endregion

// ----------------------------------------------------------------------------------------
#pragma region MonoExponential, Normalized

double MonoExponentialN_model(const double& xvalue, const vectorParam1& params)
{
    double ADCm = params(0);
	double b = xvalue;

	double r =std::exp(-b*ADCm);
	return r;
}
double MonoExponentialN_model_mtx (const std::pair<double, double>& data, const vectorParam2& params)
{
    double ADCm = params(0);
	double b = data.first;
	double r =std::exp(-b*ADCm);
    return r - data.second;
}
double MonoExponentialN_residual_MSE_mtx (const matrix<std::pair<double, double>, 0, 1>& data, const matrix<double,0,1>& params)
{
	double MSE = 0.0;
	double err;
	std::pair<double, double> pair;
    cout << " ADC=" << params(0) << " err=";
	for(long i = 0; i < data.size(); i++) {
		pair = data(i);
		err = MonoExponentialN_model_mtx(pair, params);
        cout << "[" << data(i).first << "," << data(i).second << "]" << err;
		MSE = MSE + err*err;
	}
    cout << endl;
	if((int)data.size() == 0) return -1;
	MSE = MSE/((double)data.size());
	return MSE;
}
// This function is the "residual" for a least squares problem.   It takes an input/output
// pair and compares it to the output of our model and returns the amount of error.  The idea
// is to find the set of parameters which makes the residual small on all the data pairs.
double MonoExponentialN_residual (const std::pair<double, double>& data, const vectorParam1& params)
{
	return MonoExponentialN_model(data.first, params) - data.second;
}
double MonoExponentialN_residual_weights(const std::pair<double, double>& data, const vectorParam1& params)
{
	double r = (MonoExponentialN_model(bvalues[(int)data.first], params) - data.second)*weights[(int)data.first];
	return r;
}
double MonoExponentialN_residual_MSE (const matrix<std::pair<double, double>, 0, 1>& data, const matrix<double,0,1>& params)
{
	double MSE = 0.0;
	double err;
	std::pair<double, double> pair;
	if(VERBOSE > 1) {
		for(long i = 0; i < params.size(); i++) {
			cout << " " << params(i);
		}
		cout << endl;
	}
	for(long i = 0; i < data.size(); i++) {
		pair = data(i);
		err = MonoExponentialN_residual(pair, params);
		if(VERBOSE > 1)
			cout << " " << err;
		MSE = MSE + err*err;
	}
	if(VERBOSE > 1)
		cout << endl;
	if((int)data.size() == 0) return -1;
	MSE = MSE/((double)data.size());
#if VERBOSE==2
	cout << "monoN " << params(0) << " " << params(1) << " > " << MSE << endl;
#endif
	return MSE;
}
#pragma endregion

// ----------------------------------------------------------------------------------------
#pragma region Kurt, Normalized

double KurtN_model(const double& xvalue, const vectorParam2& params)
{
    double ADCk = params(0);
	double K = params(1);
	double b = xvalue;
	double r;
	r = std::exp(-b*ADCk+1.0/6.0*b*b*ADCk*ADCk*K);
	if(!is_finite(r)) {
		if(r > 0)
			r = std::numeric_limits<double>::max();
		else
			r = std::numeric_limits<double>::min();
	}
	return r;
}
// This function is the "residual" for a least squares problem.   It takes an input/output
// pair and compares it to the output of our model and returns the amount of error.  The idea
// is to find the set of parameters which makes the residual small on all the data pairs.
double KurtN_residual (
    const std::pair<double, double>& data,
    const vectorParam2& params
)
{
    return KurtN_model(data.first, params) - data.second;
}
//Mean Squared Error of data to the fitted curve
double KurtN_residual_MSE (
    const matrix<std::pair<double, double>, 0, 1>& data,
    const matrix<double,0,1>& params
)
{
	double MSE = 0.0;
	double err;
	std::pair<double, double> pair;
	if(VERBOSE > 1) {
		for(long i = 0; i < params.size(); i++) {
			cout << " " << params(i);
		}
		cout << endl;
	}
	for(long i = 0; i < data.size(); i++) {
		pair = data(i);
		err = KurtN_residual(pair, params);
		if(VERBOSE > 1)
			cout << " " << err;
		MSE = MSE + err*err;
	}
	if(VERBOSE > 1)
		cout << endl;
	if((int)data.size() == 0) return -1;
	MSE = MSE/((double)data.size());
	return MSE;
}

#pragma endregion

// ----------------------------------------------------------------------------------------
#pragma region Stretched, Normalized

double StretchedN_model(const double& xvalue, const vectorParam2& params)
{
    double ADCs = params(0);
	double Alpha = params(1);
	double b = xvalue;
	double ADCss;
	//explicit assertion that the base of exponential
	//is allways non-negative, because negative base are not defined for this model
	if(b*ADCs < 0)
		ADCss = 0;
	else
		ADCss = std::pow(b*ADCs, Alpha);
	double r = std::exp(-ADCss);
	return r;
}
// This function is the "residual" for a least squares problem.   It takes an input/output
// pair and compares it to the output of our model and returns the amount of error.  The idea
// is to find the set of parameters which makes the residual small on all the data pairs.
double StretchedN_residual (
    const std::pair<double, double>& data,
    const vectorParam2& params
)
{
    return StretchedN_model(data.first, params) - data.second;
}
//Mean Squared Error of data to the fitted curve
double StretchedN_residual_MSE (
    const matrix<std::pair<double, double>, 0, 1>& data,
    const matrix<double,0,1>& params
)
{
	double MSE = 0.0;
	double err;
	std::pair<double, double> pair;
	for(long i = 0; i < data.size(); i++) {
		pair = data(i);
		err = StretchedN_residual(pair, params);
		MSE = MSE + err*err;
	}
	if((int)data.size() == 0) return -1;
	MSE = MSE/((double)data.size());
	return MSE;
}

#pragma endregion

// ----------------------------------------------------------------------------------------
#pragma region Biexp, Normalized

double BiexpN_model(const double& xvalue, const vectorParam3& params)
{
    double f = params(0);
    double Df = params(1);
	double Ds = params(2);
	double b = xvalue;

	double r;
	r = f*std::exp(-b*Df)+(1-f)*std::exp(-b*Ds);
	return r;
}
// This function is the "residual" for a least squares problem.   It takes an input/output
// pair and compares it to the output of our model and returns the amount of error.  The idea
// is to find the set of parameters which makes the residual small on all the data pairs.
double BiexpN_residual (
    const std::pair<double, double>& data,
    const vectorParam3& params
)
{
    return BiexpN_model(data.first, params) - data.second;
}
double BiexpN_residual_weights (
    const std::pair<double, double>& data,
    const vectorParam3& params
)
{
	return (BiexpN_model(bvalues[(int)data.first], params) - data.second)*weights[(int)data.first];
}

//Mean Squared Error of data to the fitted curve
double BiexpN_residual_MSE (
    const matrix<std::pair<double, double>, 0, 1>& data,
    const matrix<double,0,1>& params
)
{
	double MSE = 0.0;
	double err;
	std::pair<double, double> pair;
	for(long i = 0; i < data.size(); i++) {
		pair = data(i);
		err = BiexpN_residual(pair, params);
		MSE = MSE + err*err;
	}
	if((int)data.size() == 0) return -1;
	MSE = MSE/((double)data.size());
	return MSE;
}

#pragma endregion

// ----------------------------------------------------------------------------------------

#pragma region Settings functions
/**
 * Resolves initializations
 * 
 * @param limits limits to be used, each items has {min,step,max}
 * @param pinit output list, each item has {p_1_i, p_2_i, ... p_N_i}, where N is number of parameters and i is initializations index
 * @param info info struture for printouts
 * @return return string description about initializations
 */
string resolveMultipleInitializations(const std::vector<std::vector<double> >& limits, std::vector<std::vector<double> >& pinits, const executionInfo& info, int info_i) {
	std::stringstream ss;
	//init indexes list
	std::stringstream outputs;
	int no_parameters = (int)limits.size();
	std::vector<double> p_is;
	outputs << info.modelName << " ";
	int no_constants = 0;
	int i = 0;
	for(i = 0; i < no_parameters; i++) {
		//Explicitly defined value list
		if(info.parameterLimit_explicit[info_i+i]) {
			outputs << info.parameterNames[info_i+i] << "=[";
			for(unsigned long val_i = 0; val_i < limits[i].size(); val_i++) {
				outputs << limits[i][val_i];
				if(val_i < (unsigned long)limits[i].size()-1) {
					outputs << " ";
				}
			}
			outputs << "]";
		//Limits, but the limits resolve into a single value
		} else if(limits[i][0] == limits[i][1] && limits[i][1] == limits[i][2]) {
			outputs << info.parameterNames[info_i+i] << "=[" << limits[i][0] << "] ";
			no_constants++;
		//Lo, step, hi limits
		} else if(!info.parameterLimit_explicit[info_i+i]) {
			outputs << info.parameterNames[info_i+i] << "=[" << limits[i][0] << ":" << limits[i][1] << ":" << limits[i][2] << "] ";
			if(limits[i][1] <= 0) {
				ss << "iteration for parameter " << info.parameterNames[info_i+i] << " in model " << info.modelName <<  " was not positive";
				throw ss.str();
			}
		}
		p_is.push_back(limits[i][0]);
	}
	cout << info.parameterNames[0] << "=[" << limits[0][0] << ":" << limits[0][1] << ":" << limits[0][2] << "] ";
	//keep adding values until all index counters reach their maximum value
	int unfinished = 1;
	//if all values are constants, set loop to add the low limit and then exit
	if(no_constants == no_parameters) 
		unfinished = 0;
	std::vector<int> explicit_p_iter(info.parameterLimit_explicit.size());
	while(1) {
		//resolve item that we are going to add
		std::vector<double> pinit;
		for(i = 0; i < no_parameters; i++) {
			pinit.push_back(p_is[i]);
		}
		//add the item to the list
		pinits.push_back(pinit);
		if(!unfinished)
			break;
		//iterate values for the next item, starting from the last
		for(i = no_parameters-1; i >= 0; i--) {
			//increase value in steps and in explicitly defined values
			if (info.parameterLimit_explicit[info_i+i]) {
				if(explicit_p_iter[i] < (int)limits[i].size()) {
					p_is[i] = limits[i][explicit_p_iter[i]];
					explicit_p_iter[i]++;
					break;
				//if maximum of highest level is reached, then we have done all initializations
				} else {
					if(i == 0) {
						unfinished = 0;
						break;
					} else {
						//if maximum is exceeded, reset to the first value and continue increasing higher levels
						explicit_p_iter[i] = 0;
						p_is[i] = limits[i][0];
					}
				}
			} else {
				p_is[i] += limits[i][1];
				//if maximum for this level is not exceeded, exit
				if(p_is[i] <= limits[i][2]) {
					break;
				//if maximum of highest level is exceeded, then we have done all initializations
				} else {
					if(i == 0) {
						unfinished = 0;
						break;
					} else {
						//if maximum is exceeded, reset to the first value and continue increasing higher levels
						p_is[i] = limits[i][0];
					}
				}
			}
		}
	}
	#if VERBOSE > 0
		for(unsigned int i = 0; i < pinits.size(); i++) {
			for(unsigned int j = 0; j < pinits[i].size(); j++) {
				cout << pinits[i][j] << " ";
			}
			cout << endl;
		}
	#endif
	return outputs.str();
}

/**
 * Resolve initializations from configuration
 *
 * @param cr (input) dlib configuration file reader
 * @param info information about the model: (input) modelName, (output) parameterNames
 * @return limits 
 */
std::vector<std::vector<double> > resolveLimits(const config_reader& cr, executionInfo& info) {

	std::vector<string> blocks;
	std::vector<string> pblocks;
	std::vector<string> keys;
	std::stringstream ss;

	//search for model name among blocks
	cr.get_blocks(blocks);
	VERBOSE1(((int)blocks.size()) << " blocks found, searhing for " << info.modelName);
	int model_i = std::find(blocks.begin(), blocks.end(), info.modelName) - blocks.begin();
	if(model_i == (int)blocks.size()) {
		ss << "Model " << info.modelName << " was not found in settings files";
		throw ss.str();
	}
	VERBOSE1("Resolving parameter limits for " << info.modelName);
	VERBOSE1(model_i << " " << blocks[model_i]);

	//read model blocks' sub blocks for parameters
	cr.block(blocks[model_i]).get_blocks(pblocks);
	if((int)pblocks.size() == 0) {
		ss << "Model " << info.modelName << " has no parameters defined";
		throw ss.str();
	}
	VERBOSE1(((int)pblocks.size()) << " parameters found");

	std::vector<std::vector<double> > limits((int)pblocks.size());
	try {
		info.parameterNames.resize((int)pblocks.size());
		info.parameterLimit_explicit.resize((int)pblocks.size());
		for(unsigned long b_i = 0; b_i < pblocks.size(); b_i++) {
			VERBOSE1("found parameter " << pblocks[b_i]);
			int pname_i = 0;
			//add parameter number
			cr.block(blocks[model_i]).block(pblocks[b_i]).get_keys(keys);
			int key_i = std::find(keys.begin(), keys.end(), "number") - keys.begin();
			if(key_i == (int)keys.size()) {
				ss <<  "ERROR: config file parameter " << pblocks[b_i] << " does not have parameter number [1..N]";	
				throw ss.str();
			} else {
				pname_i = atof(cr.block(blocks[model_i]).block(pblocks[b_i])[keys[key_i]].c_str());
			}
			if(pname_i < 0 || pname_i > (int)info.parameterNames.size()) {
				ss <<  "ERROR: config file parameter " << pblocks[b_i] << " has invalid parameter number " << pname_i << " outside [1.." << ((int)info.parameterNames.size()) << "]";	
				throw ss.str();
			}
			info.parameterNames[pname_i-1] = pblocks[b_i];
			info.parameterLimit_explicit[pname_i-1] = false;
            VERBOSE2("parameter " << pblocks[b_i] << " done");
		}

		//go trough each parameter 
		std::vector<double> p;
		//resolve parameters in order they are defined in info
		for (int pname_i = 0; pname_i < (int)info.parameterNames.size(); pname_i++) {
            VERBOSE2("resolving parameter " << info.parameterNames[pname_i]);

			//search for parameter name
			int i = std::find(pblocks.begin(), pblocks.end(), info.parameterNames[pname_i]) - pblocks.begin();
			if(i == (int)pblocks.size()) {
				ss << "ERROR: parameter name not found " << info.parameterNames[pname_i] << " does not have parameter number [1..N]";
				throw ss.str();
			}
			//resolve limit from key values
			std::vector<double> lim(3);
			cr.block(blocks[model_i]).block(pblocks[i]).get_keys(keys);
			if((int)keys.size() == 4) {
				int key_i = std::find(keys.begin(), keys.end(), "min") - keys.begin();
				if(key_i == (int)keys.size()) {
					ss << "ERROR: config file parameter " << info.parameterNames[pname_i] << " does not have min";				
					throw ss.str();
				} else {
					lim[0] = atof(cr.block(blocks[model_i]).block(pblocks[i])[keys[key_i]].c_str());
				}
				key_i = std::find(keys.begin(), keys.end(), "step") - keys.begin();
				if(key_i == (int)keys.size()) {
					ss << "ERROR: config file parameter " << info.parameterNames[pname_i] << " does not have step";				
					throw ss.str();
				} else {
					lim[1] = atof(cr.block(blocks[model_i]).block(pblocks[i])[keys[key_i]].c_str());
				}
				key_i = std::find(keys.begin(), keys.end(), "max") - keys.begin();
				if(key_i == (int)keys.size()) {
					ss << "ERROR: config file parameter " << info.parameterNames[pname_i] << " does not have max";				
					throw ss.str();
				} else {
					lim[2] = atof(cr.block(blocks[model_i]).block(pblocks[i])[keys[key_i]].c_str());			
				}
			} else if((int)keys.size() == 2) {
				int key_i = std::find(keys.begin(), keys.end(), "value") - keys.begin();
				if(key_i == (int)keys.size()) {
					ss << "ERROR: config file parameter " << info.parameterNames[pname_i] << " does not have value field";
					throw ss.str();
				} else {
					lim.clear();
					string strtok_s = cr.block(blocks[model_i]).block(pblocks[i])[keys[key_i]];
					string delimiters = " ";
					size_t current;
					size_t next = -1;
					do
					{
						current = next + 1;
						next = strtok_s.find_first_of( delimiters, current );
						lim.push_back(atof(strtok_s.substr( current, next - current ).c_str()));
					}
					while (next != string::npos);
					//mark this parameter initialization values to be explicitly defined
					info.parameterLimit_explicit[pname_i] = true;
				}
			} else {
				ss << "ERROR: invalid number of keys in " << info.modelName << " parameter block";
				throw ss.str();
			}
			//add parameter number
			int p_index = -1;
			int key_i = std::find(keys.begin(), keys.end(), "number") - keys.begin();
			if(key_i == (int)keys.size()) {
				ss << "ERROR: config file parameter " << info.parameterNames[pname_i] << " does not have parameter number [1..N]";				
				throw ss.str();
			} else {
				p_index = atof(cr.block(blocks[model_i]).block(pblocks[i])[keys[key_i]].c_str());	
			}
			limits[p_index-1] = lim;
            VERBOSE2("resolving parameter " << info.parameterNames[pname_i] << " done");
		}
	} catch(...) {
		throw;
	}
	return limits;
}

/**
 * Read key values from configuration
 */
template <typename T> std::vector<T> readKey(config_reader& cr, string keyname) {
	cout << "resolving " << keyname << endl;
	std::vector<T> output;
	std::vector<string> keys;
	cr.get_keys(keys);
	int key_i = std::find(keys.begin(), keys.end(), keyname) - keys.begin();
	if(key_i == (int)keys.size()) {
		THROW_FATAL("ERROR: config file parameter does not have min")
	} else {
		//parse b-values with tokenizer
		dlib::tokenizer::kernel_1a_c tokzr;
		istringstream sin;
		sin.str(cr[keys[key_i]]);
		tokzr.set_stream(sin);
		int type = 0;
		string token;
		cout << keyname << "=[";
		string previous_double_str = "";
		while(1) {
			tokzr.get_token(type, token);
			if(type == dlib::tokenizer::kernel_1a_c::END_OF_LINE || type == dlib::tokenizer::kernel_1a_c::END_OF_FILE) break;
			if(type == dlib::tokenizer::kernel_1a_c::WHITE_SPACE) {
				cout << token << " ";
				if((int)previous_double_str.size() > 0) {
					double val = atof(previous_double_str.c_str());
					cout << val << " ";
					output.push_back(val);
					previous_double_str = "";
				}
			} else if(type == dlib::tokenizer::kernel_1a_c::NUMBER || 
				(type == dlib::tokenizer::kernel_1a_c::CHAR && token.compare(".") == 0)) {
				previous_double_str += token;
			}
		}
		if((int)previous_double_str.size() > 0) {
			double val = atof(previous_double_str.c_str());
			cout << val << " ";
			output.push_back(val);
		}
		cout << "]" << endl;
	}	
	return output;
} 

#pragma endregion
