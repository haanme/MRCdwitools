#pragma once

#include <dlib/optimization.h>
#include <dlib/rand.h>

using namespace dlib;
using namespace std;

/*
 * Gets environment variable for QC purposes
 * 
 * @param outputs stream where variable is written
 * @param name name of variable in current OS
 * @returns 0 for success, 1 for failure in getting the variable
 */
int getEnvVar(stringstream& outputs, std::string const & name) {
	char * pOS;	
	pOS = std::getenv ( name.c_str() );
	if(pOS!=NULL) {
		outputs << " " << pOS;
		return 0;
	} else {
		outputs << " <N/A>";
		return 1;
	}
}

#ifndef OS_STR
#define OS_STR "UNDEFINED COMPILATION OS STRING"
#endif

//Machine information for this binary
string getCompilationInfoString() {
	std::stringstream outputs;
	outputs << " (OS[";
	outputs << (OS_STR);
	outputs << "]";
#if defined(WIN32)
	outputs << " Win32";
	getEnvVar(outputs, "OS");
	outputs << " DOUBLE DEF[";
	outputs << " min:" << (DBL_MIN);
	outputs << " max:" << DBL_MAX;
	outputs << " radix:" << DBL_RADIX;
	outputs << " eps:" << DBL_EPSILON;
	outputs << "])";
#elif defined(WIN64)
	outputs << " Win64";
	outputs << " " << getenv("OS");
	outputs << " DOUBLE DEF[";
	outputs << " min:" << (DBL_MIN);
	outputs << " max:" << DBL_MAX;
	outputs << " radix:" << DBL_RADIX;
	outputs << " eps:" << DBL_EPSILON;
	outputs << "])";
#elif defined(__CYGWIN__)
	outputs << " Cygwin";
	outputs << " DOUBLE DEF[";
	outputs << " min:" << std::numeric_limits<double>::min;
	outputs << " max:" << std::numeric_limits<double>::max;
	outputs << " radix:" << std::numeric_limits<double>::radix;
	outputs << " eps:" << std::numeric_limits<double>::epsilon();
	outputs << "])";
#elif defined(unix)
	#include <limits>
	outputs << " Unix";
	outputs << " DOUBLE DEF[";
	outputs << " min:" << std::numeric_limits<double>::min;
	outputs << " max:" << std::numeric_limits<double>::max;
	outputs << " radix:" << std::numeric_limits<double>::radix;
	outputs << " radix digits:" << std::numeric_limits<double>::digits;
	outputs << " epsilon:" << std::numeric_limits<double>::epsilon();
	outputs << "])";
#elif defined(__linux__)
	outputs << " Linux";
	outputs << " DOUBLE DEF[";
	outputs << " min:" << (DBL_MIN);
	outputs << " max:" << DBL_MAX;
	outputs << " radix:" << DBL_RADIX;
	outputs << " eps:" << DBL_EPSILON;
	outputs << "])";
#elif defined(__APPLE__) && defined(__MACH__)
	outputs << " OSX";
	outputs << " DOUBLE DEF[";
	outputs << " min:" << (__DBL_MIN__);
	outputs << " max:" << __DBL_MAX__;
	outputs << " eps:" << __DBL_EPSILON__;
	outputs << "])";
#else
	outputs << " Unrecognized OS"
	outputs << " DOUBLE DEF[";
	outputs << " min:" << (__DBL_MIN__);
	outputs << " max:" << __DBL_MAX__;
	outputs << " radix:" << __DBL_RADIX__;
	outputs << " eps:" << __DBL_EPSILON__;
	outputs << "])";
#endif
	return outputs.str();
}

//vector lengths predefined because we need them already at compilation time
#define DEFINE_VECTORP(VECTOR_LENGTH, NAME) \
	typedef matrix<double,VECTOR_LENGTH,1> NAME;
DEFINE_VECTORP(1, vectorParam1)
DEFINE_VECTORP(2, vectorParam2)
DEFINE_VECTORP(3, vectorParam3)
DEFINE_VECTORP(4, vectorParam4)
DEFINE_VECTORP(5, vectorParam5)
DEFINE_VECTORP(6, vectorParam6)
DEFINE_VECTORP(7, vectorParam7)
DEFINE_VECTORP(8, vectorParam8)
DEFINE_VECTORP(9, vectorParam9)
DEFINE_VECTORP(10, vectorParam10)
DEFINE_VECTORP(11, vectorParam11)
DEFINE_VECTORP(12, vectorParam12)
DEFINE_VECTORP(13, vectorParam13)
DEFINE_VECTORP(14, vectorParam14)
DEFINE_VECTORP(15, vectorParam15)
DEFINE_VECTORP(16, vectorParam16)
DEFINE_VECTORP(17, vectorParam17)
DEFINE_VECTORP(18, vectorParam18)
DEFINE_VECTORP(19, vectorParam19)
DEFINE_VECTORP(20, vectorParam20)
#define CASE_VECTORP(VECTOR_LENGTH, VARIABLENAME) \
	case VECTOR_LENGTH: \
		matrix<double,VECTOR_LENGTH,1> VARIABLENAME;

//model names
enum ModelName {
	//Mono-exponential
	Mono = 1,
	//Kurtosis model
	Kurt = 2,
	//Stretched exponential
	Strecthed = 3,
	//Bi-exponential model
	Biexp = 4,
	//Segmented analysis presented by e. g. Marzi et al.
	MarziSegmentedBiexp	= 5,
	//Asymptotic analysis, see e.g. Pekar 1992
	AsymptoticBiexp = 6
};

//operation modes for fittings and simulations, as defines for compilations in non-VS environments
#define OP_BFGS 1
#define OP_LM 2
#define	OP_SIMULATION 3
#define OP_BFGS_NORMALIZED 4
#define OP_LM_NORMALIZED 5
//model names as strings as definitions
#define MONO "Mono"
#define T1RHO "T1rho"
#define KURT "Kurt"
#define STRETCHED "Stretched"
#define BIEXP "Biexp"
#define MARZISEGMENTEDBIEXP "MarziSegmentedBiexp"
#define CHOSEGMENTEDBIEXP "ChoSegmentedBiexp"
#define ASYMPTOTICBIEXP "AsymptoticBiexp"
#define LOGLINASYMPTOTICBIEXP "LogLinAsymptoticBiexp"
#define NGIVIM "NonGaussianIVIM"
#define LNGIVIM "LinNonGaussianIVIM"
#define SLNGIVIM "SLinNonGaussianIVIM"
#define LOGINITIVIM "LogInitBiexp"
//precision of output values written to file
#define OUTPUT_PRECISION 10
//step length between writing to screen
#define PRINTOUT_STEP 1000
#define MAX_NUMBER_OF_ITERATIONS 100
#define DELTA_STOP 1e-7
#define VERBOSE 0
#define VERBOSE1(TEXT) if(VERBOSE > 0) cout << TEXT << endl;
#define VERBOSE2(TEXT) if(VERBOSE > 1) cout << TEXT << endl;

#define OPT_VERBOSE
#if VERBOSE==1
#define OPT_VERBOSE .be_verbose()
#endif
#if VERBOSE==2
#define OPT_VERBOSE .be_verbose()
#endif
// ----------------------------------------------------------------------------------------
//execution info
struct fitInfo {
	//model name
	string name;
	// parameter names
	std::vector<string> parameterNames;
	//list of execution arguments for each parameter as {min value, step, max value}
	std::vector<std::vector<double> > multiple_initialization_values;
};

//execution info
struct executionInfo {
	//subwindow coordinates in original data
	std::vector<int> subwindow;
	//number of b-values used
	int number;
	//b-values used
	std::vector<double> bset;
	//b-value weightings
	std::vector<double> bset_weights;
	//ROI slices in original data
	std::vector<int> ROIslice;
	//ROI name
	string ROIName;
	//execution time in seconds
	long executiontime;
	// model name
	string modelName;
	// parameter names
	std::vector<string> parameterNames;
	// parameter initialization type is limits [lo,step,hi] or explicitly defined values [N1, N2, ... N3]
	std::vector<bool> parameterLimit_explicit;
	//input filename
	string input_filename;
	//description about fit
	string desc;
};

// ----------------------------------------------------------------------------------------
#define THROW_FATAL(TEXT) \
cout << TEXT << endl; \
throw TEXT; \

// ----------------------------------------------------------------------------------------
//fixed value pairs to be used in the models, actual parameters depend of currently fitted model
// ----------------------------------------------------------------------------------------
std::vector<double> fixed_parameters;
// ----------------------------------------------------------------------------------------
//Macro to define LM fit with fixed number of parameter values, implementation
#define FIT_LM_MI_IMPLEMENTATION(NAME, NUMBER_OF_PARAMETERS, NUMBER_OF_BVALUES) \
	void NAME(std::vector<double> bvalues, const std::vector<std::vector<double> >& SIdata, const std::vector<std::vector<double> >& pinit, double (*residual_fun)(const std::pair<double, double>&, const matrix<double,NUMBER_OF_PARAMETERS,1>&), std::vector<std::vector<double> >& p_vals, executionInfo& info, std::vector<double>& C_values) { \
		if((int)SIdata.size() == 0) \
			throw "SIdata length was 0"; \
		matrix<std::pair<double, double>, NUMBER_OF_BVALUES, 1> data_samples; \
		for(int b_i = 0; b_i < NUMBER_OF_BVALUES; b_i++) data_samples(b_i).first = bvalues[b_i]; \
		cout << ((int)pinit.size()) << " initializations " << NUMBER_OF_PARAMETERS << " for parameters" << endl; \
		const int no_parameters = NUMBER_OF_PARAMETERS; \
		int no_bvalues = (int)SIdata[0].size(); \
		cout << no_bvalues << " b-values" << endl; \
		double RMSE; \
		double RMSE_min; \
		double absolute_error; \
		matrix<double,NUMBER_OF_PARAMETERS,1> x; \
		matrix<double,NUMBER_OF_PARAMETERS,1> x_min; \
		p_vals.clear(); \
		clock_t startTime = clock(); \
		unsigned int p_i; \
		unsigned int c_i; \
		double C; \
		double SImax; \
		double SImin; \
		double SIrange; \
		for(unsigned int si_i = 0; si_i < SIdata.size(); si_i++) { \
			RMSE_min = std::numeric_limits<double>::max(); \
			std::fill(x_min.begin(), x_min.end(), std::numeric_limits<double>::quiet_NaN()); \
			SImax = std::numeric_limits<double>::min(); \
			SImin = std::numeric_limits<double>::max(); \
			for(int b_i = 0; b_i < no_bvalues; b_i++) { \
				data_samples(b_i).second = SIdata[si_i][b_i]; \
				if(data_samples(b_i).second < SImin) SImin = data_samples(b_i).second; \
				if(data_samples(b_i).second > SImax) SImax = data_samples(b_i).second; \
			} \
			SIrange = SImax - SImin; \
			for(c_i = 0; c_i < (int)C_values.size(); c_i++) { \
				C = SImin+C_values[c_i]*SIrange; \
				for(p_i = 0; p_i < (int)pinit.size(); p_i++) { \
					for(int param_i = 0; param_i < no_parameters-1; param_i++) { \
						x(param_i) = pinit[p_i][param_i]; \
					} \
					x(no_parameters-1) = C; \
					solve_least_squares_lm(objective_delta_stop_strategy(DELTA_STOP, MAX_NUMBER_OF_ITERATIONS)OPT_VERBOSE, residual_fun, derivative(residual_fun),	data_samples, x); \
					RMSE = 0; \
					for(int j = 0; j < no_bvalues; j++) { \
						absolute_error = residual_fun(data_samples(j), x); \
						RMSE += std::pow(absolute_error, 2.0); \
					} \
					RMSE /= (double)no_bvalues; \
					RMSE = std::sqrt(RMSE); \
					if(RMSE < RMSE_min) { \
						RMSE_min = RMSE; \
						x_min  = x; \
					} \
				} \
			} \
			std::vector<double> p_val = std::vector<double>(no_parameters+1); \
			if(si_i % PRINTOUT_STEP == 0) { \
				cout << "[" << si_i << "]:"; \
				for(p_i = 0; p_i < no_parameters; p_i++) { \
					cout << info.parameterNames[p_i] << "=" << x_min(p_i) << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
				} \
				cout << " RMSE="<< RMSE_min << std::setprecision(OUTPUT_PRECISION) << std::fixed << endl; \
			} \
			for(p_i = 0; p_i < no_parameters; p_i++) { \
				p_val[p_i] = x_min(p_i); \
			} \
			p_val[no_parameters] = RMSE_min; \
			p_vals.push_back(p_val); \
		} \
		info.executiontime = double( clock() - startTime ) / (double)CLOCKS_PER_SEC; \
		cout << "Execution finished in " << info.executiontime << " seconds" << endl; \
	} \

//Macro to define LM fit with fixed number of parameter values, implementation
#define FIT_LM_MI_IMPLEMENTATION_BIEXP(NAME, NUMBER_OF_PARAMETERS, NUMBER_OF_BVALUES) \
	void NAME(std::vector<double> bvalues, const std::vector<std::vector<double> >& SIdata, const std::vector<std::vector<double> >& pinit, double (*residual_fun)(const std::pair<double, double>&, const matrix<double,NUMBER_OF_PARAMETERS,1>&), std::vector<std::vector<double> >& p_vals, executionInfo& info, std::vector<double>& C_values) { \
		if((int)SIdata.size() == 0) \
			throw "SIdata length was 0"; \
		matrix<std::pair<double, double>, NUMBER_OF_BVALUES, 1> data_samples; \
		for(int b_i = 0; b_i < NUMBER_OF_BVALUES; b_i++) data_samples(b_i).first = bvalues[b_i]; \
		cout << ((int)pinit.size()) << " initializations for "<< NUMBER_OF_PARAMETERS << " parameters" << endl; \
		const int no_parameters = NUMBER_OF_PARAMETERS; \
		int no_bvalues = (int)SIdata[0].size(); \
		cout << no_bvalues << " b-values" << endl; \
		double RMSE; \
		double RMSE_min; \
		double absolute_error; \
		matrix<double,NUMBER_OF_PARAMETERS,1> x; \
		matrix<double,NUMBER_OF_PARAMETERS,1> x_min; \
		p_vals.clear(); \
		clock_t startTime = clock(); \
		unsigned int p_i; \
		unsigned int c_i; \
		double C; \
		double SImax; \
		double SImin; \
		double SIrange; \
		for(unsigned int si_i = 0; si_i < SIdata.size(); si_i++) { \
			RMSE_min = std::numeric_limits<double>::infinity(); \
			std::fill(x_min.begin(), x_min.end(), std::numeric_limits<double>::quiet_NaN()); \
			SImax = std::numeric_limits<double>::min(); \
			SImin = std::numeric_limits<double>::max(); \
			for(int b_i = 0; b_i < no_bvalues; b_i++) { \
				data_samples(b_i).second = SIdata[si_i][b_i]; \
				if(data_samples(b_i).second < SImin) SImin = data_samples(b_i).second; \
				if(data_samples(b_i).second > SImax) SImax = data_samples(b_i).second; \
			} \
			SIrange = SImax - SImin; \
			for(c_i = 0; c_i < (int)C_values.size(); c_i++) { \
				C = SImin+C_values[c_i]*SIrange; \
				for(p_i = 0; p_i < (int)pinit.size(); p_i++) { \
					for(int param_i = 0; param_i < no_parameters; param_i++) \
						x(param_i) = pinit[p_i][param_i]; \
					x(no_parameters-1) = C; \
					solve_least_squares(objective_delta_stop_strategy(DELTA_STOP, MAX_NUMBER_OF_ITERATIONS)OPT_VERBOSE, residual_fun, derivative(residual_fun),	data_samples, x); \
					if(x(1) < x(2)) { \
						double x_swap = x(1); \
						x(1) = x(2); \
						x(2) = x_swap; \
						x(0) = 1-x(1); \
					} \
					RMSE = 0; \
					for(int j = 0; j < no_bvalues; j++) { \
						absolute_error = residual_fun(data_samples(j), x); \
						RMSE += std::pow(absolute_error, 2.0); \
					} \
					RMSE /= (double)no_bvalues; \
					RMSE = std::sqrt(RMSE); \
					if(RMSE < RMSE_min) { \
						RMSE_min = RMSE; \
						x_min  = x; \
					} \
				} \
			} \
			std::vector<double> p_val = std::vector<double>(no_parameters+1); \
			if(si_i % PRINTOUT_STEP == 0) { \
				cout << "[" << si_i << "]:"; \
				for(p_i = 0; p_i < no_parameters; p_i++) { \
					cout << info.parameterNames[p_i] << "=" << x_min(p_i) << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
				} \
				cout << " RMSE="<< RMSE_min << std::setprecision(OUTPUT_PRECISION) << std::fixed << endl; \
			} \
			for(p_i = 0; p_i < no_parameters; p_i++) { \
				p_val[p_i] = x_min(p_i); \
			} \
			p_val[no_parameters] = RMSE_min; \
			p_vals.push_back(p_val); \
		} \
		info.executiontime = double( clock() - startTime ) / (double)CLOCKS_PER_SEC; \
		cout << "Execution finished in " << info.executiontime << " seconds" << endl; \
	} \

//Macro to define LM fit with fixed number of parameter values, implementation
#define FIT_LM_MI_IMPLEMENTATION_MARZISEGMENTEDBIEXP(NAME, NUMBER_OF_BVALUES, NUMBER_OF_TAIL_BVALUES) \
	void NAME(std::vector<double> bvalues, const std::vector<std::vector<double> >& SIdata, const std::vector<std::vector<double> >& pinit, double (*residual_fun)(const std::pair<double, double>&, const matrix<double,2,1>&), double (*residual_fun_2nd)(const std::pair<double, double>&, const matrix<double,3,1>&), std::vector<std::vector<double> >& p_vals, executionInfo& info, std::vector<double>& C_values) { \
		if((int)SIdata.size() == 0) \
			throw "SIdata length was 0"; \
		matrix<std::pair<double, double>, NUMBER_OF_BVALUES, 1> data_samples; \
		for(int b_i = 0; b_i < NUMBER_OF_BVALUES; b_i++) data_samples(b_i).first = bvalues[b_i]; \
		matrix<std::pair<double, double>, NUMBER_OF_TAIL_BVALUES, 1> data_tail_samples; \
		for(int b_i = NUMBER_OF_BVALUES-NUMBER_OF_TAIL_BVALUES; b_i < NUMBER_OF_BVALUES; b_i++) data_tail_samples(b_i-(NUMBER_OF_BVALUES-NUMBER_OF_TAIL_BVALUES)).first = bvalues[b_i]; \
		cout << ((int)pinit.size()) << " initializations for 1+3 parameters" << endl; \
		const int no_parameters = 4; \
		int no_bvalues = (int)SIdata[0].size(); \
		cout << no_bvalues << " b-values" << endl; \
		double RMSE; \
		double RMSE_min; \
		double absolute_error; \
		matrix<double,2,1> x_1st; \
		matrix<double,3,1> x_2nd; \
		matrix<double,4,1> x; \
		matrix<double,4,1> x_min; \
		p_vals.clear(); \
		clock_t startTime = clock(); \
		unsigned int p_i; \
		fixed_parameters.clear(); \
		fixed_parameters.push_back(0); \
		unsigned int c_i; \
		double C; \
		double SImax; \
		double SImin; \
		double SIrange; \
		for(unsigned int si_i = 0; si_i < (int)SIdata.size(); si_i++) { \
			RMSE_min = std::numeric_limits<double>::infinity(); \
			SImax = std::numeric_limits<double>::min(); \
			SImin = std::numeric_limits<double>::max(); \
			for(int b_i = 0; b_i < no_bvalues; b_i++) { \
				data_samples(b_i).second = SIdata[si_i][b_i]; \
				if(data_samples(b_i).second < SImin) SImin = data_samples(b_i).second; \
				if(data_samples(b_i).second > SImax) SImax = data_samples(b_i).second; \
			} \
			SIrange = SImax - SImin; \
			std::fill(x_min.begin(), x_min.end(), std::numeric_limits<double>::quiet_NaN()); \
			for(int b_i = 0; b_i < no_bvalues; b_i++) \
				data_samples(b_i).second = SIdata[si_i][b_i]; \
			for(int b_i = no_bvalues-NUMBER_OF_TAIL_BVALUES; b_i < no_bvalues; b_i++) \
				data_tail_samples(b_i-(no_bvalues-NUMBER_OF_TAIL_BVALUES)).second = SIdata[si_i][b_i]; \
			for(c_i = 0; c_i < (int)C_values.size(); c_i++) { \
				C = SImin+C_values[c_i]*SIrange; \
				for(p_i = 0; p_i < (int)pinit.size(); p_i++) { \
					x_1st(0) = pinit[p_i][2]; \
					x_1st(1) = data_tail_samples(0).second; \
					solve_least_squares_lm(objective_delta_stop_strategy(DELTA_STOP, MAX_NUMBER_OF_ITERATIONS)OPT_VERBOSE, residual_fun, derivative(residual_fun),	data_tail_samples, x_1st); \
					if(x_1st(0) < 0) x_1st(0) = 0; \
					x_2nd(0) = pinit[p_i][0]; \
					x_2nd(1) = pinit[p_i][1]; \
					fixed_parameters[0] = x_1st(0); \
					x_2nd(2) = C; \
					solve_least_squares_lm(objective_delta_stop_strategy(DELTA_STOP, MAX_NUMBER_OF_ITERATIONS)OPT_VERBOSE, residual_fun_2nd, derivative(residual_fun_2nd),	data_samples, x_2nd); \
					if(x_2nd(1) < x_1st(0)) { \
						double x_swap = x_2nd(1); \
						x_2nd(1) = x_1st(0); \
						x_2nd(1) = x_swap; \
						x_2nd(0) = 1-x_2nd(0); \
					} \
					x(0) = x_2nd(0); \
					x(1) = x_2nd(1); \
					x(2) = x_1st(0); \
					x(3) = x_2nd(2); \
					RMSE = 0; \
					for(int j = 0; j < no_bvalues; j++) { \
						absolute_error = residual_fun_2nd(data_samples(j), x_2nd); \
						RMSE += std::pow(absolute_error, 2.0); \
					} \
					RMSE /= (double)no_bvalues; \
					RMSE = std::sqrt(RMSE); \
					if(RMSE < RMSE_min) { \
						RMSE_min = RMSE; \
						x_min  = x; \
					} \
				} \
			} \
			std::vector<double> p_val = std::vector<double>(no_parameters+1); \
			if(si_i % PRINTOUT_STEP == 0) { \
				cout << "[" << si_i << "]:"; \
				for(p_i = 0; p_i < no_parameters; p_i++) { \
					cout << info.parameterNames[p_i] << "=" << x_min(p_i) << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
				} \
				cout << " RMSE="<< RMSE_min << std::setprecision(OUTPUT_PRECISION) << std::fixed << endl; \
			} \
			for(p_i = 0; p_i < no_parameters; p_i++) { \
				p_val[p_i] = x_min(p_i); \
			} \
			p_val[no_parameters] = RMSE_min; \
			p_vals.push_back(p_val); \
		} \
		info.executiontime = double( clock() - startTime ) / (double)CLOCKS_PER_SEC; \
		cout << "Execution finished in " << info.executiontime << " seconds" << endl; \
	} \

//Macros for Trust-Region reflective fitting related functions
#pragma region Trust-Region reflective macros

// ----------------------------------------------------------------------------------------
typedef matrix<double,0,1> column_vector;
//data points, defined here only once for all following models
matrix<std::pair<double, double>, 0, 1> data_samples;
matrix<std::pair<double, double>, 0, 1> data_samples2;
//function that is optimized
double (*TR_funct)(const matrix<std::pair<double, double>, 0, 1>&, const matrix<double,0,1>&);
double (*TR_funct2)(const matrix<std::pair<double, double>, 0, 1>&, const matrix<double,0,1>&);
double eps = 1e-6;
//evaluation for BFGS fitting purposes
double TR_evaluation(const column_vector& m)
{
	double r = TR_funct(data_samples, m);
	if(!is_finite(r)) {
		if(r > 0)
			r = std::numeric_limits<double>::max();
		else
			r = std::numeric_limits<double>::min();
	}
	return r;
}
double TR_evaluation2(const column_vector& m)
{
	double r = TR_funct2(data_samples2, m);
	if(!is_finite(r)) {
		if(r > 0)
			r = std::numeric_limits<double>::max();
		else
			r = std::numeric_limits<double>::min();
	}
	return r;
}
//Hessian matrix, inside TR_model which is currently unused
matrix<double> TR_hessian (const column_vector& x) {
		double h1;
		double h2;
		//temporary variables for storing neighbouring derivatives and parameter values
		column_vector x_eps_1((int)x.size());
		column_vector x_eps_2((int)x.size());
		column_vector x_eps_3((int)x.size());
		column_vector x_eps_4((int)x.size());
		column_vector funct_val(4);
		matrix<double> h_r((int)x.size(), (int)x.size());
		//calculate 1st derivative and 2nd devivative to the diagonal diagonal
		for(int i = 0; i < (int)x.size(); i++) {
			h1 = eps;
			x_eps_1 = x;
			x_eps_2 = x;
			x_eps_1(i) = x(i) + h1;
			x_eps_2(i) = x(i) - h1;
			//calculate derivative in respect to i:th parameter
			funct_val(0) = TR_funct(data_samples, x_eps_1);
			funct_val(1) = TR_funct(data_samples, x_eps_2);
			funct_val(2) = TR_funct(data_samples, x);
			//calculate 2nd derivative in respect to i:th parameter
	        h_r(i,i) = (funct_val(1)-2*funct_val(2)+funct_val(0))/(h1*h1);
#if VERBOSE==2
			cout << std::setprecision(OUTPUT_PRECISION) << std::fixed; 
			cout << "f" << i << "=" << funct_val(2) << endl; 
			cout << "h" << i << "=" << h_r(i,i) << endl; 
#endif
		}
		//fill Hessian matrix with partial derivatives
		for(int i = 0; i < (int)x.size()-1; i++) {
			for(int j = i+1; j < (int)x.size(); j++) {
				h1 = eps;
				h2 = eps;
				x_eps_1 = x;
				x_eps_2 = x;
				x_eps_3 = x;
				x_eps_4 = x;
				x_eps_1(i) = x(i) + h1;
				x_eps_1(j) = x(j) + h2;
				x_eps_2(i) = x(i) + h1;
				x_eps_2(j) = x(j) - h2;
				x_eps_3(i) = x(i) - h1;
				x_eps_3(j) = x(j) + h2;
				x_eps_4(i) = x(i) - h1;
				x_eps_4(j) = x(j) - h2;
				funct_val(0) = TR_funct(data_samples, x_eps_1);
				funct_val(1) = TR_funct(data_samples, x_eps_2);				
				funct_val(2) = TR_funct(data_samples, x_eps_3);
				funct_val(3) = TR_funct(data_samples, x_eps_4);				
		        h_r(j,i) = h_r(i,j) = (funct_val(0)-funct_val(1)-funct_val(2)+funct_val(3))/(4*h1*h2);
			}
		}
#if VERBOSE==2
		cout << std::setprecision(OUTPUT_PRECISION) << std::fixed << endl; 
		for(int i = 0; i < (int)x.size(); i++) {
			for(int j = 0; j < (int)x.size(); j++) {
				cout << h_r(i,j) << " ";
			}
			cout << endl;
		}
		cout << endl;
#endif
		return h_r;
	}
//Derivative for function
const column_vector TR_derivative (const column_vector& x) {
		double h1;
		column_vector d_r((int)x.size());
		column_vector x_eps_1((int)x.size());
		column_vector x_eps_2((int)x.size());
		column_vector funct_val(3);
		//calculate 1st derivative and 2nd devivative to the diagonal diagonal
		for(int i = 0; i < (int)x.size(); i++) {
			h1 = eps*x(i);
			x_eps_1 = x;
			x_eps_2 = x;
			x_eps_1(i) = x(i) + h1;
			x_eps_2(i) = x(i) - h1;
			//calculate derivative in respect to i:th parameter
			funct_val(0) = TR_funct(data_samples, x_eps_1);
			funct_val(1) = TR_funct(data_samples, x_eps_2);
			funct_val(2) = TR_funct(data_samples, x);
			d_r(i) = (funct_val(0)-funct_val(1))/(2*h1);
			//calculate 2nd derivative in respect to i:th parameter
#if VERBOSE==2
			cout << std::setprecision(OUTPUT_PRECISION) << std::fixed; 
			cout << "f" << i << "=" << funct_val(2) << endl; 
			cout << "d" << i << "=" << d_r(i) << endl; 
#endif
		}
		return d_r;
	}
//Derivative for function 2
const column_vector TR_derivative2 (const column_vector& x) {
		double h1;
		column_vector d_r((int)x.size());
		column_vector x_eps_1((int)x.size());
		column_vector x_eps_2((int)x.size());
		column_vector funct_val(3);
		//calculate 1st derivative and 2nd devivative to the diagonal diagonal
		for(int i = 0; i < (int)x.size(); i++) {
			h1 = eps*x(i);
			x_eps_1 = x;
			x_eps_2 = x;
			x_eps_1(i) = x(i) + h1;
			x_eps_2(i) = x(i) - h1;
			//calculate derivative in respect to i:th parameter
			funct_val(0) = TR_funct2(data_samples, x_eps_1);
			funct_val(1) = TR_funct2(data_samples, x_eps_2);
			funct_val(2) = TR_funct2(data_samples, x);
			d_r(i) = (funct_val(0)-funct_val(1))/(2*h1);
			//calculate 2nd derivative in respect to i:th parameter
#if VERBOSE==2
			cout << std::setprecision(OUTPUT_PRECISION) << std::fixed; 
			cout << "f" << i << "=" << funct_val(2) << endl; 
			cout << "d" << i << "=" << d_r(i) << endl;
#endif
		}
		return d_r;
	}
#pragma endregion
// ----------------------------------------------------------------------------------------

    template < typename funct >  double backtracking_line_search_fixed (const funct& f, double f0, double d0, double alpha, double rho, unsigned long max_iter )
    {
        DLIB_ASSERT (
            0 < rho && rho < 1 && max_iter > 0,
            "\tdouble backtracking_line_search()"
            << "\n\tYou have given invalid arguments to this function"
            << "\n\t rho:      " << rho
            << "\n\t max_iter: " << max_iter
        );

        // make sure alpha is going in the right direction.  That is, it should be opposite
        // the direction of the gradient.
        if ((d0 > 0 && alpha > 0) ||
            (d0 < 0 && alpha < 0))
        {
            alpha *= -1;
        }

        bool have_prev_alpha = false;
        double prev_alpha = 0;
        double prev_val = 0;
        unsigned long iter = 0;
        while (true)
        {
            ++iter;
            const double val = f(alpha);
            if (val <= f0 + alpha*rho*d0 || iter >= max_iter)
            {
                return alpha;
            }
            else
            {
                // Interpolate a new alpha.  We also make sure the step by which we
                // reduce alpha is not super small.
                double step;
                if (!have_prev_alpha)
                {
                    if (d0 < 0)
                        step = alpha*put_in_range(0.1,0.9, poly_min_extrap(f0, d0, val));
                    else
                        step = alpha*put_in_range(0.1,0.9, poly_min_extrap(f0, -d0, val));
                    have_prev_alpha = true;
                }
                else
                {
                    if (d0 < 0)
                        step = put_in_range(0.1*alpha,0.9*alpha, poly_min_extrap(f0, d0, alpha, val, prev_alpha, prev_val));
                    else
                        step = put_in_range(0.1*alpha,0.9*alpha, -poly_min_extrap(f0, -d0, -alpha, val, -prev_alpha, prev_val));
                }
                if(step < DELTA_STOP)
                    return DELTA_STOP;
                prev_alpha = alpha;
                prev_val = val;

                alpha = step;
            }
        }
    }
// ----------------------------------------------------------------------------------------
    template < typename search_strategy_type, typename stop_strategy_type, typename funct, typename funct_der, typename T, typename EXP1, typename EXP2 >
    double find_min_box_constrained_fixed (search_strategy_type search_strategy, stop_strategy_type stop_strategy, const funct& f, const funct_der& der, 
        T& x, const matrix_exp<EXP1>& x_lower, const matrix_exp<EXP2>& x_upper )
    {
        /*
            The implementation of this function is more or less based on the discussion in
            the paper Projected Newton-type Methods in Machine Learning by Mark Schmidt, et al.
        */

        // make sure the requires clause is not violated
        COMPILE_TIME_ASSERT(is_matrix<T>::value);
        // The starting point (i.e. x) must be a column vector.  
        COMPILE_TIME_ASSERT(T::NC <= 1);

        DLIB_ASSERT (
            is_col_vector(x) && is_col_vector(x_lower) && is_col_vector(x_upper) &&
            x.size() == x_lower.size() && x.size() == x_upper.size(),
            "\tdouble find_min_box_constrained()"
            << "\n\t The inputs to this function must be equal length column vectors."
            << "\n\t is_col_vector(x):       " << is_col_vector(x)
            << "\n\t is_col_vector(x_upper): " << is_col_vector(x_upper)
            << "\n\t is_col_vector(x_upper): " << is_col_vector(x_upper)
            << "\n\t x.size():               " << x.size()
            << "\n\t x_lower.size():         " << x_lower.size()
            << "\n\t x_upper.size():         " << x_upper.size()
        );
        DLIB_ASSERT (
            min(x_upper-x_lower) > 0,
            "\tdouble find_min_box_constrained()"
            << "\n\t You have to supply proper box constraints to this function."
            << "\n\r min(x_upper-x_lower): " << min(x_upper-x_lower)
        );


        T g, s;
        double f_value = f(x);
        g = der(x);

        DLIB_ASSERT(is_finite(f_value), "The objective function generated non-finite outputs");
        DLIB_ASSERT(is_finite(g), "The objective function generated non-finite outputs");

        // gap_eps determines how close we have to get to a bound constraint before we
        // start basically dropping it from the optimization and consider it to be an
        // active constraint.
        const double gap_eps = 1e-8;

        double last_alpha = 1;
        while(stop_strategy.should_continue_search(x, f_value, g))
        {
            s = search_strategy.get_next_direction(x, f_value, zero_bounded_variables(gap_eps, g, x, g, x_lower, x_upper));
            s = gap_step_assign_bounded_variables(gap_eps, s, x, g, x_lower, x_upper);

            double alpha = backtracking_line_search_fixed(
                        make_line_search_function(clamp_function(f,x_lower,x_upper), x, s, f_value),
                        f_value,
                        dot(g,s), // compute gradient for the line search
                        last_alpha, 
                        search_strategy.get_wolfe_rho(), 
                        search_strategy.get_max_line_search_iterations());

            // Do a trust region style thing for alpha.  The idea is that if we take a
            // small step then we are likely to take another small step.  So we reuse the
            // alpha from the last iteration unless the line search didn't shrink alpha at
            // all, in that case, we start with a bigger alpha next time.
            if (alpha == last_alpha)
                last_alpha = std::min(last_alpha*10,1.0);
            else
                last_alpha = alpha;

            // Take the search step indicated by the above line search
            x = clamp(x + alpha*s, x_lower, x_upper);
            g = der(x);

            DLIB_ASSERT(is_finite(f_value), "The objective function generated non-finite outputs");
//			for(int i = 0; i < (int)g.size(); i++) {
//				cout << g(i) << endl;
//			}
            DLIB_ASSERT(is_finite(g), "The objective function generated non-finite outputs");
            //cout << "stop_strategy.should_continue_search(" << x << ", " << f_value << ", " << g << ")" << endl;
        }

        return f_value;
    }

// ----------------------------------------------------------------------------------------
//Model class for TR fitting Currently Unused
class TR_Model 
{
public:
    typedef ::column_vector column_vector;
    typedef matrix<double> general_matrix;

	/**
	 * Evaluates function value with parameters
	 *
	 * @param x parameters that are evaluated
	 * @return function value 
	 */
	double operator() (const column_vector& x) const {
		return TR_funct(data_samples, x);
	}
	void get_derivative_and_hessian (const column_vector& x, column_vector& d, general_matrix& h) const {	
		h = TR_hessian(x);
		d = TR_derivative(x);
	}
};
// ----------------------------------------------------------------------------------------
dlib::rand rand_generator; 
//Add Rician noise to data samples
/**
 * Adds Rician noise to data
 * @param data_samples SI curve
 * @param SD Standard deviation of noise component n1 and n2
 * @param rand_generator pseudo-random number generator
 */
void AddRicianNoise(std::vector<double>& data_samples, double SD)
{
	double SI = 0;
	double r = 0;
	double X = 0;
	double Y = 0;
	for(unsigned long b_i = 0; b_i < data_samples.size(); b_i++) {
		SI = data_samples[b_i];
		//take random number in [0..1] uniform distributed
		r = rand_generator.get_random_gaussian();
		//multiple [0..1] range with Standard Deviation
		X = SI + SI*r*SD;
		Y = SI*r*SD;
		//create sample with rician noise
		data_samples[b_i] = std::sqrt(X*X+Y*Y);
	}
}

// Simulations with model parameters
#define SIMU_IMPLEMENTATION(NAME, NUMBER_OF_PARAMETERS) \
	void NAME(std::vector<double> bvalues, std::vector<std::vector<double> >& SIdata, const std::vector<std::vector<double> >& pinit, const std::vector<std::vector<double> >& simulation_pinit, double (*model_fun)(const double&, const matrix<double,NUMBER_OF_PARAMETERS,1>&), executionInfo& info) { \
		cout << ((int)pinit.size()) << " initializations " << NUMBER_OF_PARAMETERS << " for parameters" << endl; \
		const int no_parameters = NUMBER_OF_PARAMETERS; \
		int no_bvalues = (int)bvalues.size(); \
		cout << no_bvalues << " b-values" << endl; \
		column_vector x(NUMBER_OF_PARAMETERS); \
		SIdata.clear(); \
		clock_t startTime = clock(); \
		unsigned int p_i; \
		std::vector<double> data_samples_no_noise; \
		std::vector<double> data_samples; \
		int iterations = -1; \
		double SD = -1; \
		double C; \
		int PSTEP = (int)pinit.size()/100+1; \
		std::stringstream ss; \
		std::time_t rawtime; \
		std::tm* timeinfo; \
		char buffer [80]; \
		std::time(&rawtime); \
		timeinfo = std::localtime(&rawtime); \
		std::strftime(buffer,80,"%Y-%m-%d-%H-%M-%S",timeinfo); \
		ss << string(buffer); \
		rand_generator.set_seed(ss.str()); \
		for(p_i = 0; p_i < (int)pinit.size(); p_i++) { \
			for(int param_i = 0; param_i < no_parameters-1; param_i++) { \
				x(param_i) = pinit[p_i][param_i]; \
			} \
			C = pinit[p_i][no_parameters-1]; \
			x(no_parameters-1) = C; \
			for(unsigned int simu_p_i = 0; simu_p_i < simulation_pinit.size(); simu_p_i++) { \
				SD = simulation_pinit[simu_p_i][0]; \
				iterations = simulation_pinit[simu_p_i][1]; \
				data_samples_no_noise = std::vector<double>(); \
				try { \
					for(unsigned long b_i = 0; b_i < bvalues.size(); b_i++) { \
						data_samples_no_noise.push_back(model_fun(bvalues[b_i], x)); \
					} \
				} catch(...) { \
					cout << "exception occurred in simulation" << endl; \
					continue; \
				} \
				for(int it = 0; it < iterations; it++) { \
					data_samples = std::vector<double>(); \
					for(unsigned int data_i = 0; data_i < data_samples_no_noise.size(); data_i++) \
						data_samples.push_back(data_samples_no_noise[data_i]); \
					try { \
						AddRicianNoise(data_samples, SD); \
					} catch(...) { \
						cout << "exception occurred in simulation " << endl; \
						continue; \
					} \
					SIdata.push_back(data_samples); \
				} \
			} \
			if(p_i % PSTEP == 0) { \
				cout << "[" << p_i << "]:"; \
				for(int param_i = 0; param_i < no_parameters; param_i++) { \
					cout << info.parameterNames[param_i] << "=" << x(param_i) << std::setprecision(OUTPUT_PRECISION) << " " << std::fixed; \
				} \
				cout << " SD=" << SD << " it=" << iterations << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
				cout << endl; \
			} \
		} \
		cout << "[" << ((int)pinit.size())-1 << "](last):"; \
		for(int param_i = 0; param_i < no_parameters; param_i++) { \
			cout << info.parameterNames[param_i] << "=" << x(param_i) << std::setprecision(OUTPUT_PRECISION) << " " << std::fixed; \
		} \
		cout << " SD=" << SD << " it=" << iterations << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
		cout << endl; \
		info.executiontime = double( clock() - startTime ) / (double)CLOCKS_PER_SEC; \
		cout << "Execution finished in " << info.executiontime << " seconds" << endl; \
	} \

// ----------------------------------------------------------------------------------------
/**
 * Mean Squared Error of Simple Linear Regression
 * @param data_samples SI curve
 * @param slope slope of regression line  
 * @param intercept intercept of regression line
 * @return Mean Squared Error of line to data points
 */
double SimpleLinearRegression_MSE(const matrix<std::pair<double, double>, 0, 1>& data_samples, double slope, double intercept)
{
	double MSE = 0.0;
	double err;
	for(long i = 0; i < (int)data_samples.size(); i++) {
		err = (data_samples(i).first*slope+intercept)-data_samples(i).second;
		MSE = MSE + err*err;
	}
	if((int)data_samples.size() == 0) return -1;
	MSE = MSE/((double)data_samples.size());
	return MSE;
}

/**
 * Simple Linear Regression of Data
 * @param data_samples SI curve
 * @param slope slope of regression line  
 * @param intercept intercept of regression line
 */
void SimpleLinearRegression(const matrix<std::pair<double, double>, 0, 1>& data_samples, double& slope, double& intercept)
{
	double x_mean = 0;
	double y_mean = 0;
	double B_Cov = 0;
	double B_Var = 0;
	double x_dev = 0;
	for(long b_i = 0; b_i < data_samples.size(); b_i++) {
		x_mean += data_samples(b_i).first;
		y_mean += data_samples(b_i).second;
	}
	x_mean /= (double)data_samples.size();
	y_mean /= (double)data_samples.size();	
	for(long b_i = 0; b_i < data_samples.size(); b_i++) {
		x_dev = (data_samples(b_i).first - x_mean);
		B_Cov += x_dev*(data_samples(b_i).second - y_mean);
		B_Var += x_dev*x_dev;
	}
	if(B_Var <= 0) {
		slope = std::numeric_limits<double>::quiet_NaN();
		intercept = std::numeric_limits<double>::quiet_NaN();
		return;
	}
	slope = B_Cov/B_Var;
	intercept = y_mean - slope*x_mean;
}

//Log-linearized version of fitting
#define FIT_CONST_MI_LIN_IMPLEMENTATION(NAME, NUMBER_OF_PARAMETERS, NUMBER_OF_BVALUES) \
	void NAME(std::vector<double> bvalues, const std::vector<std::vector<double> >& SIdata, const std::vector<std::vector<double> >& pinit, double (*residual_MSE_fun)(const matrix<std::pair<double, double>, 0, 1>&, const matrix<double,0,1>&), std::vector<std::vector<double> >& p_vals, executionInfo& info, std::vector<double>& C_values, const column_vector& x_lo, const column_vector& x_hi) { \
		if((int)SIdata.size() == 0) { \
			throw "SIdata length was 0"; \
		} \
		data_samples = matrix<std::pair<double, double>, 0, 1>(NUMBER_OF_BVALUES); \
		double B_0; \
		for(int b_i = 0; b_i < NUMBER_OF_BVALUES; b_i++) { \
			data_samples(b_i).first = bvalues[b_i]; \
		} \
		cout << ((int)pinit.size()) << " initializations " << NUMBER_OF_PARAMETERS << " for parameters" << endl; \
		const int no_parameters = NUMBER_OF_PARAMETERS; \
		int no_bvalues = (int)SIdata[0].size(); \
		cout << no_bvalues << " b-values" << endl; \
		double RMSE; \
		double RMSE_min; \
		double MSE_error; \
		TR_funct = residual_MSE_fun; \
		column_vector x(NUMBER_OF_PARAMETERS); \
		column_vector x_min(NUMBER_OF_PARAMETERS); \
		p_vals.clear(); \
		clock_t startTime = clock(); \
		unsigned int p_i; \
		unsigned int c_i; \
		double C; \
		double SImax; \
		double SImin; \
		double SIrange; \
		for(unsigned int si_i = 0; si_i < SIdata.size(); si_i++) { \
			RMSE_min = std::numeric_limits<double>::max(); \
			std::fill(x_min.begin(), x_min.end(), std::numeric_limits<double>::quiet_NaN()); \
			SImax = std::numeric_limits<double>::min(); \
			SImin = std::numeric_limits<double>::max(); \
			B_0 = SIdata[si_i][0]; \
			if(B_0 == 0) { \
				B_0 = 1; \
			} \
			for(int b_i = 0; b_i < no_bvalues; b_i++) { \
				data_samples(b_i).second = std::log(SIdata[si_i][b_i]/B_0); \
				if(data_samples(b_i).second < SImin) SImin = data_samples(b_i).second; \
				if(data_samples(b_i).second > SImax) SImax = data_samples(b_i).second; \
			} \
			SIrange = SImax - SImin; \
			for(c_i = 0; c_i < (int)C_values.size(); c_i++) { \
				C = SImin+C_values[c_i]*SIrange; \
				for(p_i = 0; p_i < (int)pinit.size(); p_i++) { \
					for(int param_i = 0; param_i < no_parameters-1; param_i++) { \
						x(param_i) = pinit[p_i][param_i]; \
					} \
					x(no_parameters-1) = C; \
					try { \
						find_min_box_constrained_fixed(bfgs_search_strategy(), objective_delta_stop_strategy(DELTA_STOP, MAX_NUMBER_OF_ITERATIONS)OPT_VERBOSE, TR_evaluation, derivative(TR_evaluation), x, x_lo, x_hi); \
					} catch(...) { \
						cout << "exception occurred in optimization" << endl; \
						continue; \
					} \
					MSE_error = TR_evaluation(x); \
					RMSE = std::sqrt(MSE_error); \
					if(RMSE < RMSE_min) { \
						RMSE_min = RMSE; \
						x_min  = x; \
					} \
				} \
			} \
			x_min(no_parameters-1) = std::exp(x_min(no_parameters-1)*B_0); \
			std::vector<double> p_val = std::vector<double>(no_parameters+1); \
			if(si_i % PRINTOUT_STEP == 0) { \
				cout << "[" << si_i << "]:"; \
				for(p_i = 0; p_i < no_parameters; p_i++) { \
					cout << info.parameterNames[p_i] << "=" << x_min(p_i) << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
				} \
				cout << " RMSE="<< RMSE_min << std::setprecision(OUTPUT_PRECISION) << std::fixed << endl; \
			} \
			for(p_i = 0; p_i < no_parameters; p_i++) { \
				p_val[p_i] = x_min(p_i); \
			} \
			for(int b_i = 0; b_i < no_bvalues; b_i++) { \
				data_samples(b_i).second = SIdata[si_i][b_i]; \
			} \
			MSE_error = TR_evaluation(x_min); \
			p_val[no_parameters] = std::sqrt(MSE_error); \
			p_vals.push_back(p_val); \
		} \
		info.executiontime = double( clock() - startTime ) / (double)CLOCKS_PER_SEC; \
		cout << "Execution finished in " << info.executiontime << " seconds" << endl; \
	} \

//Linearized version of fitting with simple linear regression, model can have only two parameters: slope and intercept
#define FIT_CONST_MI_SIMPLE_LIN_IMPLEMENTATION(NAME, NUMBER_OF_BVALUES, NUMBER_OF_TAIL_BVALUES) \
	void NAME(std::vector<double> bvalues, const std::vector<std::vector<double> >& SIdata, const std::vector<std::vector<double> >& pinit, double (*residual_MSE_fun)(const matrix<std::pair<double, double>, 0, 1>&, const matrix<double,0,1>&), std::vector<std::vector<double> >& p_vals, executionInfo& info, const column_vector& x_lo, const column_vector& x_hi) { \
		if((int)SIdata.size() == 0) { \
			throw "SIdata length was 0"; \
		} \
		data_samples = matrix<std::pair<double, double>, 0, 1>(NUMBER_OF_BVALUES); \
		double B_0; \
		for(int b_i = 0; b_i < NUMBER_OF_BVALUES; b_i++) { \
			data_samples(b_i).first = bvalues[b_i]; \
		} \
		cout << ((int)pinit.size()) << " initializations for " << 2 << " parameters" << endl; \
		cout << ((int)pinit.size()) << " using " << NUMBER_OF_TAIL_BVALUES << " b-values from tail" << endl; \
		const int no_parameters = 2; \
		int no_bvalues = NUMBER_OF_TAIL_BVALUES; \
		cout << no_bvalues << " b-values" << endl; \
		double RMSE_min; \
		double MSE_error; \
		TR_funct = residual_MSE_fun; \
		column_vector x(2); \
		column_vector x_min(2); \
		p_vals.clear(); \
		clock_t startTime = clock(); \
		unsigned int p_i; \
		double f; \
		double D; \
		for(unsigned int si_i = 0; si_i < SIdata.size(); si_i++) { \
			RMSE_min = std::numeric_limits<double>::max(); \
			std::fill(x_min.begin(), x_min.end(), std::numeric_limits<double>::quiet_NaN()); \
			B_0 = SIdata[si_i][0]; \
			if(B_0 == 0) { \
				B_0 = 1; \
			} \
			for(int b_i = 0; b_i < no_bvalues; b_i++) { \
				data_samples(b_i).second = std::log(SIdata[si_i][b_i]/B_0); \
			} \
			try { \
				SimpleLinearRegression(data_samples, D, f); \
				RMSE_min = SimpleLinearRegression_MSE(data_samples, D, f); \
				x_min(0) = -f; \
				x_min(1) = D; \
			} catch(...) { \
				cout << "exception occurred in optimization" << endl; \
				continue; \
			} \
			std::vector<double> p_val = std::vector<double>(no_parameters+1); \
			if(si_i % PRINTOUT_STEP == 0) { \
				cout << "[" << si_i << "]:"; \
				for(p_i = 0; p_i < no_parameters; p_i++) { \
					cout << info.parameterNames[p_i] << "=" << x_min(p_i) << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
				} \
				cout << " RMSE="<< RMSE_min << std::setprecision(OUTPUT_PRECISION) << std::fixed << endl; \
			} \
			for(p_i = 0; p_i < no_parameters; p_i++) { \
				p_val[p_i] = x_min(p_i); \
			} \
			for(int b_i = 0; b_i < no_bvalues; b_i++) { \
				data_samples(b_i).second = SIdata[si_i][b_i]; \
			} \
			MSE_error = TR_evaluation(x_min); \
			p_val[no_parameters] = std::sqrt(MSE_error); \
			p_vals.push_back(p_val); \
		} \
		info.executiontime = double( clock() - startTime ) / (double)CLOCKS_PER_SEC; \
		cout << "Execution finished in " << info.executiontime << " seconds" << endl; \
	} \

#define FIT_CONST_MI_IMPLEMENTATION(NAME, NUMBER_OF_PARAMETERS, NUMBER_OF_BVALUES) \
	void NAME(std::vector<double> bvalues, const std::vector<std::vector<double> >& SIdata, const std::vector<std::vector<double> >& pinit, double (*residual_MSE_fun)(const matrix<std::pair<double, double>, 0, 1>&, const matrix<double,0,1>&), std::vector<std::vector<double> >& p_vals, executionInfo& info, std::vector<double>& C_values, const column_vector& x_lo, const column_vector& x_hi) { \
		if((int)SIdata.size() == 0) { \
			throw "SIdata length was 0"; \
		} \
		data_samples = matrix<std::pair<double, double>, 0, 1>(NUMBER_OF_BVALUES); \
		double B_0; \
		for(int b_i = 0; b_i < NUMBER_OF_BVALUES; b_i++) { \
			data_samples(b_i).first = bvalues[b_i]; \
		} \
		cout << ((int)pinit.size()) << " initializations " << NUMBER_OF_PARAMETERS << " for parameters" << endl; \
		const int no_parameters = NUMBER_OF_PARAMETERS; \
		int no_bvalues = (int)SIdata[0].size(); \
		cout << no_bvalues << " b-values" << endl; \
		double RMSE; \
		double RMSE_min; \
		double MSE_error; \
		TR_funct = residual_MSE_fun; \
		column_vector x(NUMBER_OF_PARAMETERS); \
		column_vector x_min(NUMBER_OF_PARAMETERS); \
		p_vals.clear(); \
		clock_t startTime = clock(); \
		unsigned int p_i; \
		unsigned int c_i; \
		double C; \
		double SImax; \
		double SImin; \
		double SIrange; \
		for(unsigned int si_i = 0; si_i < SIdata.size(); si_i++) { \
			RMSE_min = std::numeric_limits<double>::max(); \
			std::fill(x_min.begin(), x_min.end(), std::numeric_limits<double>::quiet_NaN()); \
			SImax = std::numeric_limits<double>::min(); \
			SImin = std::numeric_limits<double>::max(); \
			B_0 = SIdata[si_i][0]; \
			if(B_0 == 0) { \
				B_0 = 1; \
			} \
			for(int b_i = 0; b_i < no_bvalues; b_i++) { \
				data_samples(b_i).second = SIdata[si_i][b_i]/B_0; \
				if(data_samples(b_i).second < SImin) SImin = data_samples(b_i).second; \
				if(data_samples(b_i).second > SImax) SImax = data_samples(b_i).second; \
			} \
			SIrange = SImax - SImin; \
			for(c_i = 0; c_i < (int)C_values.size(); c_i++) { \
				C = SImin+C_values[c_i]*SIrange; \
				for(p_i = 0; p_i < (int)pinit.size(); p_i++) { \
					for(int param_i = 0; param_i < no_parameters-1; param_i++) { \
						x(param_i) = pinit[p_i][param_i]; \
					} \
					x(no_parameters-1) = C; \
					try { \
						find_min_box_constrained_fixed(bfgs_search_strategy(), objective_delta_stop_strategy(DELTA_STOP, MAX_NUMBER_OF_ITERATIONS)OPT_VERBOSE, TR_evaluation, derivative(TR_evaluation), x, x_lo, x_hi); \
					} catch(...) { \
						cout << "exception occurred in optimization" << endl; \
						continue; \
					} \
					MSE_error = TR_evaluation(x); \
					RMSE = std::sqrt(MSE_error); \
					if(RMSE < RMSE_min) { \
						RMSE_min = RMSE; \
						x_min  = x; \
					} \
				} \
			} \
			x_min(no_parameters-1) = x_min(no_parameters-1)*B_0; \
			std::vector<double> p_val = std::vector<double>(no_parameters+1); \
			if(si_i % PRINTOUT_STEP == 0) { \
				cout << "[" << si_i << "]:"; \
				for(p_i = 0; p_i < no_parameters; p_i++) { \
					cout << info.parameterNames[p_i] << "=" << x_min(p_i) << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
				} \
				cout << " RMSE="<< RMSE_min << std::setprecision(OUTPUT_PRECISION) << std::fixed << endl; \
			} \
			for(p_i = 0; p_i < no_parameters; p_i++) { \
				p_val[p_i] = x_min(p_i); \
			} \
			p_val[no_parameters] = RMSE_min*B_0; \
			p_vals.push_back(p_val); \
		} \
		info.executiontime = double( clock() - startTime ) / (double)CLOCKS_PER_SEC; \
		cout << "Execution finished in " << info.executiontime << " seconds" << endl; \
	} \

//Normalized fit
#define FIT_CONST_MI_NORMALIZED_IMPLEMENTATION(NAME, NUMBER_OF_PARAMETERS, NUMBER_OF_BVALUES) \
	void NAME(std::vector<double> bvalues, const std::vector<std::vector<double> >& SIdata, const std::vector<std::vector<double> >& pinit, double (*residual_MSE_fun)(const matrix<std::pair<double, double>, 0, 1>&, const matrix<double,0,1>&), std::vector<std::vector<double> >& p_vals, executionInfo& info, const column_vector& x_lo, const column_vector& x_hi) { \
		if((int)SIdata.size() == 0) { \
			throw "SIdata length was 0"; \
		} \
		data_samples = matrix<std::pair<double, double>, 0, 1>(NUMBER_OF_BVALUES); \
		double B_0; \
		for(int b_i = 0; b_i < NUMBER_OF_BVALUES; b_i++) { \
			data_samples(b_i).first = bvalues[b_i]; \
		} \
		cout << ((int)pinit.size()) << " initializations " << NUMBER_OF_PARAMETERS << " for parameters" << endl; \
		const int no_parameters = NUMBER_OF_PARAMETERS; \
		int no_bvalues = (int)SIdata[0].size(); \
		cout << no_bvalues << " b-values" << endl; \
		double RMSE; \
		double RMSE_min; \
		double MSE_error; \
		TR_funct = residual_MSE_fun; \
		column_vector x(NUMBER_OF_PARAMETERS); \
		column_vector x_min(NUMBER_OF_PARAMETERS); \
		p_vals.clear(); \
		clock_t startTime = clock(); \
		unsigned int p_i; \
		for(unsigned int si_i = 0; si_i < SIdata.size(); si_i++) { \
			RMSE_min = std::numeric_limits<double>::max(); \
			std::fill(x_min.begin(), x_min.end(), std::numeric_limits<double>::quiet_NaN()); \
			B_0 = SIdata[si_i][0]; \
			if(B_0 == 0) { \
				B_0 = 1; \
			} \
			for(int b_i = 0; b_i < no_bvalues; b_i++) { \
				data_samples(b_i).second = SIdata[si_i][b_i]/B_0; \
			} \
			for(p_i = 0; p_i < (int)pinit.size(); p_i++) { \
				for(int param_i = 0; param_i < no_parameters; param_i++) { \
					x(param_i) = pinit[p_i][param_i]; \
				} \
				try { \
					find_min_box_constrained_fixed(bfgs_search_strategy(), objective_delta_stop_strategy(DELTA_STOP, MAX_NUMBER_OF_ITERATIONS)OPT_VERBOSE, TR_evaluation, derivative(TR_evaluation), x, x_lo, x_hi); \
				} catch(...) { \
					cout << "exception occurred in optimization" << endl; \
					continue; \
				} \
				MSE_error = TR_evaluation(x); \
				RMSE = std::sqrt(MSE_error); \
				if(RMSE < RMSE_min) { \
					RMSE_min = RMSE; \
					x_min  = x; \
				} \
			} \
			std::vector<double> p_val = std::vector<double>(no_parameters+1); \
			if(si_i % PRINTOUT_STEP == 0) { \
				cout << "[" << si_i << "]:"; \
				for(p_i = 0; p_i < no_parameters; p_i++) { \
					cout << info.parameterNames[p_i] << "=" << x_min(p_i) << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
				} \
				cout << " RMSE="<< RMSE_min << std::setprecision(OUTPUT_PRECISION) << std::fixed << endl; \
			} \
			for(p_i = 0; p_i < no_parameters; p_i++) { \
				p_val[p_i] = x_min(p_i); \
			} \
			p_val[no_parameters] = RMSE_min*B_0; \
			p_vals.push_back(p_val); \
		} \
		info.executiontime = double( clock() - startTime ) / (double)CLOCKS_PER_SEC; \
		cout << "Execution finished in " << info.executiontime << " seconds" << endl; \
	} \

#define FIT_CONST_MI_NORMALIZED_IMPLEMENTATION_BIEXP(NAME, NUMBER_OF_PARAMETERS, NUMBER_OF_BVALUES) \
	void NAME(std::vector<double> bvalues, const std::vector<std::vector<double> >& SIdata, const std::vector<std::vector<double> >& pinit, double (*residual_MSE_fun)(const matrix<std::pair<double, double>, 0, 1>&, const matrix<double,0,1>&), std::vector<std::vector<double> >& p_vals, executionInfo& info, const column_vector& x_lo, const column_vector& x_hi) { \
		if((int)SIdata.size() == 0) { \
			throw "SIdata length was 0"; \
		} \
		data_samples = matrix<std::pair<double, double>, 0, 1>(NUMBER_OF_BVALUES); \
		double B_0; \
		for(int b_i = 0; b_i < NUMBER_OF_BVALUES; b_i++) { \
			data_samples(b_i).first = bvalues[b_i]; \
		} \
		cout << ((int)pinit.size()) << " initializations for "<< NUMBER_OF_PARAMETERS << " parameters" << endl; \
		const int no_parameters = NUMBER_OF_PARAMETERS; \
		int no_bvalues = (int)SIdata[0].size(); \
		cout << no_bvalues << " b-values" << endl; \
		double RMSE; \
		double RMSE_min; \
		double MSE_error; \
		TR_funct = residual_MSE_fun; \
		column_vector x(NUMBER_OF_PARAMETERS); \
		column_vector x_min(NUMBER_OF_PARAMETERS); \
		p_vals.clear(); \
		clock_t startTime = clock(); \
		unsigned int p_i; \
		double x_swap; \
		for(unsigned int si_i = 0; si_i < (int)SIdata.size(); si_i++) { \
			RMSE_min = std::numeric_limits<double>::max(); \
			std::fill(x_min.begin(), x_min.end(), std::numeric_limits<double>::quiet_NaN()); \
			B_0 = SIdata[si_i][0]; \
			if(B_0 == 0) { \
				B_0 = 1; \
			} \
			for(int b_i = 0; b_i < no_bvalues; b_i++) { \
				data_samples(b_i).second = SIdata[si_i][b_i]/B_0; \
			} \
			for(p_i = 0; p_i < (int)pinit.size(); p_i++) { \
				for(int param_i = 0; param_i < no_parameters; param_i++) { \
					x(param_i) = pinit[p_i][param_i]; \
				} \
				try { \
					find_min_box_constrained_fixed(bfgs_search_strategy(), objective_delta_stop_strategy(DELTA_STOP, MAX_NUMBER_OF_ITERATIONS)OPT_VERBOSE, TR_evaluation, derivative(TR_evaluation), x, x_lo, x_hi); \
				} catch(...) { \
					cout << "exception occurred in optimization" << endl; \
					continue; \
				} \
				if(x(1) < x(2)) { \
					x_swap = x(1); \
					x(1) = x(2); \
					x(2) = x_swap; \
					x(0) = 1-x(0); \
				} \
				MSE_error = TR_evaluation(x); \
				RMSE = std::sqrt(MSE_error); \
				if(RMSE < RMSE_min) { \
					RMSE_min = RMSE; \
					x_min  = x; \
				} \
			} \
			std::vector<double> p_val = std::vector<double>(no_parameters+1); \
			if(si_i % PRINTOUT_STEP == 0) { \
				cout << "[" << si_i << "]:"; \
				for(p_i = 0; p_i < no_parameters; p_i++) { \
					cout << info.parameterNames[p_i] << "=" << x_min(p_i) << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
				} \
				cout << " RMSE="<< RMSE_min << std::setprecision(OUTPUT_PRECISION) << std::fixed << endl; \
			} \
			for(p_i = 0; p_i < no_parameters; p_i++) { \
				p_val[p_i] = x_min(p_i); \
			} \
			p_val[no_parameters] = RMSE_min*B_0; \
			p_vals.push_back(p_val); \
		} \
		info.executiontime = double( clock() - startTime ) / (double)CLOCKS_PER_SEC; \
		cout << "Execution finished in " << info.executiontime << " seconds" << endl; \
	} \

#define FIT_CONST_MI_IMPLEMENTATION_BIEXP(NAME, NUMBER_OF_PARAMETERS, NUMBER_OF_BVALUES) \
	void NAME(std::vector<double> bvalues, const std::vector<std::vector<double> >& SIdata, const std::vector<std::vector<double> >& pinit, double (*residual_MSE_fun)(const matrix<std::pair<double, double>, 0, 1>&, const matrix<double,0,1>&), std::vector<std::vector<double> >& p_vals, executionInfo& info, std::vector<double>& C_values, const column_vector& x_lo, const column_vector& x_hi) { \
		if((int)SIdata.size() == 0) { \
			throw "SIdata length was 0"; \
		} \
		data_samples = matrix<std::pair<double, double>, 0, 1>(NUMBER_OF_BVALUES); \
		double B_0; \
		for(int b_i = 0; b_i < NUMBER_OF_BVALUES; b_i++) { \
			data_samples(b_i).first = bvalues[b_i]; \
		} \
		cout << ((int)pinit.size()) << " initializations for "<< NUMBER_OF_PARAMETERS << " parameters" << endl; \
		const int no_parameters = NUMBER_OF_PARAMETERS; \
		int no_bvalues = (int)SIdata[0].size(); \
		cout << no_bvalues << " b-values" << endl; \
		double RMSE; \
		double RMSE_min; \
		double MSE_error; \
		TR_funct = residual_MSE_fun; \
		column_vector x(NUMBER_OF_PARAMETERS); \
		column_vector x_min(NUMBER_OF_PARAMETERS); \
		p_vals.clear(); \
		clock_t startTime = clock(); \
		unsigned int p_i; \
		unsigned int c_i; \
		double C; \
		double SImax; \
		double SImin; \
		double SIrange; \
		double x_swap; \
		for(unsigned int si_i = 0; si_i < (int)SIdata.size(); si_i++) { \
			RMSE_min = std::numeric_limits<double>::max(); \
			std::fill(x_min.begin(), x_min.end(), std::numeric_limits<double>::quiet_NaN()); \
			SImax = std::numeric_limits<double>::min(); \
			SImin = std::numeric_limits<double>::max(); \
			B_0 = SIdata[si_i][0]; \
			if(B_0 == 0) { \
				B_0 = 1; \
			} \
			for(int b_i = 0; b_i < no_bvalues; b_i++) { \
				data_samples(b_i).second = SIdata[si_i][b_i]/B_0; \
				if(data_samples(b_i).second < SImin) SImin = data_samples(b_i).second; \
				if(data_samples(b_i).second > SImax) SImax = data_samples(b_i).second; \
			} \
			SIrange = SImax - SImin; \
			for(c_i = 0; c_i < (int)C_values.size(); c_i++) { \
				C = SImin+C_values[c_i]*SIrange; \
				for(p_i = 0; p_i < (int)pinit.size(); p_i++) { \
					for(int param_i = 0; param_i < no_parameters; param_i++) { \
						x(param_i) = pinit[p_i][param_i]; \
					} \
					x(no_parameters-1) = C; \
					try { \
						find_min_box_constrained_fixed(bfgs_search_strategy(), objective_delta_stop_strategy(DELTA_STOP, MAX_NUMBER_OF_ITERATIONS)OPT_VERBOSE, TR_evaluation, derivative(TR_evaluation), x, x_lo, x_hi); \
					} catch(...) { \
						cout << "exception occurred in optimization" << endl; \
						continue; \
					} \
					if(x(1) < x(2)) { \
						x_swap = x(1); \
						x(1) = x(2); \
						x(2) = x_swap; \
						x(0) = 1-x(0); \
					} \
					MSE_error = TR_evaluation(x); \
					RMSE = std::sqrt(MSE_error); \
					if(RMSE < RMSE_min) { \
						RMSE_min = RMSE; \
						x_min  = x; \
					} \
				} \
			} \
			x_min(no_parameters-1) = x_min(no_parameters-1)*B_0; \
			std::vector<double> p_val = std::vector<double>(no_parameters+1); \
			if(si_i % PRINTOUT_STEP == 0) { \
				cout << "[" << si_i << "]:"; \
				for(p_i = 0; p_i < no_parameters; p_i++) { \
					cout << info.parameterNames[p_i] << "=" << x_min(p_i) << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
				} \
				cout << " RMSE="<< RMSE_min << std::setprecision(OUTPUT_PRECISION) << std::fixed << endl; \
			} \
			for(p_i = 0; p_i < no_parameters; p_i++) { \
				p_val[p_i] = x_min(p_i); \
			} \
			p_val[no_parameters] = RMSE_min*B_0; \
			p_vals.push_back(p_val); \
		} \
		info.executiontime = double( clock() - startTime ) / (double)CLOCKS_PER_SEC; \
		cout << "Execution finished in " << info.executiontime << " seconds" << endl; \
	} \

//Macro to define LM fit with fixed number of parameter values, implementation
#define FIT_CONST_MI_IMPLEMENTATION_MARZISEGMENTEDBIEXP(NAME, NUMBER_OF_BVALUES, NUMBER_OF_TAIL_BVALUES) \
	void NAME(std::vector<double> bvalues, const std::vector<std::vector<double> >& SIdata, const std::vector<std::vector<double> >& pinit, const std::vector<std::vector<double> >& pinit2, double (*residual_MSE_fun)(const matrix<std::pair<double, double>, 0, 1>&, const matrix<double, 0, 1>&), double (*residual_MSE_fun_2nd)(const matrix<std::pair<double, double>, 0, 1>&, const matrix<double, 0, 1>&), std::vector<std::vector<double> >& p_vals, executionInfo& info, std::vector<double>& C_values, const column_vector& x_lo1, const column_vector& x_hi1, const column_vector& x_lo2, const column_vector& x_hi2) { \
		if((int)SIdata.size() == 0) { \
			throw "SIdata length was 0"; \
		} \
		data_samples = matrix<std::pair<double, double>, 0, 1>(NUMBER_OF_TAIL_BVALUES); \
		data_samples2 = matrix<std::pair<double, double>, 0, 1>(NUMBER_OF_BVALUES); \
		double B_0; \
		for(int b_i = 0; b_i < NUMBER_OF_BVALUES; b_i++) { \
			data_samples2(b_i).first = bvalues[b_i]; \
		} \
		for(int b_i = (NUMBER_OF_BVALUES-NUMBER_OF_TAIL_BVALUES); b_i < NUMBER_OF_BVALUES; b_i++) { \
			data_samples(b_i-(NUMBER_OF_BVALUES-NUMBER_OF_TAIL_BVALUES)).first = bvalues[b_i]; \
		} \
		cout << ((int)pinit.size()) << " initializations for 1+3 parameters" << endl; \
		const int no_parameters = 4; \
		int no_bvalues = (int)SIdata[0].size(); \
		cout << no_bvalues << " b-values" << endl; \
		double RMSE; \
		double RMSE_min; \
		double RMSE_min_1st; \
		double MSE_error; \
		TR_funct = residual_MSE_fun; \
		TR_funct2 = residual_MSE_fun_2nd; \
		matrix<double,2,1> x_1st; \
		matrix<double,2,1> x_min_1st; \
		matrix<double,3,1> x_2nd; \
		matrix<double,4,1> x; \
		matrix<double,4,1> x_min; \
		p_vals.clear(); \
		clock_t startTime = clock(); \
		unsigned int p_i; \
		fixed_parameters.clear(); \
		fixed_parameters.push_back(0); \
		unsigned int c_i; \
		double C; \
		double SImax; \
		double SImin; \
		double SIrange; \
		for(unsigned int si_i = 0; si_i < (int)SIdata.size(); si_i++) { \
			RMSE_min = std::numeric_limits<double>::infinity(); \
			SImax = std::numeric_limits<double>::min(); \
			SImin = std::numeric_limits<double>::max(); \
			B_0 = SIdata[si_i][0]; \
			if(B_0 == 0) { \
				B_0 = 1; \
			} \
			for(int b_i = no_bvalues-NUMBER_OF_TAIL_BVALUES; b_i < no_bvalues; b_i++) { \
				data_samples(b_i-(no_bvalues-NUMBER_OF_TAIL_BVALUES)).second = SIdata[si_i][b_i]/B_0; \
				if(data_samples(b_i-(no_bvalues-NUMBER_OF_TAIL_BVALUES)).second < SImin) SImin = data_samples(b_i-(no_bvalues-NUMBER_OF_TAIL_BVALUES)).second; \
				if(data_samples(b_i-(no_bvalues-NUMBER_OF_TAIL_BVALUES)).second > SImax) SImax = data_samples(b_i-(no_bvalues-NUMBER_OF_TAIL_BVALUES)).second; \
			} \
			for(int b_i = 0; b_i < no_bvalues; b_i++) { \
				data_samples2(b_i).second = SIdata[si_i][b_i]/B_0; \
			} \
			SIrange = SImax - SImin; \
			std::fill(x_min_1st.begin(), x_min_1st.end(), std::numeric_limits<double>::quiet_NaN()); \
			RMSE_min_1st = std::numeric_limits<double>::infinity(); \
			for(c_i = 0; c_i < (int)C_values.size(); c_i++) { \
				C = SImin+C_values[c_i]*SIrange; \
				for(p_i = 0; p_i < (int)pinit.size(); p_i++) { \
					x_1st(0) = pinit[p_i][0]; \
					x_1st(1) = C; \
					find_min_box_constrained_fixed(bfgs_search_strategy(), objective_delta_stop_strategy(DELTA_STOP, MAX_NUMBER_OF_ITERATIONS)OPT_VERBOSE, TR_evaluation, derivative(TR_evaluation), x_1st, x_lo1, x_hi1); \
					MSE_error = TR_evaluation(x_1st); \
					RMSE = std::sqrt(MSE_error); \
					if(x_1st(0) < 0) x_1st(0) = 0; \
					if(RMSE < RMSE_min_1st) { \
						RMSE_min_1st = RMSE; \
						x_min_1st  = x_1st; \
					} \
				} \
			} \
			RMSE_min = std::numeric_limits<double>::infinity(); \
			std::fill(x_min.begin(), x_min.end(), std::numeric_limits<double>::quiet_NaN()); \
			for(c_i = 0; c_i < (int)C_values.size(); c_i++) { \
				C = SImin+C_values[c_i]*SIrange; \
				for(p_i = 0; p_i < (int)pinit2.size(); p_i++) { \
					x_2nd(0) = pinit2[p_i][0]; \
					x_2nd(1) = pinit2[p_i][1]; \
					fixed_parameters[0] = x_min_1st(0); \
					x_2nd(2) = C; \
					find_min_box_constrained_fixed(bfgs_search_strategy(), objective_delta_stop_strategy(DELTA_STOP, MAX_NUMBER_OF_ITERATIONS)OPT_VERBOSE, TR_evaluation2, derivative(TR_evaluation2), x_2nd, x_lo2, x_hi2); \
					x(0) = x_2nd(0); \
					x(1) = x_2nd(1); \
					x(2) = x_min_1st(0); \
					x(3) = x_2nd(2); \
					MSE_error = TR_evaluation2(x_2nd); \
					RMSE = std::sqrt(MSE_error); \
					if(RMSE < RMSE_min) { \
						RMSE_min = RMSE; \
						x_min  = x; \
					} \
				} \
			} \
			std::vector<double> p_val = std::vector<double>(no_parameters+1); \
			x_min(no_parameters-1) = x_min(no_parameters-1)*B_0; \
			if(si_i % PRINTOUT_STEP == 0) { \
				cout << "[" << si_i << "]:"; \
				for(p_i = 0; p_i < no_parameters; p_i++) { \
					cout << info.parameterNames[p_i] << "=" << x_min(p_i) << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
				} \
				cout << " RMSE="<< RMSE_min << std::setprecision(OUTPUT_PRECISION) << std::fixed << endl; \
			} \
			for(p_i = 0; p_i < no_parameters; p_i++) { \
				p_val[p_i] = x_min(p_i); \
			} \
			p_val[no_parameters] = RMSE_min*B_0; \
			p_vals.push_back(p_val); \
		} \
		info.executiontime = double( clock() - startTime ) / (double)CLOCKS_PER_SEC; \
		cout << "Execution finished in " << info.executiontime << " seconds" << endl; \
	} \

#define FIT_CONST_MI_IMPLEMENTATION_CHOSEGMENTEDBIEXP(NAME, NUMBER_OF_BVALUES, NUMBER_OF_TAIL_BVALUES) \
	void NAME(std::vector<double> bvalues, const std::vector<std::vector<double> >& SIdata, const std::vector<std::vector<double> >& pinit, const std::vector<std::vector<double> >& pinit2, double (*residual_MSE_fun)(const matrix<std::pair<double, double>, 0, 1>&, const matrix<double, 0, 1>&), double (*residual_MSE_fun_2nd)(const matrix<std::pair<double, double>, 0, 1>&, const matrix<double, 0, 1>&), std::vector<std::vector<double> >& p_vals, executionInfo& info, std::vector<double>& C_values, const column_vector& x_lo1, const column_vector& x_hi1, const column_vector& x_lo2, const column_vector& x_hi2) { \
		if((int)SIdata.size() == 0) { \
			throw "SIdata length was 0"; \
		} \
		data_samples = matrix<std::pair<double, double>, 0, 1>(NUMBER_OF_TAIL_BVALUES); \
		data_samples2 = matrix<std::pair<double, double>, 0, 1>(NUMBER_OF_BVALUES); \
		double B_0; \
		for(int b_i = 0; b_i < NUMBER_OF_BVALUES; b_i++) { \
			data_samples2(b_i).first = bvalues[b_i]; \
		} \
		for(int b_i = (NUMBER_OF_BVALUES-NUMBER_OF_TAIL_BVALUES); b_i < NUMBER_OF_BVALUES; b_i++) { \
			data_samples(b_i-(NUMBER_OF_BVALUES-NUMBER_OF_TAIL_BVALUES)).first = bvalues[b_i]; \
		} \
		cout << ((int)pinit.size()) << " initializations for 1+3 parameters" << endl; \
		const int no_parameters = 4; \
		int no_bvalues = (int)SIdata[0].size(); \
		cout << no_bvalues << " b-values" << endl; \
		double RMSE; \
		double RMSE_min; \
		double RMSE_min_1st; \
		double MSE_error; \
		TR_funct = residual_MSE_fun; \
		TR_funct2 = residual_MSE_fun_2nd; \
		matrix<double,2,1> x_1st; \
		matrix<double,2,1> x_min_1st; \
		matrix<double,2,1> x_2nd; \
		matrix<double,4,1> x; \
		matrix<double,4,1> x_min; \
		p_vals.clear(); \
		clock_t startTime = clock(); \
		unsigned int p_i; \
		fixed_parameters.clear(); \
		fixed_parameters.push_back(0); \
		fixed_parameters.push_back(0); \
		unsigned int c_i; \
		double C; \
		double SImax; \
		double SImin; \
		double SIrange; \
		for(unsigned int si_i = 0; si_i < (int)SIdata.size(); si_i++) { \
			RMSE_min = std::numeric_limits<double>::infinity(); \
			SImax = std::numeric_limits<double>::min(); \
			SImin = std::numeric_limits<double>::max(); \
			B_0 = SIdata[si_i][0]; \
			if(B_0 == 0) { \
				B_0 = 1; \
			} \
			for(int b_i = no_bvalues-NUMBER_OF_TAIL_BVALUES; b_i < no_bvalues; b_i++) { \
				data_samples(b_i-(no_bvalues-NUMBER_OF_TAIL_BVALUES)).second = SIdata[si_i][b_i]/B_0; \
				if(data_samples(b_i-(no_bvalues-NUMBER_OF_TAIL_BVALUES)).second < SImin) SImin = data_samples(b_i-(no_bvalues-NUMBER_OF_TAIL_BVALUES)).second; \
				if(data_samples(b_i-(no_bvalues-NUMBER_OF_TAIL_BVALUES)).second > SImax) SImax = data_samples(b_i-(no_bvalues-NUMBER_OF_TAIL_BVALUES)).second; \
			} \
			for(int b_i = 0; b_i < no_bvalues; b_i++) { \
				data_samples2(b_i).second = SIdata[si_i][b_i]/B_0; \
			} \
			SIrange = SImax - SImin; \
			std::fill(x_min_1st.begin(), x_min_1st.end(), std::numeric_limits<double>::quiet_NaN()); \
			RMSE_min_1st = std::numeric_limits<double>::infinity(); \
			for(c_i = 0; c_i < (int)C_values.size(); c_i++) { \
				C = SImin+C_values[c_i]*SIrange; \
				for(p_i = 0; p_i < (int)pinit.size(); p_i++) { \
					x_1st(0) = pinit[p_i][0]; \
					x_1st(1) = C; \
					find_min_box_constrained_fixed(bfgs_search_strategy(), objective_delta_stop_strategy(DELTA_STOP, MAX_NUMBER_OF_ITERATIONS)OPT_VERBOSE, TR_evaluation, derivative(TR_evaluation), x_1st, x_lo1, x_hi1); \
					MSE_error = TR_evaluation(x_1st); \
					RMSE = std::sqrt(MSE_error); \
					if(x_1st(0) < 0) x_1st(0) = 0; \
					if(RMSE < RMSE_min_1st) { \
						RMSE_min_1st = RMSE; \
						x_min_1st  = x_1st; \
					} \
				} \
			} \
			RMSE_min = std::numeric_limits<double>::infinity(); \
			std::fill(x_min.begin(), x_min.end(), std::numeric_limits<double>::quiet_NaN()); \
			for(c_i = 0; c_i < (int)C_values.size(); c_i++) { \
				C = SImin+C_values[c_i]*SIrange; \
				for(p_i = 0; p_i < (int)pinit2.size(); p_i++) { \
					fixed_parameters[0] = (1.0-x_min_1st(1)); \
					if(fixed_parameters[0] < 0.001) { \
						fixed_parameters[0] = 0.001; \
					} \
					x_2nd(0) = pinit2[p_i][0]; \
					fixed_parameters[1] = x_min_1st(0); \
					x_2nd(1) = C; \
					find_min_box_constrained_fixed(bfgs_search_strategy(), objective_delta_stop_strategy(DELTA_STOP, MAX_NUMBER_OF_ITERATIONS)OPT_VERBOSE, TR_evaluation2, derivative(TR_evaluation2), x_2nd, x_lo2, x_hi2); \
					x(0) = (1.0-x_min_1st(1)); \
					x(1) = x_2nd(0); \
					x(2) = x_min_1st(0); \
					x(3) = x_2nd(1); \
					MSE_error = TR_evaluation2(x_2nd); \
					RMSE = std::sqrt(MSE_error); \
					if(RMSE < RMSE_min) { \
						RMSE_min = RMSE; \
						x_min  = x; \
					} \
				} \
			} \
			std::vector<double> p_val = std::vector<double>(no_parameters+1); \
			x_min(no_parameters-1) = x_min(no_parameters-1)*B_0; \
			if(si_i % PRINTOUT_STEP == 0) { \
				cout << "[" << si_i << "]:"; \
				cout << "f(analytical)=" << x_min(0) << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
				cout << info.parameterNames[0] << "=" << x_min(1) << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
				cout << info.parameterNames[1] << "=" << x_min(2) << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
				cout << info.parameterNames[2] << "=" << x_min(3) << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
				cout << " RMSE="<< RMSE_min << std::setprecision(OUTPUT_PRECISION) << std::fixed << endl; \
			} \
			for(p_i = 0; p_i < no_parameters; p_i++) { \
				p_val[p_i] = x_min(p_i); \
			} \
			p_val[no_parameters] = RMSE_min*B_0; \
			p_vals.push_back(p_val); \
		} \
		info.executiontime = double( clock() - startTime ) / (double)CLOCKS_PER_SEC; \
		cout << "Execution finished in " << info.executiontime << " seconds" << endl; \
	} \

// ----------------------------------------------------------------------------------------
//Asymptotic fitting
// [1] J. Pekar, C.T. Moonen, P.C. van Zijl, 
// On the precision of diffusion/perfusion imaging by gradient sensitization., 
// Magn. Reson. Med. 23 (1992) 122–9.
#define FIT_MI_IMPLEMENTATION_ASYMPTOTICBIEXP(NAME, NUMBER_OF_BVALUES, NUMBER_OF_TAIL_BVALUES) \
	void NAME(std::vector<double> bvalues, const std::vector<std::vector<double> >& SIdata, const std::vector<std::vector<double> >& pinit, double (*residual_MSE_fun)(const matrix<std::pair<double, double>, 0, 1>&, const matrix<double,0,1>&), std::vector<std::vector<double> >& p_vals, executionInfo& info, std::vector<double>& C_values, const column_vector& x_lo, const column_vector& x_hi) { \
		if((int)SIdata.size() == 0) { \
			throw "SIdata length was 0"; \
		} \
		data_samples = matrix<std::pair<double, double>, 0, 1>(NUMBER_OF_TAIL_BVALUES); \
		double B_0; \
		for(int b_i = NUMBER_OF_BVALUES-NUMBER_OF_TAIL_BVALUES; b_i < NUMBER_OF_BVALUES; b_i++) { \
			data_samples(b_i-(NUMBER_OF_BVALUES-NUMBER_OF_TAIL_BVALUES)).first = bvalues[b_i]; \
		} \
		cout << ((int)pinit.size()) << " initializations for 2 parameters (+ 1 derived parameter not initialized) " << endl; \
		const int no_parameters = 2; \
		int no_bvalues = (int)SIdata[0].size(); \
		cout << no_bvalues << " b-values" << endl; \
		double RMSE; \
		double RMSE_min; \
		double MSE_error; \
		TR_funct = residual_MSE_fun; \
		column_vector x(no_parameters); \
		column_vector x_min(no_parameters); \
		p_vals.clear(); \
		clock_t startTime = clock(); \
		unsigned int p_i; \
		fixed_parameters.clear(); \
		unsigned int c_i; \
		double C; \
		double SImax; \
		double SImin; \
		double SIrange; \
		for(unsigned int si_i = 0; si_i < SIdata.size(); si_i++) { \
			RMSE_min = std::numeric_limits<double>::infinity(); \
			SImax = std::numeric_limits<double>::min(); \
			SImin = std::numeric_limits<double>::max(); \
			B_0 = SIdata[si_i][0]; \
			if(B_0 == 0) { \
				B_0 = 1; \
			} \
			for(int b_i = no_bvalues-NUMBER_OF_TAIL_BVALUES; b_i < no_bvalues; b_i++) { \
				data_samples(b_i-(no_bvalues-NUMBER_OF_TAIL_BVALUES)).second = SIdata[si_i][b_i]/B_0; \
				if(SIdata[si_i][b_i]/B_0 < SImin) SImin = data_samples(b_i-(no_bvalues-NUMBER_OF_TAIL_BVALUES)).second; \
				if(SIdata[si_i][b_i]/B_0 > SImax) SImax = data_samples(b_i-(no_bvalues-NUMBER_OF_TAIL_BVALUES)).second; \
			} \
			SIrange = SImax - SImin; \
			std::fill(x_min.begin(), x_min.end(), std::numeric_limits<double>::quiet_NaN()); \
			for(c_i = 0; c_i < (int)C_values.size(); c_i++) { \
				C = SImin+C_values[c_i]*SIrange; \
				for(p_i = 0; p_i < (int)pinit.size(); p_i++) { \
					x(0) = pinit[p_i][1]; \
					x(1) = C; \
					try { \
						find_min_box_constrained_fixed(bfgs_search_strategy(), objective_delta_stop_strategy(DELTA_STOP, MAX_NUMBER_OF_ITERATIONS)OPT_VERBOSE, TR_evaluation, derivative(TR_evaluation), x, x_lo, x_hi); \
					} catch(...) { \
						cout << "exception occurred in optimization" << endl; \
						continue; \
					} \
					MSE_error = TR_evaluation(x); \
					RMSE = std::sqrt(MSE_error); \
					if(RMSE < RMSE_min) { \
						RMSE_min = RMSE; \
						x_min  = x; \
					} \
				} \
			} \
			std::vector<double> p_val = std::vector<double>(no_parameters+2); \
			if(B_0 == 0) { \
				p_val[0] = std::numeric_limits<double>::quiet_NaN(); \
			} else { \
				p_val[0] = (B_0-x_min(1)*B_0)/B_0; \
			} \
			p_val[1] = x_min(0); \
			p_val[2] = x_min(1)*B_0; \
			p_val[3] = RMSE_min*B_0; \
			if(si_i % PRINTOUT_STEP == 0) { \
				cout << "[" << si_i << "]:"; \
				cout << "f(analytical)=" << p_val[0] << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
				cout << info.parameterNames[1] << "=" << p_val[1] << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
				cout << "C=" << p_val[2] << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
				cout << "RMSE=" << p_val[3] << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
				cout << endl; \
			} \
			p_vals.push_back(p_val); \
		} \
		info.executiontime = double( clock() - startTime ) / (double)CLOCKS_PER_SEC; \
		cout << "Execution finished in " << info.executiontime << " seconds" << endl; \
	} \

//Log-linearized fit according to 
// [1] Y. Pang, B. Turkbey, M. Bernardo, J. Kruecker, S. Kadoury, M.J. Merino, et al., 
// Intravoxel incoherent motion MR imaging for prostate cancer: an evaluation of perfusion 
// fraction and diffusion coefficient derived from different b-value combinations., Magn. Reson. Med. 69 (2013) 553–62. 
// doi:10.1002/mrm.24277. 
#define FIT_MI_IMPLEMENTATION_LOGLINASYMPTOTICBIEXP(NAME, NUMBER_OF_BVALUES, NUMBER_OF_TAIL_BVALUES) \
	void NAME(std::vector<double> bvalues, const std::vector<std::vector<double> >& SIdata, double (*residual_MSE_fun)(const matrix<std::pair<double, double>, 0, 1>&, const matrix<double,0,1>&), std::vector<std::vector<double> >& p_vals, executionInfo& info) { \
		if((int)SIdata.size() == 0) { \
			throw "SIdata length was 0"; \
		} \
		data_samples = matrix<std::pair<double, double>, 0, 1>(NUMBER_OF_TAIL_BVALUES); \
		double B_0; \
		for(int b_i = NUMBER_OF_BVALUES-NUMBER_OF_TAIL_BVALUES; b_i < NUMBER_OF_BVALUES; b_i++) { \
			data_samples(b_i-(NUMBER_OF_BVALUES-NUMBER_OF_TAIL_BVALUES)).first = bvalues[b_i]; \
		} \
		cout << 1 << " initializations for 2 parameters " << endl; \
		int no_bvalues = (int)SIdata[0].size(); \
		cout << no_bvalues << " b-values" << endl; \
		double RMSE_min; \
		double MSE_error; \
		TR_funct = residual_MSE_fun; \
		column_vector x_min(2); \
		p_vals.clear(); \
		clock_t startTime = clock(); \
		fixed_parameters.clear(); \
		int p_i; \
		double SImax; \
		double SImin; \
		double SIrange; \
		double f; \
		double D; \
		for(unsigned int si_i = 0; si_i < SIdata.size(); si_i++) { \
			RMSE_min = std::numeric_limits<double>::infinity(); \
			SImax = std::numeric_limits<double>::min(); \
			SImin = std::numeric_limits<double>::max(); \
			B_0 = SIdata[si_i][0]; \
			if(B_0 == 0) { \
				B_0 = 1; \
			} \
			for(int b_i = 0; b_i < no_bvalues; b_i++) { \
				if(SIdata[si_i][b_i]/B_0 < SImin) SImin = std::log(SIdata[si_i][b_i]/B_0); \
				if(SIdata[si_i][b_i]/B_0 > SImax) SImax = std::log(SIdata[si_i][b_i]/B_0); \
			} \
			SIrange = SImax - SImin; \
			for(int b_i = no_bvalues-NUMBER_OF_TAIL_BVALUES; b_i < no_bvalues; b_i++) { \
				data_samples(b_i-(no_bvalues-NUMBER_OF_TAIL_BVALUES)).second = std::log(SIdata[si_i][b_i]/B_0); \
			} \
			std::fill(x_min.begin(), x_min.end(), std::numeric_limits<double>::quiet_NaN()); \
			try { \
				SimpleLinearRegression(data_samples, D, f); \
				x_min(0) = -f; \
				x_min(1) = -D; \
			} catch(...) { \
				cout << "exception occurred in optimization" << endl; \
				continue; \
			} \
			MSE_error = TR_evaluation(x_min); \
			RMSE_min = std::sqrt(MSE_error); \
			std::vector<double> p_val = std::vector<double>(3); \
			p_val[0] = x_min(0); \
			p_val[1] = x_min(1); \
			p_val[2] = RMSE_min; \
			if(si_i % PRINTOUT_STEP == 0) { \
				cout << "[" << si_i << "]:"; \
				for(p_i = 0; p_i < 2; p_i++) { \
					cout << info.parameterNames[p_i] << "=" << p_val[p_i] << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
				} \
				cout << " RMSE="<< RMSE_min << std::setprecision(OUTPUT_PRECISION) << std::fixed << endl; \
			} \
			p_vals.push_back(p_val); \
		} \
		info.executiontime = double( clock() - startTime ) / (double)CLOCKS_PER_SEC; \
		cout << "Execution finished in " << info.executiontime << " seconds" << endl; \
	} \

// ----------------------------------------------------------------------------------------
#define FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED(NAME, NUMBER_OF_BVALUES, NUMBER_OF_TAIL_BVALUES) \
void NAME(std::vector<double> bvalues, const std::vector<std::vector<double> >& SIdata, const std::vector<std::vector<double> >& pinit, const std::vector<std::vector<double> >& pinit2, double (*residual_fun)(const std::pair<double, double>&, const matrix<double,1,1>&), double (*residual_fun_2nd)(const std::pair<double, double>&, const matrix<double,2,1>&), std::vector<std::vector<double> >& p_vals, executionInfo& info) { \
    if((int)SIdata.size() == 0) \
        throw "SIdata length was 0"; \
    matrix<std::pair<double, double>, NUMBER_OF_BVALUES, 1> data_samples; \
    for(int b_i = 0; b_i < NUMBER_OF_BVALUES; b_i++) data_samples(b_i).first = bvalues[b_i]; \
    matrix<std::pair<double, double>, NUMBER_OF_TAIL_BVALUES, 1> data_tail_samples; \
    for(int b_i = NUMBER_OF_BVALUES-NUMBER_OF_TAIL_BVALUES; b_i < NUMBER_OF_BVALUES; b_i++) data_tail_samples(b_i-(NUMBER_OF_BVALUES-NUMBER_OF_TAIL_BVALUES)).first = bvalues[b_i]; \
    cout << ((int)pinit.size()) << " initializations for 1+2 parameters" << endl; \
    const int no_parameters = 3; \
    int no_bvalues = (int)SIdata[0].size(); \
    cout << no_bvalues << " b-values" << endl; \
    double RMSE; \
    double RMSE_min; \
    double absolute_error; \
    matrix<double,1,1> x_1st; \
    matrix<double,1,1> x_min_1st; \
    matrix<double,2,1> x_2nd; \
    matrix<double,3,1> x; \
    matrix<double,3,1> x_min; \
    p_vals.clear(); \
    clock_t startTime = clock(); \
    unsigned int p_i; \
    fixed_parameters.clear(); \
    fixed_parameters.push_back(0); \
    for(unsigned int si_i = 0; si_i < (int)SIdata.size(); si_i++) { \
        for(int b_i = 0; b_i < no_bvalues; b_i++) \
            data_samples(b_i).second = SIdata[si_i][b_i]/SIdata[si_i][0]; \
        for(int b_i = no_bvalues-NUMBER_OF_TAIL_BVALUES; b_i < no_bvalues; b_i++) \
            data_tail_samples(b_i-(no_bvalues-NUMBER_OF_TAIL_BVALUES)).second = SIdata[si_i][b_i]/SIdata[si_i][no_bvalues-NUMBER_OF_TAIL_BVALUES]; \
        RMSE_min = std::numeric_limits<double>::infinity(); \
        std::fill(x_min_1st.begin(), x_min_1st.end(), std::numeric_limits<double>::quiet_NaN()); \
        for(p_i = 0; p_i < (int)pinit.size(); p_i++) { \
            if(SIdata[si_i][no_bvalues-1] == 0 || (SIdata[si_i][no_bvalues-NUMBER_OF_TAIL_BVALUES]/SIdata[si_i][no_bvalues-1]) <= 0.0) \
                x_1st(0) = 0.0; \
            else \
                x_1st(0) = (std::log(SIdata[si_i][no_bvalues-NUMBER_OF_TAIL_BVALUES]/SIdata[si_i][no_bvalues-1])/(bvalues[no_bvalues-1]-bvalues[no_bvalues-NUMBER_OF_TAIL_BVALUES]))*pinit[p_i][0]; \
            cout << "Ds="<< p_i << " " << x_1st(0) << std::setprecision(OUTPUT_PRECISION) << std::fixed << " init" << endl; \
            solve_least_squares_lm(objective_delta_stop_strategy(DELTA_STOP, MAX_NUMBER_OF_ITERATIONS)OPT_VERBOSE, residual_fun, derivative(residual_fun),	data_tail_samples, x_1st); \
            RMSE = 0; \
            for(int j = 0; j < no_bvalues; j++) { \
                absolute_error = residual_fun(data_samples(j), x_1st); \
                RMSE += std::pow(absolute_error, 2.0); \
            } \
            RMSE /= (double)no_bvalues; \
            RMSE = std::sqrt(RMSE); \
            if(RMSE < RMSE_min) { \
                RMSE_min = RMSE; \
                x_min_1st = x_1st; \
            } \
            cout << "Ds=" << x_min_1st(0) << std::setprecision(OUTPUT_PRECISION) << std::fixed << endl; \
        } \
        RMSE_min = std::numeric_limits<double>::infinity(); \
        std::fill(x_min.begin(), x_min.end(), std::numeric_limits<double>::quiet_NaN()); \
        for(p_i = 0; p_i < (int)pinit2.size(); p_i++) { \
            x_2nd(0) = (1.0-x_min_1st(0))*pinit2[p_i][0]; \
            if(SIdata[si_i][no_bvalues-NUMBER_OF_TAIL_BVALUES] <= 0.0) \
                x_2nd(1) = 0.0; \
            else \
                x_2nd(1) = (std::log(SIdata[si_i][no_bvalues-NUMBER_OF_TAIL_BVALUES])/(bvalues[no_bvalues-NUMBER_OF_TAIL_BVALUES]-bvalues[0]))*pinit2[p_i][1]; \
            solve_least_squares_lm(objective_delta_stop_strategy(DELTA_STOP, MAX_NUMBER_OF_ITERATIONS)OPT_VERBOSE, residual_fun_2nd, derivative(residual_fun_2nd),	data_samples, x_2nd); \
            x(0) = x_2nd(0); \
            x(1) = x_2nd(1); \
            x(2) = x_1st(0); \
            RMSE = 0; \
            for(int j = 0; j < no_bvalues; j++) { \
                absolute_error = residual_fun_2nd(data_samples(j), x_2nd); \
                RMSE += std::pow(absolute_error, 2.0); \
            } \
            RMSE /= (double)no_bvalues; \
            RMSE = std::sqrt(RMSE); \
            if(RMSE < RMSE_min) { \
                RMSE_min = RMSE; \
                x_min = x; \
            } \
            cout << "f=" << x(0) << "Df=" << x(1) << std::setprecision(OUTPUT_PRECISION) << std::fixed << endl; \
        } \
        std::vector<double> p_val = std::vector<double>(no_parameters+1); \
        if(si_i % PRINTOUT_STEP == 0) { \
            cout << "[" << si_i << "]:"; \
            for(p_i = 0; p_i < no_parameters; p_i++) { \
                cout << info.parameterNames[p_i] << "=" << x_min(p_i) << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
            } \
            cout << " RMSE="<< RMSE_min << std::setprecision(OUTPUT_PRECISION) << std::fixed << endl; \
        } \
        for(p_i = 0; p_i < no_parameters; p_i++) { \
            p_val[p_i] = x_min(p_i); \
        } \
        p_val[no_parameters] = RMSE_min; \
        p_vals.push_back(p_val); \
        info.executiontime = double( clock() - startTime ) / (double)CLOCKS_PER_SEC; \
        cout << "Execution finished in " << info.executiontime << " seconds" << endl; \
    } \
} \

// ----------------------------------------------------------------------------------------
//Asymptotic fitting
// [1] J. Pekar, C.T. Moonen, P.C. van Zijl,
// On the precision of diffusion/perfusion imaging by gradient sensitization.,
// Magn. Reson. Med. 23 (1992) 122–9.
//Macro to define LM fit with fixed number of parameter values, implementation
#define FIT_MI_IMPLEMENTATION_LOGINIT_NORMALIZED_BFGS(NAME, NUMBER_OF_BVALUES, NUMBER_OF_TAIL_BVALUES) \
void NAME(std::vector<double> bvalues, const std::vector<std::vector<double> >& SIdata, const std::vector<std::vector<double> >& pinit, const std::vector<std::vector<double> >& pinit2, double (*residual_MSE_fun)(const matrix<std::pair<double, double>, 0, 1>&, const matrix<double, 0, 1>&), double (*residual_MSE_fun_2nd)(const matrix<std::pair<double, double>, 0, 1>&, const matrix<double, 0, 1>&), std::vector<std::vector<double> >& p_vals, executionInfo& info, std::vector<double>& C_values, const column_vector& x_lo1, const column_vector& x_hi1, const column_vector& x_lo2, const column_vector& x_hi2) { \
    if((int)SIdata.size() == 0) { \
        throw "SIdata length was 0"; \
    } \
    data_samples = matrix<std::pair<double, double>, 0, 1>(NUMBER_OF_TAIL_BVALUES); \
    data_samples2 = matrix<std::pair<double, double>, 0, 1>(NUMBER_OF_BVALUES); \
    double B_0; \
    for(int b_i = 0; b_i < NUMBER_OF_BVALUES; b_i++) { \
        data_samples2(b_i).first = bvalues[b_i]; \
    } \
    for(int b_i = (NUMBER_OF_BVALUES-NUMBER_OF_TAIL_BVALUES); b_i < NUMBER_OF_BVALUES; b_i++) { \
        data_samples(b_i-(NUMBER_OF_BVALUES-NUMBER_OF_TAIL_BVALUES)).first = bvalues[b_i]; \
    } \
    cout << ((int)pinit.size()) << " initializations for 1+3 parameters" << endl; \
    const int no_parameters = 4; \
    int no_bvalues = (int)SIdata[0].size(); \
    cout << no_bvalues << " b-values" << endl; \
    double RMSE; \
    double RMSE_min; \
    double RMSE_min_1st; \
    double MSE_error; \
    TR_funct = residual_MSE_fun; \
    TR_funct2 = residual_MSE_fun_2nd; \
    matrix<double,2,1> x_1st; \
    matrix<double,2,1> x_min_1st; \
    matrix<double,3,1> x_2nd; \
    matrix<double,4,1> x; \
    matrix<double,4,1> x_min; \
    p_vals.clear(); \
    clock_t startTime = clock(); \
    unsigned int p_i; \
    fixed_parameters.clear(); \
    fixed_parameters.push_back(0); \
    unsigned int c_i; \
    double C; \
    double SImax; \
    double SImin; \
    double SIrange; \
    for(unsigned int si_i = 0; si_i < (int)SIdata.size(); si_i++) { \
        RMSE_min = std::numeric_limits<double>::infinity(); \
        SImax = std::numeric_limits<double>::min(); \
        SImin = std::numeric_limits<double>::max(); \
        B_0 = SIdata[si_i][0]; \
        if(B_0 == 0) { \
            B_0 = 1; \
        } \
        for(int b_i = no_bvalues-NUMBER_OF_TAIL_BVALUES; b_i < no_bvalues; b_i++) { \
            data_samples(b_i-(no_bvalues-NUMBER_OF_TAIL_BVALUES)).second = SIdata[si_i][b_i]/B_0; \
            if(data_samples(b_i-(no_bvalues-NUMBER_OF_TAIL_BVALUES)).second < SImin) SImin = data_samples(b_i-(no_bvalues-NUMBER_OF_TAIL_BVALUES)).second; \
            if(data_samples(b_i-(no_bvalues-NUMBER_OF_TAIL_BVALUES)).second > SImax) SImax = data_samples(b_i-(no_bvalues-NUMBER_OF_TAIL_BVALUES)).second; \
        } \
        for(int b_i = 0; b_i < no_bvalues; b_i++) { \
            data_samples2(b_i).second = SIdata[si_i][b_i]/B_0; \
        } \
        SIrange = SImax - SImin; \
        std::fill(x_min_1st.begin(), x_min_1st.end(), std::numeric_limits<double>::quiet_NaN()); \
        RMSE_min_1st = std::numeric_limits<double>::infinity(); \
        for(c_i = 0; c_i < (int)C_values.size(); c_i++) { \
            C = SImin+C_values[c_i]*SIrange; \
            for(p_i = 0; p_i < (int)pinit.size(); p_i++) { \
                if(SIdata[si_i][no_bvalues-1] == 0 || (SIdata[si_i][no_bvalues-NUMBER_OF_TAIL_BVALUES]/SIdata[si_i][no_bvalues-1]) <= 0.0) \
                    x_1st(0) = 0.0; \
                else \
                    x_1st(0) = (std::log(SIdata[si_i][no_bvalues-NUMBER_OF_TAIL_BVALUES]/SIdata[si_i][no_bvalues-1])/(bvalues[no_bvalues-1]-bvalues[no_bvalues-NUMBER_OF_TAIL_BVALUES]))*pinit[p_i][0]; \
                x_1st(1) = C; \
                find_min_box_constrained_fixed(bfgs_search_strategy(), objective_delta_stop_strategy(DELTA_STOP, MAX_NUMBER_OF_ITERATIONS)OPT_VERBOSE, TR_evaluation, derivative(TR_evaluation), x_1st, x_lo1, x_hi1); \
                MSE_error = TR_evaluation(x_1st); \
                RMSE = std::sqrt(MSE_error); \
                if(x_1st(0) < 0) x_1st(0) = 0; \
                if(RMSE < RMSE_min_1st) { \
                    RMSE_min_1st = RMSE; \
                    x_min_1st  = x_1st; \
                } \
            } \
        } \
        RMSE_min = std::numeric_limits<double>::infinity(); \
        std::fill(x_min.begin(), x_min.end(), std::numeric_limits<double>::quiet_NaN()); \
        for(c_i = 0; c_i < (int)C_values.size(); c_i++) { \
            C = SImin+C_values[c_i]*SIrange; \
            for(p_i = 0; p_i < (int)pinit2.size(); p_i++) { \
                x_2nd(0) = (1.0-x_min_1st(0))*pinit2[p_i][0]; \
                if(SIdata[si_i][no_bvalues-NUMBER_OF_TAIL_BVALUES] <= 0.0) \
                    x_2nd(1) = 0.0; \
                else \
                    x_2nd(1) = (std::log(SIdata[si_i][no_bvalues-NUMBER_OF_TAIL_BVALUES])/(bvalues[no_bvalues-NUMBER_OF_TAIL_BVALUES]-bvalues[0]))*pinit2[p_i][1]; \
                fixed_parameters[0] = x_min_1st(0); \
                x_2nd(2) = C; \
                find_min_box_constrained_fixed(bfgs_search_strategy(), objective_delta_stop_strategy(DELTA_STOP, MAX_NUMBER_OF_ITERATIONS)OPT_VERBOSE, TR_evaluation2, derivative(TR_evaluation2), x_2nd, x_lo2, x_hi2); \
                x(0) = x_2nd(0); \
                x(1) = x_2nd(1); \
                x(2) = x_min_1st(0); \
                x(3) = x_2nd(2); \
                MSE_error = TR_evaluation2(x_2nd); \
                RMSE = std::sqrt(MSE_error); \
                if(RMSE < RMSE_min) { \
                    RMSE_min = RMSE; \
                    x_min  = x; \
                } \
            } \
        } \
        std::vector<double> p_val = std::vector<double>(no_parameters+1); \
        x_min(no_parameters-1) = x_min(no_parameters-1)*B_0; \
        if(si_i % PRINTOUT_STEP == 0) { \
            cout << "[" << si_i << "]:"; \
            for(p_i = 0; p_i < no_parameters; p_i++) { \
                cout << info.parameterNames[p_i] << "=" << x_min(p_i) << " " << std::setprecision(OUTPUT_PRECISION) << std::fixed; \
            } \
            cout << " RMSE="<< RMSE_min << std::setprecision(OUTPUT_PRECISION) << std::fixed << endl; \
        } \
        for(p_i = 0; p_i < no_parameters; p_i++) { \
            p_val[p_i] = x_min(p_i); \
        } \
        p_val[no_parameters] = RMSE_min*B_0; \
        p_vals.push_back(p_val); \
    } \
    info.executiontime = double( clock() - startTime ) / (double)CLOCKS_PER_SEC; \
    cout << "Execution finished in " << info.executiontime << " seconds" << endl; \
} \

