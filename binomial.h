// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__BINOMIAL__H
#define __BAT__BINOMIAL__H

#include <BAT/BCModel.h>

#include <string>
#include <vector>
#include <math.h>
#include <time.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>

// This is a binomial header file.
// Model source code is located in file binomial/binomial.cxx

// ---------------------------------------------------------
class binomial : public BCModel
{

public:

   	// Constructor
    	binomial(const std::string& model_name, const std::vector<std::string> var_names, size_t fit_num_params, const std::vector<double> fit_params, const std::vector<double> fit_errors, const int fit_order, const std::vector<double> var_minvalue, const std::vector<double> var_maxvalue, std::vector<double> ref_values, std::vector<double> ref_values_errors, std::vector<int> weights);

   	// Destructor
    	~binomial();

   	// Overload LogLikelihood to implement model
    	double LogLikelihood(const std::vector<double>& pars);
	double fit_func(const std::vector<double>& parameters, const int& index);
	double fit_err_func(const std::vector<double>& parameters, const int& index);
	void setFactors(const int& index);
	void setFactorsErrors(const int& index);
	void updateHistory();
	void binvaluedistribution();
	void chi2plot();
	void chainhistoryplot();
	void chi2distribution();

  	// Overload LogAprioriProbability if not using built-in 1D priors
     	//double LogAPrioriProbability(const std::vector<double> & pars);

    // Overload CalculateObservables if using observables
    // void CalculateObservables(const std::vector<double> & pars);
std::vector<std::vector<std::vector<double>>> MCMC_history;
private:

	size_t fit_num_params;
	std::vector<double> fit_params, fit_errors, var_minvalue, var_maxvalue, ref_values, ref_values_errors, offset, err_offset;
	std::vector<std::string> var_names;
	std::vector<int> weights, index_shifts;
	std::vector<std::vector<double>>  factors1, factors2, factors3, factors4, err_factors1, err_factors2, err_factors3, err_factors4;
	int fit_order, iteration_number = 0;
	

};
// ---------------------------------------------------------

#endif
