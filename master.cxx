#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include "binomial.h"
#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include "datareading.h"
#include "refreading.h"
#include "weightsreading.h"
#include <omp.h>
#include "TCanvas.h"

using namespace std;


int main(int argc, char** argv) {


/**
 * @dr: getting the information from the interpolationfile created by prof2-ipol
 * @weights: this vector contains all the weights for the bins of the interpolation
 * @rr: this object stores the information about the reference bin values
 * @argv: this is the name of the interpolationfile and maybe the name of the weightingfile
 */
	if(argc > 1)
	{
		
		datareading dr(argv[1]);
		
		vector<int> weights;
		//if there is a weightingfile given, the weightsreading-object will be created
		//@wr stores all the weights afterwards and sorts them in the same order as dr has them
		if(argc > 2)
		{
			weightsreading wr(argv[2]);
			wr.sort(dr);	
			weights = wr.weights;
		}
		else
		{
			//if no weightingfile is given, all the weights will be set to 1
			for(size_t i = 0; i < dr.analysis.size(); i++)
				weights.push_back(1);
		}
		//the constructor of refreading changes the names of dr.analysis by appending ".root", therefore it must (!) be called after the weights were read	
		refreading rr(dr);
		
/**
 * the following parts create and use a BAT-Object based on the informations of the previous data
 */
		
		//setting up the omp enviroment
		//the chains are distributed in threads, as well as the for-loop in fit_func(). therefore its nested.
		//the dynamic-part is needed so that threads are permanent available and don't have to be deleted and recreated (improves runtime significantly)
		omp_set_dynamic(1);
		omp_set_nested(1);

 		// set nicer style for drawing than the ROOT default
    		BCAux::SetStyle();
		
    		// open log file
    		BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    		// create new binomial object
    		binomial m("Name_Me", dr.var_names, dr.fit_num_params, dr.fit_params, dr.fit_error, dr.fit_order, dr.var_minvalue, dr.var_maxvalue, rr.ref_values, rr.ref_values_errors, weights);

    		// set precision
    		m.SetPrecision(BCEngineMCMC::kMedium);
		
		//vector<double> start = {0.64,1.27,0.90,0.289,0.138,0.404};
		//vector<vector<double>> startingpoints;
		//for(size_t i = 0; i < m.GetNChains(); i++)
		//	startingpoints.push_back(start);

		//m.SetInitialPositions(startingpoints);

		//m.SetMinimumEfficiency(0.5);
		//m.SetMaximumEfficiency(0.75);
		//m.SetNChains(4);
		//m.SetScaleFactorLowerLimit(0.01);

		//the number of iterations should assure, that a convergence can be reached and therefore an optimal mainrun can be performed
		m.SetNIterationsPreRunMax(100 * m.GetNIterationsPreRunMax());
    		BCLog::OutSummary("Test model created");

		//m.Normalize();

		// run MCMC, marginalizing posterior
    		m.MarginalizeAll(BCIntegrate::kMargMetropolis);
		
    		// run mode finding; by default using Minuit
  		m.FindMode(m.GetBestFitParameters());
		
  		// draw all marginalized distributions into a PDF file
  		m.PrintAllMarginalized(m.GetSafeName() + "_plots.pdf");

  		// print summary plots
   		m.PrintParameterPlot(m.GetSafeName() + "_parameters.pdf");
   		m.PrintCorrelationPlot(m.GetSafeName() + "_correlation.pdf");
   		m.PrintCorrelationMatrix(m.GetSafeName() + "_correlationMatrix.pdf");
   		m.PrintKnowledgeUpdatePlots(m.GetSafeName() + "_update.pdf");

  		// print results of the analysis into a text file
  		m.PrintSummary();

  		// close log file
   		BCLog::OutSummary("Exiting");
  		BCLog::CloseLog();
		
		m.chainhistoryplot();
		//m.binvaluedistribution();
		//m.chi2plot();
	}
	else
		cout << "no filenames specified" << endl;

	return 0;
}
