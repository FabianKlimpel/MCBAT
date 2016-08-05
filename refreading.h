#ifndef __REFREADING__H
#define __REFREADING__H

#include <iostream>
#include <string>
#include <stdlib.h>
#include <vector>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "datareading.h"

using namespace std;

/**
* class for reading the reference data of a root file
*/
class refreading
{
public:
	refreading(string analysis, string oberservable);
	refreading(datareading &dr);

/**
 * variables to store the reference data
 * @num_bins: stores the number of bins for an observable
 * @ref_values, @ref_errors: storage for the value in each bin of every histogram and respectively their errors
 */
	
	vector<double> ref_values, ref_values_errors;

private:

	int ref_num_bins;
	
	
};

#endif
