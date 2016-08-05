#ifndef __DATAREADING__H
#define __DATAREADING__H

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>

using namespace std;

/**
* class for reading the interpolation file
* the reading process is divided in 4 void functions
* the public parameters are for the external access by the BCModel constructor
*/
class datareading
{
public:
	datareading(char* interpolation);
	void parameter_reading();
	void var_min_max_reading();
	void fit_params_reading();
	void analysis_name();
	int get_num_histos();


	size_t dimension = 0, fit_num_params = 0;
	vector<double> fit_params, fit_error, var_minvalue, var_maxvalue;
	int fit_order = 0;
	vector<string> var_names, analysis, observable; 
	vector<size_t> changing_histo;
	
private:

	size_t pos_begin, pos_end, length;
	string line, parameternames;
};

#endif
