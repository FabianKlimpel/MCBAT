// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "binomial.h"
#include <vector>
#include <string>
#include <BAT/BCMath.h>
#include <omp.h>
/**
 * interpolation function for 4th order interpolation
 * @fit_params: fitted parameters of the interpolation
 * @start_index: Due to the construction of one vector, that consists of all interpolationparameters,
 * this parameter represents the first index of the parameters of a certain bin
 * @var_minvalue , @var_maxvalue: these represent the boarders of the hypercube for the values
 * they are needed, because the values were shifted for the interpolation by the center of the hypercube
 */
double fit_func_redundant(const std::vector<double>& fit_params, const size_t start_index, const std::vector<double>& var_minvalue, const std::vector<double>& var_maxvalue, const std::vector<double>& parameters, const std::vector<int>& index_shifts){

	int loop_counter1 = 0, loop_counter2 = 0, loop_counter3 = 0, loop_counter4 = 0;
	size_t i, j, k, m;
	size_t size = parameters.size();

	//0. order of interpolation
	double result = fit_params[start_index];

	
	for(m = size; m >= 1; m--)
	{
		//1. order of interpolation
		result += fit_params[start_index + 1 + loop_counter1] 
				* (parameters[m - 1] - var_minvalue[m - 1]) / (var_maxvalue[m - 1] - var_minvalue[m - 1]);
		loop_counter1++;
		
		
		for(k = size; k >= m; k--)
		{
			//2. order of interpolation
			result += fit_params[start_index + 1 + index_shifts[1] + loop_counter2] 
					* (parameters[m - 1] - var_minvalue[m - 1]) / (var_maxvalue[m - 1] - var_minvalue[m - 1]) 
					* (parameters[k - 1] - var_minvalue[k - 1]) / (var_maxvalue[k - 1] - var_minvalue[k - 1]);
			loop_counter2++;

			
			for(j = size; j >= k; j--)
			{
				//3. order of interpolation
				result += fit_params[start_index + 1 + index_shifts[1] + index_shifts[2] + loop_counter3] 
						* (parameters[j - 1] - var_minvalue[j - 1]) / (var_maxvalue[j - 1] - var_minvalue[j - 1]) 
						* (parameters[k - 1] - var_minvalue[k - 1]) / (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
						* (parameters[m - 1] - var_minvalue[m - 1]) / (var_maxvalue[m - 1] - var_minvalue[m - 1]);
				loop_counter3++;
	
				
				for(i = size; i >= j; i--)
				{
					//4. order of interpolation
					result += fit_params[start_index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 
							* (parameters[i - 1] - var_minvalue[i - 1]) / (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							* (parameters[j - 1] - var_minvalue[j - 1]) / (var_maxvalue[j - 1] - var_minvalue[j - 1])
							* (parameters[k - 1] - var_minvalue[k - 1]) / (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							* (parameters[m - 1] - var_minvalue[m - 1]) / (var_maxvalue[m - 1] - var_minvalue[m - 1]);
					loop_counter4++;

				}
			}
		}
	}
	return result;
}


double fit_err_func_redundant(const std::vector<double>& fit_errors, const size_t start_index, const std::vector<double>& var_minvalue, const std::vector<double>& var_maxvalue, const std::vector<double>& parameters, const std::vector<int>& index_shifts){

	int loop_counter1 = 0, loop_counter2 = 0, loop_counter3 = 0, loop_counter4 = 0;
	size_t i, j, k, m;
	size_t size = parameters.size();

	//0. order of interpolation
	double result = fit_errors[start_index] * fit_errors[start_index];

	
	for(m = size; m >= 1; m--)
	{
		//1. order of interpolation
		result += fit_errors[start_index + 1 + loop_counter1] 
				* (parameters[m - 1] - var_minvalue[m - 1]) / (var_maxvalue[m - 1] - var_minvalue[m - 1]) 
				* fit_errors[start_index + 1 + loop_counter1] 
				* (parameters[m - 1] - var_minvalue[m - 1]) / (var_maxvalue[m - 1] - var_minvalue[m - 1]);
		loop_counter1++;

		
		for(k = size; k >= m; k--)
		{
			//2. order of interpolation
			result += fit_errors[start_index + 1 + index_shifts[1] + loop_counter2] 
					* (parameters[m - 1] - var_minvalue[m - 1]) / (var_maxvalue[m - 1] - var_minvalue[m - 1]) 
					* (parameters[k - 1] - var_minvalue[k - 1]) / (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
					* fit_errors[start_index + 1 + index_shifts[1] + loop_counter2] 
					* (parameters[m - 1] - var_minvalue[m - 1]) / (var_maxvalue[m - 1] - var_minvalue[m - 1]) 
					* (parameters[k - 1] - var_minvalue[k - 1]) / (var_maxvalue[k - 1] - var_minvalue[k - 1]);
			loop_counter2++;

		
			for(j = size; j >= k; j--)
			{
				//3. order of interpolation
				result += fit_errors[start_index + 1 + index_shifts[1] + index_shifts[2] + loop_counter3] 
						* (parameters[j - 1] - var_minvalue[j - 1]) / (var_maxvalue[j - 1] - var_minvalue[j - 1]) 
						* (parameters[k - 1] - var_minvalue[k - 1]) / (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
						* (parameters[m - 1] - var_minvalue[m - 1]) / (var_maxvalue[m - 1] - var_minvalue[m - 1]) 
						* fit_errors[start_index + 1 + index_shifts[1] + index_shifts[2] + loop_counter3] 
						* (parameters[j - 1] - var_minvalue[j - 1]) / (var_maxvalue[j - 1] - var_minvalue[j - 1]) 
						* (parameters[k - 1] - var_minvalue[k - 1]) / (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
						* (parameters[m - 1] - var_minvalue[m - 1]) / (var_maxvalue[m - 1] - var_minvalue[m - 1]);
				loop_counter3++;
				
			
				for(i = size; i >= j; i--)
				{
					//4. order of interpolation
					result += fit_errors[start_index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 
							* (parameters[i - 1] - var_minvalue[i - 1]) / (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							* (parameters[j - 1] - var_minvalue[j - 1]) / (var_maxvalue[j - 1] - var_minvalue[j - 1])
							* (parameters[k - 1] - var_minvalue[k - 1]) / (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							* (parameters[m - 1] - var_minvalue[m - 1]) / (var_maxvalue[m - 1] - var_minvalue[m - 1]) 
							* fit_errors[start_index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 
							* (parameters[i - 1] - var_minvalue[i - 1]) / (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							* (parameters[j - 1] - var_minvalue[j - 1]) / (var_maxvalue[j - 1] - var_minvalue[j - 1])
							* (parameters[k - 1] - var_minvalue[k - 1]) / (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							* (parameters[m - 1] - var_minvalue[m - 1]) / (var_maxvalue[m - 1] - var_minvalue[m - 1]);
					loop_counter4++;

				}
			}
		}
	}

	return result;
}

double binomial::fit_func(const std::vector<double>& parameters, const int& index){

	int loop_counter1 = 0, loop_counter2 = 0, loop_counter3 = 0, loop_counter4 = 0;
	size_t i, j, k, m;
	size_t size = parameters.size();

	//0. order of interpolation
	double result = offset[index];

	if(fit_order >= 1)
	for(m = size; m >= 1; m--)
	{
		//1. order of interpolation
		result += factors1[index][loop_counter1] * parameters[m - 1];

		loop_counter1++;
		
		if(fit_order >= 2)
		for(k = size; k >= m; k--)
		{
			//2. order of interpolation
			result += factors2[index][loop_counter2] * parameters[m - 1] * parameters[k - 1];
			result += factors2[index][loop_counter2 + 1] * parameters[m - 1];
			result += factors2[index][loop_counter2 + 2] * parameters[k - 1];
					
			loop_counter2 += 3;

			if(fit_order >= 3)
			for(j = size; j >= k; j--)
			{
				//3. order of interpolation
				result += factors3[index][loop_counter3] * parameters[j - 1] * parameters[k - 1] * parameters[m - 1];
				result += factors3[index][loop_counter3 + 1] * parameters[j - 1] * parameters[k - 1];
				result += factors3[index][loop_counter3 + 2] * parameters[j - 1] * parameters[m - 1];				
				result += factors3[index][loop_counter3 + 3] * parameters[j - 1];
				result += factors3[index][loop_counter3 + 4] * parameters[k - 1] * parameters[m - 1];
				result += factors3[index][loop_counter3 + 5] * parameters[k - 1];				
				result += factors3[index][loop_counter3 + 6] * parameters[m - 1];

				loop_counter3 += 7;

				if(fit_order == 4)
				for(i = size; i >= j; i--)
				{
					//4. order of interpolation
					result += factors4[index][loop_counter4] * parameters[i - 1] * parameters[j - 1] * parameters[k - 1] * parameters[m - 1];
					result += factors4[index][loop_counter4 + 1] * parameters[j - 1] * parameters[k - 1] * parameters[m - 1];
					result += factors4[index][loop_counter4 + 2] * parameters[i - 1] * parameters[j - 1] * parameters[k - 1];
					result += factors4[index][loop_counter4 + 3] * parameters[j - 1] * parameters[k - 1];
					result += factors4[index][loop_counter4 + 4] * parameters[i - 1] * parameters[j - 1] * parameters[m - 1];
					result += factors4[index][loop_counter4 + 5] * parameters[j - 1] * parameters[m - 1];
					result += factors4[index][loop_counter4 + 6] * parameters[i - 1] * parameters[j - 1];
					result += factors4[index][loop_counter4 + 7] * parameters[j - 1];
					result += factors4[index][loop_counter4 + 8] * parameters[i - 1] * parameters[k - 1] * parameters[m - 1];
					result += factors4[index][loop_counter4 + 9] * parameters[k - 1] * parameters[m - 1];
					result += factors4[index][loop_counter4 + 10] * parameters[i - 1] * parameters[k - 1];
					result += factors4[index][loop_counter4 + 11] * parameters[k - 1];
					result += factors4[index][loop_counter4 + 12] * parameters[i - 1] * parameters[m - 1];
					result += factors4[index][loop_counter4 + 13] * parameters[m - 1];
					result += factors4[index][loop_counter4 + 14] * parameters[i - 1];

					loop_counter4 += 15;

				}
			}
		}
	}
	return result;
}

void binomial::setFactors(const int& index){

	int loop_counter1 = 0, loop_counter2 = 0, loop_counter3 = 0, loop_counter4 = 0;
	size_t i, j, k, m;
	size_t size = var_names.size();
	std::vector<double> tmp_factors1, tmp_factors2, tmp_factors3, tmp_factors4;
	double tmp_offset;

	//0. order of interpolation
	tmp_offset = fit_params[index];

	if(fit_order >= 1)
	for(m = size; m >= 1; m--)
	{
		//1. order of interpolation
		tmp_offset -= fit_params[index + 1 + loop_counter1] 
				* var_minvalue[m - 1] / (var_maxvalue[m - 1] - var_minvalue[m - 1]);
		tmp_factors1.push_back(fit_params[index + 1 + loop_counter1] 
					/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));
		loop_counter1++;
		
		if(fit_order >= 2)
		for(k = size; k >= m; k--)
		{
			//2. order of interpolation
			tmp_offset += fit_params[index + 1 + index_shifts[1] + loop_counter2] 
					* var_minvalue[m - 1] / (var_maxvalue[m - 1] - var_minvalue[m - 1]) 
					* var_minvalue[k - 1] / (var_maxvalue[k - 1] - var_minvalue[k - 1]);

			tmp_factors2.push_back(fit_params[index + 1 + index_shifts[1] + loop_counter2] 
						/ (var_maxvalue[m - 1] - var_minvalue[m - 1]) 
						/ (var_maxvalue[k - 1] - var_minvalue[k - 1]));
			tmp_factors2.push_back(-fit_params[index + 1 + index_shifts[1] + loop_counter2] 
						/ (var_maxvalue[m - 1] - var_minvalue[m - 1]) 
						* var_minvalue[k - 1]
						/ (var_maxvalue[k - 1] - var_minvalue[k - 1]));
			tmp_factors2.push_back(-fit_params[index + 1 + index_shifts[1] + loop_counter2] 
						* var_minvalue[m - 1]
						/ (var_maxvalue[m - 1] - var_minvalue[m - 1]) 
						/ (var_maxvalue[k - 1] - var_minvalue[k - 1]));
			loop_counter2++;

			if(fit_order >= 3)
			for(j = size; j >= k; j--)
			{
				//3. order of interpolation
				tmp_offset -= fit_params[index + 1 + index_shifts[1] + index_shifts[2] + loop_counter3] 
						* var_minvalue[j - 1] / (var_maxvalue[j - 1] - var_minvalue[j - 1]) 
						* var_minvalue[k - 1] / (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
						* var_minvalue[m - 1] / (var_maxvalue[m - 1] - var_minvalue[m - 1]);

				tmp_factors3.push_back(fit_params[index + 1 + index_shifts[1] + index_shifts[2] + loop_counter3] 
						/ (var_maxvalue[j - 1] - var_minvalue[j - 1]) 
						/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
						/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

				tmp_factors3.push_back(-fit_params[index + 1 + index_shifts[1] + index_shifts[2] + loop_counter3] 
						/ (var_maxvalue[j - 1] - var_minvalue[j - 1]) 
						/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
						* var_minvalue[m - 1] 
						/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

				tmp_factors3.push_back(-fit_params[index + 1 + index_shifts[1] + index_shifts[2] + loop_counter3] 
						/ (var_maxvalue[j - 1] - var_minvalue[j - 1]) 
						* var_minvalue[k - 1] 
						/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 						
						/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

				tmp_factors3.push_back(fit_params[index + 1 + index_shifts[1] + index_shifts[2] + loop_counter3] 					
						/ (var_maxvalue[j - 1] - var_minvalue[j - 1]) 
						* var_minvalue[k - 1] 
						/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 	
						* var_minvalue[m - 1]					
						/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

				tmp_factors3.push_back(-fit_params[index + 1 + index_shifts[1] + index_shifts[2] + loop_counter3] 
						* var_minvalue[j - 1]					
						/ (var_maxvalue[j - 1] - var_minvalue[j - 1]) 					 
						/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 						
						/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

				tmp_factors3.push_back(fit_params[index + 1 + index_shifts[1] + index_shifts[2] + loop_counter3] 
						* var_minvalue[j - 1] 					
						/ (var_maxvalue[j - 1] - var_minvalue[j - 1]) 
						/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 	
						* var_minvalue[m - 1]					
						/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

				tmp_factors3.push_back(fit_params[index + 1 + index_shifts[1] + index_shifts[2] + loop_counter3] 
						* var_minvalue[j - 1] 					
						/ (var_maxvalue[j - 1] - var_minvalue[j - 1]) 
						* var_minvalue[k - 1]
						/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 								
						/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

				loop_counter3++;

				if(fit_order == 4)
				for(i = size; i >= j; i--)
				{
					//4. order of interpolation
					tmp_offset += fit_params[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 
							* var_minvalue[i - 1] / (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							* var_minvalue[j - 1] / (var_maxvalue[j - 1] - var_minvalue[j - 1])
							* var_minvalue[k - 1] / (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							* var_minvalue[m - 1] / (var_maxvalue[m - 1] - var_minvalue[m - 1]);

					tmp_factors4.push_back(fit_params[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(-fit_params[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 
							* var_minvalue[i - 1] 
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(-fit_params[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 	 
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							* var_minvalue[m - 1]
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(fit_params[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 
							* var_minvalue[i - 1]	 
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							* var_minvalue[m - 1]
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(-fit_params[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 	 
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							* var_minvalue[k - 1]
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(fit_params[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 
							* var_minvalue[i - 1]	 
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							* var_minvalue[k - 1]
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(fit_params[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4]  
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							* var_minvalue[k - 1]
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							* var_minvalue[m - 1]	
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(-fit_params[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4]  
							* var_minvalue[i - 1]
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							* var_minvalue[k - 1]
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							* var_minvalue[m - 1]	
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(-fit_params[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 	 
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							* var_minvalue[j - 1]
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(fit_params[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 
							* var_minvalue[i - 1]	 
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							* var_minvalue[j - 1]
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(fit_params[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 	 
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							* var_minvalue[j - 1]
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							* var_minvalue[m - 1]
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(-fit_params[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4]  
							* var_minvalue[i - 1]
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							* var_minvalue[j - 1]
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							* var_minvalue[m - 1]	
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(fit_params[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 	 
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							* var_minvalue[j - 1]
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							* var_minvalue[k - 1]
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(-fit_params[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4]  
							* var_minvalue[i - 1]
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							* var_minvalue[j - 1]
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							* var_minvalue[k - 1]
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(-fit_params[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4]  
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							* var_minvalue[j - 1]
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							* var_minvalue[k - 1]
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							* var_minvalue[m - 1]
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					loop_counter4++;

				}
			}
		}
	}
	offset.push_back(tmp_offset);
	factors1.push_back(tmp_factors1);
	factors2.push_back(tmp_factors2);
	factors3.push_back(tmp_factors3);
	factors4.push_back(tmp_factors4);
}



double binomial::fit_err_func(const std::vector<double>& parameters, const int& index){

	int loop_counter1 = 0, loop_counter2 = 0, loop_counter3 = 0, loop_counter4 = 0;
	size_t i, j, k, m;
	size_t size = parameters.size();

	//0. order of interpolation
	double result = err_offset[index] * err_offset[index], tmp_result;

	if(fit_order >= 1)
	for(m = size; m >= 1; m--)
	{
		//1. order of interpolation
		tmp_result = err_factors1[index][loop_counter1];
		tmp_result += err_factors1[index][loop_counter1 + 1] * parameters[m - 1];
		result += tmp_result * tmp_result;

		loop_counter1 += 2;
		
		if(fit_order >= 2)
		for(k = size; k >= m; k--)
		{
			//2. order of interpolation
			tmp_result = err_factors2[index][loop_counter2];
			tmp_result += err_factors2[index][loop_counter2 + 1] * parameters[m - 1] * parameters[k - 1];
			tmp_result += err_factors2[index][loop_counter2 + 2] * parameters[m - 1];
			tmp_result += err_factors2[index][loop_counter2 + 3] * parameters[k - 1];
			result += tmp_result * tmp_result;
		
			loop_counter2 += 4;

			if(fit_order >= 3)
			for(j = size; j >= k; j--)
			{
				//3. order of interpolation
				tmp_result = err_factors3[index][loop_counter3];
				tmp_result += err_factors3[index][loop_counter3 + 1] * parameters[j - 1] * parameters[k - 1] * parameters[m - 1];
				tmp_result += err_factors3[index][loop_counter3 + 2] * parameters[j - 1] * parameters[k - 1];
				tmp_result += err_factors3[index][loop_counter3 + 3] * parameters[j - 1] * parameters[m - 1];				
				tmp_result += err_factors3[index][loop_counter3 + 4] * parameters[j - 1];
				tmp_result += err_factors3[index][loop_counter3 + 5] * parameters[k - 1] * parameters[m - 1];
				tmp_result += err_factors3[index][loop_counter3 + 6] * parameters[k - 1];				
				tmp_result += err_factors3[index][loop_counter3 + 7] * parameters[m - 1];
				result += tmp_result * tmp_result;

				loop_counter3 += 8;

				if(fit_order == 4)
				for(i = size; i >= j; i--)
				{
					//4. order of interpolation
					tmp_result = err_factors4[index][loop_counter4];
					tmp_result += err_factors4[index][loop_counter4 + 1] * parameters[i - 1] * parameters[j - 1] * parameters[k - 1] * parameters[m - 1];
					tmp_result += err_factors4[index][loop_counter4 + 2] * parameters[j - 1] * parameters[k - 1] * parameters[m - 1];
					tmp_result += err_factors4[index][loop_counter4 + 3] * parameters[i - 1] * parameters[j - 1] * parameters[k - 1];
					tmp_result += err_factors4[index][loop_counter4 + 4] * parameters[j - 1] * parameters[k - 1];
					tmp_result += err_factors4[index][loop_counter4 + 5] * parameters[i - 1] * parameters[j - 1] * parameters[m - 1];
					tmp_result += err_factors4[index][loop_counter4 + 6] * parameters[j - 1] * parameters[m - 1];
					tmp_result += err_factors4[index][loop_counter4 + 7] * parameters[i - 1] * parameters[j - 1];
					tmp_result += err_factors4[index][loop_counter4 + 8] * parameters[j - 1];
					tmp_result += err_factors4[index][loop_counter4 + 9] * parameters[i - 1] * parameters[k - 1] * parameters[m - 1];
					tmp_result += err_factors4[index][loop_counter4 + 10] * parameters[k - 1] * parameters[m - 1];
					tmp_result += err_factors4[index][loop_counter4 + 11] * parameters[i - 1] * parameters[k - 1];
					tmp_result += err_factors4[index][loop_counter4 + 12] * parameters[k - 1];
					tmp_result += err_factors4[index][loop_counter4 + 13] * parameters[i - 1] * parameters[m - 1];
					tmp_result += err_factors4[index][loop_counter4 + 14] * parameters[m - 1];
					tmp_result += err_factors4[index][loop_counter4 + 15] * parameters[i - 1];
					result += tmp_result * tmp_result;

					loop_counter4 += 16;

				}
			}
		}
	}
	return result;
}


void binomial::setFactorsErrors(const int& index){

	int loop_counter1 = 0, loop_counter2 = 0, loop_counter3 = 0, loop_counter4 = 0;
	size_t i, j, k, m;
	size_t size = var_names.size();
	std::vector<double> tmp_factors1, tmp_factors2, tmp_factors3, tmp_factors4;
	

	//0. order of interpolation
	err_offset.push_back(fit_errors[index]);

	if(fit_order >= 1)
	for(m = size; m >= 1; m--)
	{
		//1. order of interpolation
		tmp_factors1.push_back(-fit_errors[index + 1 + loop_counter1] 
				* var_minvalue[m - 1] / (var_maxvalue[m - 1] - var_minvalue[m - 1]));
		tmp_factors1.push_back(fit_errors[index + 1 + loop_counter1] 
					/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));
		loop_counter1++;
		
		if(fit_order >= 2)
		for(k = size; k >= m; k--)
		{
			//2. order of interpolation
			tmp_factors2.push_back(fit_errors[index + 1 + index_shifts[1] + loop_counter2] 
					* var_minvalue[m - 1] / (var_maxvalue[m - 1] - var_minvalue[m - 1]) 
					* var_minvalue[k - 1] / (var_maxvalue[k - 1] - var_minvalue[k - 1]));

			tmp_factors2.push_back(fit_errors[index + 1 + index_shifts[1] + loop_counter2] 
						/ (var_maxvalue[m - 1] - var_minvalue[m - 1]) 
						/ (var_maxvalue[k - 1] - var_minvalue[k - 1]));
			tmp_factors2.push_back(-fit_errors[index + 1 + index_shifts[1] + loop_counter2] 
						/ (var_maxvalue[m - 1] - var_minvalue[m - 1]) 
						* var_minvalue[k - 1]
						/ (var_maxvalue[k - 1] - var_minvalue[k - 1]));
			tmp_factors2.push_back(-fit_errors[index + 1 + index_shifts[1] + loop_counter2] 
						* var_minvalue[m - 1]
						/ (var_maxvalue[m - 1] - var_minvalue[m - 1]) 
						/ (var_maxvalue[k - 1] - var_minvalue[k - 1]));
			loop_counter2++;

			if(fit_order >= 3)
			for(j = size; j >= k; j--)
			{
				//3. order of interpolation
				tmp_factors3.push_back(-fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + loop_counter3] 
						* var_minvalue[j - 1] / (var_maxvalue[j - 1] - var_minvalue[j - 1]) 
						* var_minvalue[k - 1] / (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
						* var_minvalue[m - 1] / (var_maxvalue[m - 1] - var_minvalue[m - 1]));

				tmp_factors3.push_back(fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + loop_counter3] 
						/ (var_maxvalue[j - 1] - var_minvalue[j - 1]) 
						/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
						/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

				tmp_factors3.push_back(-fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + loop_counter3] 
						/ (var_maxvalue[j - 1] - var_minvalue[j - 1]) 
						/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
						* var_minvalue[m - 1] 
						/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

				tmp_factors3.push_back(-fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + loop_counter3] 
						/ (var_maxvalue[j - 1] - var_minvalue[j - 1]) 
						* var_minvalue[k - 1] 
						/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 						
						/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

				tmp_factors3.push_back(fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + loop_counter3] 					
						/ (var_maxvalue[j - 1] - var_minvalue[j - 1]) 
						* var_minvalue[k - 1] 
						/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 	
						* var_minvalue[m - 1]					
						/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

				tmp_factors3.push_back(-fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + loop_counter3] 
						* var_minvalue[j - 1]					
						/ (var_maxvalue[j - 1] - var_minvalue[j - 1]) 					 
						/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 						
						/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

				tmp_factors3.push_back(fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + loop_counter3] 
						* var_minvalue[j - 1] 					
						/ (var_maxvalue[j - 1] - var_minvalue[j - 1]) 
						/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 	
						* var_minvalue[m - 1]					
						/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

				tmp_factors3.push_back(fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + loop_counter3] 
						* var_minvalue[j - 1] 					
						/ (var_maxvalue[j - 1] - var_minvalue[j - 1]) 
						* var_minvalue[k - 1]
						/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 								
						/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

				loop_counter3++;

				if(fit_order == 4)
				for(i = size; i >= j; i--)
				{
					//4. order of interpolation
					tmp_factors4.push_back(fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 
							* var_minvalue[i - 1] / (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							* var_minvalue[j - 1] / (var_maxvalue[j - 1] - var_minvalue[j - 1])
							* var_minvalue[k - 1] / (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							* var_minvalue[m - 1] / (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(-fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 
							* var_minvalue[i - 1] 
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(-fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 	 
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							* var_minvalue[m - 1]
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 
							* var_minvalue[i - 1]	 
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							* var_minvalue[m - 1]
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(-fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 	 
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							* var_minvalue[k - 1]
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 
							* var_minvalue[i - 1]	 
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							* var_minvalue[k - 1]
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4]  
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							* var_minvalue[k - 1]
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							* var_minvalue[m - 1]	
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(-fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4]  
							* var_minvalue[i - 1]
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							* var_minvalue[k - 1]
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							* var_minvalue[m - 1]	
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(-fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 	 
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							* var_minvalue[j - 1]
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 
							* var_minvalue[i - 1]	 
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							* var_minvalue[j - 1]
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 	 
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							* var_minvalue[j - 1]
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							* var_minvalue[m - 1]
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(-fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4]  
							* var_minvalue[i - 1]
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							* var_minvalue[j - 1]
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							* var_minvalue[m - 1]	
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4] 	 
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							* var_minvalue[j - 1]
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							* var_minvalue[k - 1]
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(-fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4]  
							* var_minvalue[i - 1]
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							* var_minvalue[j - 1]
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							* var_minvalue[k - 1]
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					tmp_factors4.push_back(-fit_errors[index + 1 + index_shifts[1] + index_shifts[2] + index_shifts[3] + loop_counter4]  
							/ (var_maxvalue[i - 1] - var_minvalue[i - 1]) 
							* var_minvalue[j - 1]
							/ (var_maxvalue[j - 1] - var_minvalue[j - 1])
							* var_minvalue[k - 1]
							/ (var_maxvalue[k - 1] - var_minvalue[k - 1]) 
							* var_minvalue[m - 1]
							/ (var_maxvalue[m - 1] - var_minvalue[m - 1]));

					loop_counter4++;

				}
			}
		}
	}

	err_factors1.push_back(tmp_factors1);
	err_factors2.push_back(tmp_factors2);
	err_factors3.push_back(tmp_factors3);
	err_factors4.push_back(tmp_factors4);
}

/**
* simple function that delivers a flat prior based on the actual sampled hypercube
* @var_minvalue & @var_maxvalue : these are the vectors that represent the hypercube
* @var_names : this is used for debugging in case of a 0-sized intervall
*/
double flat_prior(const std::vector<double>& parameters, const std::vector<double> var_minvalue, const std::vector<double> var_maxvalue, const std::vector<std::string> var_names)
{

	double result = 0;
	for(size_t i = 0; i < var_minvalue.size(); i++)
		if(parameters[i] > var_maxvalue[i] || parameters[i] < var_minvalue[i])
			return 0;
		else
			if(var_maxvalue[i] != var_minvalue[i])
				result *= 1/(var_maxvalue[i] - var_minvalue[i]);
			else
			{			
				std::cout << "error: no interval for parameter " << var_names[i] << std::endl;
				return 0;
			}
	return result;
}

// ---------------------------------------------------------
binomial::binomial(const std::string& model_name, const std::vector<std::string> var_names, size_t fit_num_params, const std::vector<double> fit_params, const std::vector<double> fit_errors, const int fit_order, const std::vector<double> var_minvalue, const std::vector<double> var_maxvalue, std::vector<double> ref_values, std::vector<double> ref_values_errors, std::vector<int> weights)
    : BCModel(model_name), fit_num_params(fit_num_params), fit_params(fit_params), fit_errors(fit_errors), var_minvalue(var_minvalue), var_maxvalue(var_maxvalue), ref_values(ref_values), ref_values_errors(ref_values_errors), weights(weights)
{

	size_t i, j, k, m, loop_counter = 0;

	for(i = 0; i < var_names.size(); i++)
    	 AddParameter(var_names[i], var_minvalue[i], var_maxvalue[i]);
     	
	this->var_names = var_names;

	this->fit_order = fit_order;

	SetPriorConstantAll();
	
	//TODO: sehr schlecht, aber sollte arbeit erledigen
	index_shifts.push_back(0);
	index_shifts.push_back(var_names.size());

	for(i = var_names.size(); i >= 1; i--)
		for(j = var_names.size(); j >= i; j--)
			loop_counter++;

	index_shifts.push_back(loop_counter);
	loop_counter = 0;

	for(k = var_names.size(); k >= 1; k--)
		for(j = var_names.size(); j >= k; j--)
			for(i = var_names.size(); i >= j; i--)
				loop_counter++;

	index_shifts.push_back(loop_counter);
	loop_counter = 0;

	for(m = var_names.size(); m >= 1; m--)
		for(k = var_names.size(); k >= m; k--)
			for(j = var_names.size(); j >= k; j--)
				for(i = var_names.size(); i >= j; i--)
					loop_counter++;

	index_shifts.push_back(loop_counter);

	for(i = 0; i < ref_values.size(); i++)
	{
		setFactors(i * fit_num_params);
		setFactorsErrors(i * fit_num_params);
	}
}

// ---------------------------------------------------------
binomial::~binomial()
{
    // destructor
}

// ---------------------------------------------------------
double binomial::LogLikelihood(const std::vector<double>& parameters)
{
	double logprob = 0;
	


	if(fMCMCPhase == BCEngineMCMC::kMainRun)
		updateHistory();

     	//Gaussian distribution
	#pragma omp parallel for reduction(+:logprob) schedule(guided) //num_threads(8) //private(begin,end) 
	for(size_t i = 0; i < ref_values.size(); i++)
	{
	
		logprob += weights[i] * BCMath::LogGaus(fit_func(parameters, i), ref_values[i], sqrt(ref_values_errors[i] * ref_values_errors[i] + fit_err_func(parameters, i)));
		
	}

	
	// return the log of the conditional probability p(data|pars).
    	return logprob;
}

void binomial::updateHistory(){
	#pragma omp critical (MCMC_history_update)
	{
	
	if(iteration_number != GetCurrentIteration())
	{
		for(int i = 0; i < GetCurrentIteration() - iteration_number; i++)
			MCMC_history.push_back(Getx());
		iteration_number = GetCurrentIteration();	
	}

	}
}

void binomial::binlogprobdistribution(){

	std::vector<double> mode = FindMode(GetBestFitParameters());
	
	double maximum = 0, minimum = 0;
	for(size_t i = 0; i < ref_values.size(); i++)
	{
		if(BCMath::LogGaus(fit_func(mode, i), ref_values[i], sqrt(ref_values_errors[i] * ref_values_errors[i] + fit_err_func(mode, i))) > maximum)
			maximum = BCMath::LogGaus(fit_func(mode, i), ref_values[i], sqrt(ref_values_errors[i] * ref_values_errors[i] + fit_err_func(mode, i)));
		if(BCMath::LogGaus(fit_func(mode, i), ref_values[i], sqrt(ref_values_errors[i] * ref_values_errors[i] + fit_err_func(mode, i))) < minimum)
			minimum = BCMath::LogGaus(fit_func(mode, i), ref_values[i], sqrt(ref_values_errors[i] * ref_values_errors[i] + fit_err_func(mode, i)));
	}


	TCanvas* tc;
	TH1D* th1d;

	tc = new TCanvas("binlogprobdistribution.pdf");
	tc->cd();
	th1d = new TH1D("binlogprogdistribution.pdf", "binlogprobdistribution.pdf", (int) (1.5 * sqrt(ref_values.size())), minimum, maximum);

	for(size_t i = 0; i < ref_values.size(); i++)
		th1d->Fill(BCMath::LogGaus(fit_func(mode, i), ref_values[i], sqrt(ref_values_errors[i] * ref_values_errors[i] + fit_err_func(mode, i))));

	th1d->Draw();
	gPad->Update();
	tc->Print("binlogprobdistribution.pdf");
	tc->Clear();
}

void binomial::chi2plot(){

	std::vector<double> mode = FindMode(GetBestFitParameters());
	
	double maximum = 0;
	for(size_t i = 0; i < ref_values.size(); i++)
		if(BCMath::LogGaus(fit_func(mode, i), ref_values[i], sqrt(ref_values_errors[i] * ref_values_errors[i] + fit_err_func(mode, i))) > maximum)
			maximum = BCMath::LogGaus(fit_func(mode, i), ref_values[i], sqrt(ref_values_errors[i] * ref_values_errors[i] + fit_err_func(mode, i)));


	TCanvas* tc;
	TH1D* th1d;
	

	tc = new TCanvas("chi2plot.pdf");
	tc->cd();
	th1d = new TH1D("chi2plot.pdf", "chi2plot.pdf", (int) (1.5 * sqrt(ref_values.size())), 0, maximum);

	for(size_t i = 0; i < ref_values.size(); i++)
		th1d->Fill((fit_func(mode, i) -  ref_values[i]) * (fit_func(mode, i) -  ref_values[i]) / (ref_values_errors[i] * ref_values_errors[i] + fit_err_func(mode, i)));

	
	th1d->Draw();
	
	gPad->Update();
	tc->Print("chi2plot.pdf");
	tc->Clear();
}

void binomial::chainhistoryplot(){

	// print the evolution of each MCMC
	TCanvas* tc;
	TGraph* tg;

	for(size_t i = 0; i < GetNChains(); i++)
		for(size_t j = 0; j < var_names.size(); j++)
		{
			tc = new TCanvas((var_names[j] + "_" + std::to_string(i) + "_history.pdf").c_str());
			tc->cd();
			tg = new TGraph(MCMC_history.size());
			tg->SetMarkerSize(0);
			for(size_t k = 0; k < MCMC_history.size(); k++)
				tg->SetPoint(k, MCMC_history[k][i][j], k);
			tg->Draw();
			gPad->Update();
			tc->Print((var_names[j] + "_" + std::to_string(i) + "_history.pdf").c_str());
			tc->Clear();
		}

}

double chi2function(double x, int ndof){

	return pow(x, (ndof / 2) - 1) * exp(-x / 2) / (exp2(ndof / 2) * tgamma(ndof / 2));

}

void binomial::chi2distribution(){

	TCanvas* tc = new TCanvas("chi2distribution.pdf");
	tc->cd();

	double sum = 0;
	std::vector<double> mode = FindMode(GetBestFitParameters());

	for(size_t i = 0; i < ref_values.size(); i++)
		sum += (fit_func(mode, i) -  ref_values[i]) * (fit_func(mode, i) -  ref_values[i]) / (ref_values_errors[i] * ref_values_errors[i] + fit_err_func(mode, i));

	int ndof = ref_values.size() - GetNParameters();
	
	TF1* tf1 = new TF1("chi2","pow(x, ([0] / 2) - 1) * exp(-x / 2) / (exp2([0] / 2) * tgamma([0] / 2))", 0, sum);
	tf1->SetParameter(0,ndof);
	tf1->Draw();

	gPad->Update();
	tc->Print("chi2distribution.pdf");
	tc->Clear();
	
}


//double binomial::LogLikelihood(const std::vector<double>& parameters)
//{
//	double logprob = 0;
	
     	//Gaussian distribution
	
	//std::clock_t begin, end;
//	#pragma omp parallel for reduction(+:logprob)  schedule(dynamic) //private(begin,end) num_threads(2)
//	for(size_t i = 0; i < ref_values.size() / 8; i = i+8)
//	{

		//begin = clock();
		
//		logprob += weights[i] * BCMath::LogGaus(fit_func(parameters, i), ref_values[i], sqrt(ref_values_errors[i] * ref_values_errors[i] + fit_err_func(parameters, i))) 
//			+ weights[i + 1] * BCMath::LogGaus(fit_func(parameters, i + 1), ref_values[i + 1], sqrt(ref_values_errors[i + 1] * ref_values_errors[i + 1] + fit_err_func(parameters, i + 1)))
//			+ weights[i + 2] * BCMath::LogGaus(fit_func(parameters, i + 2), ref_values[i + 2], sqrt(ref_values_errors[i + 2] * ref_values_errors[i + 2] + fit_err_func(parameters, i + 2)))
//			+ weights[i + 3] * BCMath::LogGaus(fit_func(parameters, i + 3), ref_values[i + 3], sqrt(ref_values_errors[i + 3] * ref_values_errors[i + 3] + fit_err_func(parameters, i + 3)))
//			+ weights[i + 4] * BCMath::LogGaus(fit_func(parameters, i + 4), ref_values[i + 4], sqrt(ref_values_errors[i + 4] * ref_values_errors[i + 4] + fit_err_func(parameters, i + 4)))
//			+ weights[i + 5] * BCMath::LogGaus(fit_func(parameters, i + 5), ref_values[i + 5], sqrt(ref_values_errors[i + 5] * ref_values_errors[i + 5] + fit_err_func(parameters, i + 5)))
//			+ weights[i + 6] * BCMath::LogGaus(fit_func(parameters, i + 6), ref_values[i + 6], sqrt(ref_values_errors[i + 6] * ref_values_errors[i + 6] + fit_err_func(parameters, i + 6)))
//			+ weights[i + 7] * BCMath::LogGaus(fit_func(parameters, i + 7), ref_values[i + 7], sqrt(ref_values_errors[i + 7] * ref_values_errors[i + 7] + fit_err_func(parameters, i + 7))); 
		//end = clock();
		//std::cout << begin << "\t" << end << "\t" << end - begin << "\t" << "\t" << (end - begin) / 8 << "\t" << omp_get_thread_num() << std::endl;//omp_get_thread_num() << "---";
//	}

//	logprob += weights[1112] * BCMath::LogGaus(fit_func(parameters, 1112), ref_values[1112], sqrt(ref_values_errors[1112] * ref_values_errors[1112] + fit_err_func(parameters, 1112))); 
//	logprob += weights[1113] * BCMath::LogGaus(fit_func(parameters, 1113), ref_values[1113], sqrt(ref_values_errors[1113] * ref_values_errors[1113] + fit_err_func(parameters, 1113))); 
//	logprob += weights[1114] * BCMath::LogGaus(fit_func(parameters, 1114), ref_values[1114], sqrt(ref_values_errors[1114] * ref_values_errors[1114] + fit_err_func(parameters, 1114))); 

	
	// return the log of the conditional probability p(data|pars).
//    	return logprob;
//}

// ---------------------------------------------------------
//double binomial::LogAPrioriProbability(const std::vector<double>& parameters)
// {  
     // return the log of the prior probability p(pars)
     // If you use built-in priors, leave this function commented out.
//	return flat_prior(parameters, var_minvalue, var_maxvalue, var_names);
// }



// ---------------------------------------------------------
// double binomial::CalculateObservables(const std::vector<double>& pars)
// {
//     // Calculate and store obvserables. For example:
//     GetObservable(0) = pow(pars[0], 2);
// }
