#include "datareading.h"

using namespace std;

datareading::datareading(char* interpolation){

	ifstream ifile(interpolation);

	if (ifile.is_open()) 
	{
		cout << "reading interpolation file ...";
		getline(ifile, line);

		//row-wise stepping through the file
		while (getline(ifile, line)) //row-wise filereading
		{
	
			//TODO:if abragen rausziehen
			length = line.length();
			parameter_reading();
			var_min_max_reading();
			fit_params_reading();
			analysis_name();
		}

		ifile.close();
		cout << "complete" << endl;
	}
	else
		cout << "unable to read interpolation file" << endl;
}

void datareading::parameter_reading(){
	//just extracting the line with the parameter names
	//this will be used after the amount of parameters is known
	//this is done due to the formation of the file
	if (line.substr(0, 10) == "ParamNames")
	{
		parameternames = line;
	}

	//reading the dimension and the parameter names
	if (line.substr(0, 9) == "Dimension")
	{
		pos_begin = 10;
		//cout << "Dimension: ";
		pos_end = line.find(" ", pos_begin + 1);

		dimension = atoi(line.substr(pos_begin, pos_end - pos_begin).c_str());

		//cout << dimension << endl;
		//cout << "Parameter names: ";

		pos_begin = 11;
		while(var_names.size() < dimension)
		{
		pos_end = parameternames.find(" ", pos_begin + 1);
		var_names.push_back(parameternames.substr(pos_begin + 1, pos_end - pos_begin));
		//cout << parameternames.substr(pos_begin, pos_end - pos_begin);
		pos_begin = pos_end;
		}
			//cout << endl;
	}

}

void datareading::var_min_max_reading(){

	if (line.substr(0, 12) == "MinParamVals")
	{
		pos_begin = 14;
		//cout << "Minimum values: ";
		while(var_minvalue.size() < dimension)
		{
			pos_end = line.find(" ", pos_begin + 1);
			var_minvalue.push_back(atof(line.substr(pos_begin, pos_end - pos_begin).c_str()));
			//cout << line.substr(pos_begin, pos_end - pos_begin);
			pos_begin = pos_end;
		}

		//cout << endl;
	}

	if (line.substr(0, 12) == "MaxParamVals")
	{
		pos_begin = 14;
		//cout << "Maximum values: ";

		while(var_maxvalue.size() < dimension)
		{
			pos_end = line.find(" ", pos_begin + 1);
			var_maxvalue.push_back(atof(line.substr(pos_begin, pos_end - pos_begin).c_str()));
			//cout << line.substr(pos_begin, pos_end - pos_begin);
			pos_begin = pos_end;
		}

		//cout << endl;
  	}

}


void datareading::fit_params_reading(){

	if (line.substr(2, 3) == "val") //searching for interpolationparameters
	{

		if(fit_order == 0)
		{
		//reading the weightning of the bin
		pos_begin = 8;
		//cout << "bin weighting: ";

		pos_end = line.find(" ", pos_begin + 1);

		fit_order = atoi(line.substr(pos_begin, pos_end - pos_begin).c_str());
		//cout << fit_order << endl;
		}

		//reading the value of every parameter in every bin
		pos_begin = 10;
		//cout << "parameter value: ";

		while (pos_begin + 1 < length) 
		{
			pos_end = line.find(" ", pos_begin + 1);

			//cout << line.substr(pos_begin, pos_end - pos_begin) << " ";
			fit_params.push_back(atof(line.substr(pos_begin, pos_end - pos_begin).c_str()));
			
			pos_begin = pos_end;

		}
		//cout << endl;

		//calculating number of parameters in use
		//it will be done once per file, therefore the 'if'-clause
		if (fit_num_params == 0) 
		{
			
			fit_num_params = fit_params.size();
			//cout << "numer of parameters: " << fit_num_params << endl;
		}

	}

	//reading the binwise error of the interpolation
	if (line.substr(2, 3) == "err")
	{
	
		pos_begin = 10;
		//cout << "error: ";

		//cout << atof(line.substr(pos_begin, pos_end - pos_begin).c_str()) << endl;

		while (pos_begin + 1 < length) 
		{
			pos_end = line.find(" ", pos_begin + 1);

			//cout << line.substr(pos_begin, pos_end - pos_begin) << " ";
			fit_error.push_back(atof(line.substr(pos_begin, pos_end - pos_begin).c_str()));
			
			pos_begin = pos_end;
		}
	}

}

void datareading::analysis_name(){

	if (line.substr(0, 1) == "/")
	{
		pos_begin = 1;
		
		pos_end = line.find("/", pos_begin + 1);
		analysis.push_back(line.substr(pos_begin, pos_end - pos_begin));

		//cout << line.substr(pos_begin, pos_end - pos_begin) << endl;
		pos_begin = pos_end + 1;

		pos_end = line.find("#", pos_begin);
		observable.push_back(line.substr(pos_begin, pos_end - pos_begin));
			
		//cout << line.substr(pos_begin, pos_end - pos_begin) << endl;

		if(observable.size() == 1)
			changing_histo.push_back(0);

		if(observable.size() > 1)
			if(observable.back() != observable[observable.size() - 2])
				changing_histo.push_back(observable.size() - 1);
	}

}

int datareading::get_num_histos(){

	if(analysis.size() == 0)
		return 0;

	string tmp_analysis = analysis[0];
	string tmp_observable = observable[0];
	int num_histos = 1;

	for(size_t i = 1; i < analysis.size(); i++)
	{
		if(analysis[i] != tmp_analysis || observable[i] != tmp_observable)
		{
			tmp_analysis = analysis[i];
			tmp_observable = observable[i];
			num_histos++;
		}
	}

	return num_histos;

}


