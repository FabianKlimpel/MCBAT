#include "refreading.h"

using namespace std;

//Constructor that reads an inputfile and stores the values and errors of the reference data
refreading::refreading(string analysis, string observable){

	//reading the root file and converting it to a data type,
	//that allows reading the properties of the plots
	cout << "initialize reference data reading ...";
	analysis.append(".root");
	TFile tf(analysis.c_str());
	TGraphAsymmErrors* tgae = (TGraphAsymmErrors*) tf.FindObjectAny(observable.c_str());
	cout << "complete" << endl;

	cout << "start reading ...";
	//reading the number of bins and the content of each
	ref_num_bins = tgae->GetN();

	//TODO:check, ob fehler stets symmetrisch sind
	//simple push_back on values and errors
	for(int i = 0; i < ref_num_bins; i++)
	{
		ref_values.push_back(tgae->GetY()[i]);
		ref_values_errors.push_back(tgae->GetErrorY(i));
	}
	cout << "complete" << endl;
	
}

//Constructor that reads multiple inputfiles and stores the values and errors of the reference data
//the names of the needed inputfiles is derived of a datareading object
refreading::refreading(datareading &dr){

	int num_histos = dr.get_num_histos();

	//walk over all histograms named in dr
	for(int i = 0; i < num_histos; i++)
	{

	//reading the root file and converting it to a data type,
	//that allows reading the properties of the plots
	cout << "initialize reference data reading from " << dr.analysis[dr.changing_histo[i]] << " : " << dr.observable[dr.changing_histo[i]] << " ...";
	dr.analysis[dr.changing_histo[i]].append(".root");
	TFile tf(dr.analysis[dr.changing_histo[i]].c_str());
	TGraphAsymmErrors* tgae = (TGraphAsymmErrors*) tf.FindObjectAny(dr.observable[dr.changing_histo[i]].c_str());
	cout << "complete" << endl;

	cout << "start reading ...";
	//reading the number of bins and the content of each
	ref_num_bins = tgae->GetN();

	//TODO:check, ob fehler stets symmetrisch sind
	//simple push_back on values and errors
	for(int j = 0; j < ref_num_bins; j++)
	{
		ref_values.push_back(tgae->GetY()[j]);
		ref_values_errors.push_back(tgae->GetErrorY(j));
	}
	cout << "complete" << endl;
	
	}
	cout << "all reference data loaded" << endl;
}
