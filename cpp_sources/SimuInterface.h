#ifndef SIMUINTERFACE_H
#define SIMUINTERFACE_H

#include <sstream>  // stringstream is input and output


#include "NeuroNet.h"

using namespace std;

class SimuInterface{
public:
	SimuInterface();
	NeuroNet network; // use container?
	
	// Import Network Setup Data
	string in_filename; // path+name
	ifstream inputfile; // current input file (.ygin or .ygin_syn)

	void simulate();
	
	// output data
	string output_suffix; // output filename extension (.ygout)
	string out_filename; // without suffix
	string gen_out_filename(); // generate unique file name using time stamp


	bool import_restart_HDF5(string in_filename_input);
	void export_restart_HDF5();
	// string gen_restart_filename();
	
	void output_results_HDF5();
	bool import_HDF5(string in_filename);

};

#endif
