#ifndef SIMUINTERFACE_H
#define SIMUINTERFACE_H

#include <sstream>  // stringstream is input and output


#include "NeuroNet.h"

using namespace std;

class SimuInterface{
public:
	SimuInterface();
	NeuroNet network; // use container?
	
	// Format
	char delim; // delim used to delimit the entries in the same line in files, note that the last entry of each line also has a delim
	char indicator; // indicator of data-info line, always the first char in a line, followed by infomation about following data, say, name of the data variable, population index, etc
	char commentor; // indicator of comment lines
	
	// Import Network Setup Data
	string in_filename; // path+name
	ifstream inputfile; // current input file (.ygin or .ygin_syn)
	
	bool import(string in_filename);
	void simulate();
	
	// output data
	string output_suffix; // output filename extension (.ygout)
	string out_filename; // without suffix
	string gen_out_filename(); // generate unique file name using time stamp

	void output_results();
#ifdef HDF5
	void output_results_HDF5();
	bool import_HDF5(string in_filename);
	void read_vector_HDF5(const H5File & file, const string & name, vector<int> & v_tmp);
	void read_vector_HDF5(const H5File & file, const string & name, vector<bool> & v_tmp);
	void read_vector_HDF5(const H5File & file, const string & name, vector<double> & v_tmp);
	bool group_exist_HDF5(const H5File & file, const string & name);
	bool group_exist_HDF5(const string & filename, const string & name);
	//void read_matrix_HDF5(const H5File & file, const string & name, vector< vector <double> > & m_tmp);
	template < typename Type > Type read_scalar_HDF5(const H5File & file, const string & name);
#endif
	// Helper functions
	template < typename Type > Type read_next_entry(istringstream &line_ss);
	template < typename T, typename A > void read_next_line_as_vector( vector<T,A> &vec );

};

#endif
