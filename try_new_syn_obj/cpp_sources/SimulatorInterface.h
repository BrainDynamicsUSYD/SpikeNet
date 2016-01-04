#ifndef SIMULATORINTERFACE_H
#define SIMULATORINTERFACE_H


#include <iostream> // ofstream: Stream class to write on files, ifstream : Stream class to read from files, istringstream is for input, ostringstream for output
#include <fstream> // fstream : Stream class to both read and write from / to files
#include <sstream>  // stringstream is input and output
#include <string> // "" for string, '' for char


#include <vector>

#include "NeuronNetwork.h"

using namespace std;

class SimulatorInterface{
public:
	SimulatorInterface();
	NeuronNetwork network; // use container?
	
	// Import Network Setup Data
	string in_filename; // path+name
	ifstream inputfile; // current input file (.ygin or .ygin_syn)
	bool import(string in_filename);
	void simulate();
	
	// output data
	string gen_out_filename(); // generate unique file name using time stamp
	string out_filename;
	ofstream output_file;
	void output_results();

	// Format
	string output_suffix; // output filename extension (.ygout)
	char delim; // delim used to delimit the entries in the same line in files, note that the last entry of each line also has a delim
	char indicator; // indicator of data-info line, always the first char in a line, followed by infomation about following data, say, name of the data variable, population index, etc
	char commentor; // indicator of comment lines
	
	// Helper functions
	template < typename Type > Type read_next_entry(istringstream &line_ss);
	template < typename T, typename A > void read_next_line_as_vector( vector<T,A> &vec );
	void write2file(vector<int>& v);

};

#endif
