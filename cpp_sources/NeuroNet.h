// header guard at start of header file
#ifndef NEURONET_H
#define NEURONET_H

#include "ChemSyn.h" 


using namespace std;

class NeuroNet{
public:
	NeuroNet();
	NeuroNet(vector<int> N_array, double dt, int step_tot, char delim, char indicator); // parameterised constructor for heterogeneous coupling

	vector<int> N_array; // number of neurons in each pupolation
	double dt; // (ms)
	int step_tot; // total number of simulation steps
	int Num_pop; // number of populations

	vector<NeuroPop*> NeuroPopArray; // array of neuron populations
	vector<ChemSyn*> ChemSynArray; // array of chemical synapses: inter-/intra-population connections
	bool runaway_killed;
	int step_killed; // initialised as -1
	
	void update(int step_current); // update the network to current time step, use "virtual" if want override by derived class
	
	char delim;
	char indicator;
	void output_results(ofstream& output_file);
	void output_results(H5File& file_HDF5);
	


};

inline NeuroNet::NeuroNet(){};
#endif // End guard at bottom of header file

