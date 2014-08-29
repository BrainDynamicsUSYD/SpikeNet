// header guard at start of header file
#ifndef NEURONNETWORK_H
#define NEURONNETWORK_H

#include <vector>
#include "Neurons.h"
#include "ChemicalSynapses.h" 
//#include "ElectricalSynapses.h"
using namespace std;

class NeuronNetwork{
public:
	NeuronNetwork();
	NeuronNetwork(vector<int> N_array, double dt, int step_tot); // parameterised constructor for heterogeneous coupling
	

	
	vector<int> N_array; // number of neurons in each pupolation
	double dt; // (ms)
	int step_tot; // total number of simulation steps
	int Num_pop; // number of populations

	vector<Neurons> NeuronPopArray; // array of neuron populations
	vector<ChemicalSynapses> ChemicalSynapsesArray; // array of chemical synapses: inter-/intra-population connections
	bool runaway_killed;
	int step_killed; // initialised as -1
	
	void update(int step_current); // update the network to current time step, use "virtual" if want override by derived class


};

inline NeuronNetwork::NeuronNetwork(){};
#endif // End guard at bottom of header file

