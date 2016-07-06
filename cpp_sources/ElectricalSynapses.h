#ifndef ELECTRICALSYNAPSES_H
#define ELECTRICALSYNAPSES_H

//#include <vector>

#include "Neurons.h"

using namespace std;

/*CurrentLy tranmission delay is not modelled since that will require us to record membrane potential history.
Also, the connection strength is assumed to be the same for all the electrical couplings*/
class ElectricalSynapses{
public:
	ElectricalSynapses();
	ElectricalSynapses(double dt);

	void update_action_potential(Neurons &pop, int step_current);
	void update_electrical_coupling(Neurons &pre_pop, Neurons &post_pop);

public:
		
	double
		dt, // simulation time step
		extent; // radius of the homogeneous cicular coupling topology

	// square-shape spike
	// impluse-like shape to model the contribution to post-synaptic excitatory current, neglecting the actual shape of the action potential
	double
		// Gap junction conductance
		g_GJ; // no good experimental data for reference
		// compared to Leak conductance divided by capacitance (ms^-1) g_lk = 0.05, or g_lk^-1 = 20 ms, eqivalent to 50e-6 Ohm^-1/cm^2 (50e-6 S/cm^2) with C=10^-6 F/cm^2 ?????
	double
		spike_peak, // (mV), height of spike, artificial value only to model EPSC/IPSC
		spike_duration; // (ms), duration of spike
	int
		spike_steps; // spike witdth in number of simulation steps
	int
		boundary; // 1 = "Periodic", 2 = "Free"
	vector<vector<int>>
		C; // post synaptic connection matrix of gap junction


};

inline ElectricalSynapses::ElectricalSynapses(){};

#endif
