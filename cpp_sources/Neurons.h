#ifndef NEURONS_H
#define NEURONS_H
#include <vector>
#include <random> // <boost/random.hpp>
#include <string> // "" for string, '' for char

#include <iostream> // cout/cin, ofstream: Stream class to write on files, ifstream : Stream class to read from files, istringstream is for input, ostringstream for output

class NeuronNetwork; //forward declaration, better than #include "NeuronNetwork.h", if you do not need to access the internal of the class
class ChemicalSynapses;
// class ElectricalSynapses;
class SimulatorInterface;
using namespace std;

// Leaky integrate-and-fire point neuron populations
// Two types of synapses are supported:
// 1) chemical synapses, both excitatory and inhibitory, conductance-based, with transmission delay, post-synaptic conductance described by alpha-function
// 2) electrical synapses (the most common type is modelled, i.e., gap junction)
class Neurons
{
public:
	Neurons(); //default constructor
	Neurons(int pop_ind, int N_input, double dt_input, int step_tot); // parameterised constructor
	friend class NeuronNetwork; // Let NeuronNetwork access its private members
	friend class ChemicalSynapses; // considers ChemicalSynapses as its friend so that ChemicalSynapses can access its private members.
	//friend class ElectricalSynapses;
	friend class SimulatorInterface;

	void init(); // initialise neurons, called by constructor after parameter assignment
	void set_para(string para, char delim); // set parameters if not using default ones


	string dump_para(char delim); // dump all the parameter values used
	void output_results(ofstream& output_file, char delim, char indicator);
	void write2file(ofstream& output_file, char delim, vector< vector<int> >& v);
	void write2file(ofstream& output_file, char delim, vector< vector<double> >& v);
	void write2file(ofstream& output_file, char delim, vector<int>& v);
	void write2file(ofstream& output_file, char delim, vector<double>& v);
	
	void start_V_mean_std_record();
	void record_V_mean_std(); //
	

	void random_V(double firing_probability); // Generate random initial condition for V



	void update_spikes(int step_current); // Find the firing neurons, record them, reset their potential and update nonref
	// Following member(s) should not be inherited
private:
	void update_V(int step_current); // Update potential

public:
	void set_gaussian_I_ext(double mean, double std);
	
	void add_sampling(vector<int> sample_neurons, vector<bool> sample_type, vector<bool> sample_time_points); 
	

	void sample_data(int step_current);

	void init_runaway_killer(double min_ms, double Hz, double Hz_ms); // kill the simulation when runaway activity of the network is detected: 
	// mean number of refractory neurons over previous steps "runaway_steps" in any population exceeding "mean_num_ref"
	void runaway_check(int step_current);

protected:
	// Space and time
	int // actually we can use "unsigned int" here
		pop_ind,
		N; // total number of neurons in the population
	double
		dt; // simulation time step (ms)
	int
		step_tot; // total number of simulation steps

	// Intrinsic neuron properties of the population
	double
		tau_ref;     // Refractory time (ms)
	double
		Cm, // membrane capacitance (uF=1000nF)
		// Potential constants (mV)
		V_rt, // Reset, usually same as leak reversal, use -75 to model relative refractory period??
		V_lk, // Leak reversal
		V_th, //  Threshold
		// Leak conductance 
		g_lk; // (nS)



	// Bookkeeping
	vector<double>
		V,   // Membrane potential
		I_leak, // leaky current
		I_AMPA, // current due to AMPA chemical synapses
		I_GABA, // current due to AMPA chemical synapses
		I_NMDA, // current due to AMPA chemical synapses
		I_GJ, // current due to gap junction (GJ)
		I_ext; // external input current (usually noise)
	int
		ref_steps; // number of simulation steps that a neuron remains refractory after firing
	vector<int>
		ref_step_left; // neurons' state
		// ref_left = 0 for non-refractory, ref_left > 0 for time steps left in refraction 
	vector<int>
		spikes_current; // index vector of current spiking neurons
	vector<int>
		spike_hist_tot, // entire history of spikes from all populations, only recorded when necessary; 
		// a compact data structure, requiring vector "num_spikes_pop" to unpack it.
		num_spikes_pop, // number of spikes at each time step
		num_ref_pop; // number of refractory neurons at each time step

	bool
		V_mean_std_record;
	vector<double>
		V_mean,
		V_std;
		


	// parameters for Generate Gaussian random external current
	double 
		I_ext_mean,
		I_ext_std;

		
		
		
	// Data sampling
	vector<int> 
		sample_neurons; // neuron indices
	vector<bool> 
		sample_type;// boolean vector indicating which data to be sampled
	// must correspond to [V,I_leak,I_AMPA,I_GABA,I_NMDA,I_GJ,I_ext]
	vector<bool> 
		sample_time_points; // logical vector as long as time vector

	vector< vector< vector<double> > >
		sample; // types of data x sampled neurons x time points



	

	// random number generator
	int 
		my_seed;
	typedef default_random_engine
		base_generator_type; // A typedef is used so that base generator type can be changed
      	base_generator_type 
		gen;
	

	// runaway-killer parameters
	bool killer_license; // you need a license to kill
	bool runaway_killed; // true for killed
	int step_killed; // if killing, record when
	int Hz_steps; // the number of steps used to calculate population firing rate (Hz)
	double runaway_Hz; // if the population firing rate is higher than this, kill the simulation
	int min_steps; // minimum number of steps the simulation should run before killed
	int min_pop_size; // No women, no kids;
	
	



}; //class declaration must end with a semi-colon.



inline Neurons::Neurons(){}; // default constructor

#endif
