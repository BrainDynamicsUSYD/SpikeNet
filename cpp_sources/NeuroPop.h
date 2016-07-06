#ifndef NEUROPOP_H
#define NEUROPOP_H
#include <vector>
#include <random> // <boost/random.hpp>
#include <string> // "" for string, '' for char
#include <iostream> // cout/cin, ofstream: Stream class to write on files, ifstream : Stream class to read from files, istringstream is for input, ostringstream for output
#include <fstream> // fstream : Stream class to both read and write from / to files
class NeuronNetwork; //forward declaration, better than #include "NeuronNetwork.h", if you do not need to access the internal of the class
class ChemSyn;
// class ElectricalSynapses;
class SimulatorInterface;
using namespace std;

// Leaky integrate-and-fire point neuron populations
// Two types of synapses are supported:
// 1) chemical synapses, both excitatory and inhibitory, conductance-based, with transmission delay, post-synaptic conductance described by alpha-function
// 2) electrical synapses (the most common type is modelled, i.e., gap junction)
class NeuroPop{
public:
	NeuroPop(); /// default constructor
	NeuroPop(int pop_ind, int N_input, double dt_input, int step_tot, char delim, char indicator); /// parameterised constructor

	void init(); // initialise neurons, called by constructor after parameter assignment
	void set_para(string para); /// set parameters if not using default ones

	void output_results(ofstream& output_file);
	void output_sampled_data_real_time(int step_current);

	
	void recv_I(vector<double>& I_add, int pop_ind_pre, int syn_type);
	const vector< int >& get_spikes_current();
	const vector< double >& get_V();
	const bool & get_runaway_killed();
	
	void start_stats_record();

	void start_LFP_record(vector< vector<bool> >& LFP_neurons);
	void record_LFP();

	void random_V(double firing_probability); /// Generate random initial condition for V. This function is deprecated!
	void set_init_condition(double r_V0, double p_fire); /// Uniform random distribution [V_rt, V_rt + (V_th - V_rt)*r_V0] and then randomly set neurons to fire according to p_fire

	void update_spikes(int step_current); /// Find the firing neurons, record them, reset their potential and update nonref
	// Following member(s) should not be inherited
	void update_V(int step_current); // Update potential
	void set_gaussian_I_ext(vector<double>& mean, vector<double>& std);
	void set_gaussian_g_ext(vector<double>& mean, vector<double>& std);
	
	void add_sampling(vector<int>& sample_neurons, vector<bool>& sample_type, vector<bool>& sample_time_points); 
	void add_sampling_real_time(vector<int>& sample_neurons_input, vector<bool>& sample_type_input, vector<bool>& sample_time_points_input, string samp_file_name);

	void init_runaway_killer(double min_ms, double Hz, double Hz_ms); /// kill the simulation when runaway activity of the network is detected: 
	// mean number of refractory neurons over previous steps "runaway_steps" in any population exceeding "mean_num_ref"
	void runaway_check(int step_current);
	void add_perturbation(int step_perturb);
	void add_spike_freq_adpt(); /// add spike-frequency adaptation
	
private:
	void generate_I_ext(int step_current);
	void record_stats(int step_current); 
	void write2file(ofstream& output_file, vector< vector<int> >& v);
	void write2file(ofstream& output_file, vector< vector<double> >& v);
	void write2file(ofstream& output_file, vector<int>& v);
	void write2file(ofstream& output_file, vector<double>& v);
	string dump_para(); /// dump all the parameter values used
	char delim;
	char indicator;
	void sample_data(int step_current);

protected:
	ofstream 
		samp_file; /// the output file stream for sampled time series
	string 
		samp_file_name; /// the file name for sampled time series
	
	// Space and time
	int // actually we can use "unsigned int" here
		pop_ind, /// population index
		N; /// total number of neurons in the population
	double
		dt; /// simulation time step (ms)
	int
		step_tot; /// total number of simulation steps

	// Intrinsic neuron properties of the population
	double
		tau_ref;  /// absolute refractory time (ms)
	double
		Cm, /// membrane capacitance (uF=1000nF)
		// Potential constants (mV)
		V_rt, /// reset potential
		V_lk, /// leaky reversal
		V_th, ///  firing threshold
		// Leak conductance 
		g_lk; /// leaky conductance (nS)



	// Bookkeeping
	vector<double>
		V,   /// neuron membrane potential
		I_leak, /// leaky current
		I_input, /// sum of all input currents (I_leak not included)
		I_AMPA, /// current due to AMPA chemical synapses
		I_GABA, /// current due to AMPA chemical synapses
		I_NMDA, /// current due to AMPA chemical synapses
		I_GJ, /// current due to gap junction (GJ)
		I_ext; /// external input current (usually noise)
	int
		ref_steps; /// number of simulation steps that a neuron remains refractory after firing
	vector<int>
		ref_step_left; /// current number of refractory steps left for the neurons
		// ref_left = 0 for non-refractory, ref_left > 0 for time steps left in refraction 
	vector<int>
		spikes_current; /// index vector of current spiking neurons
	vector<int>
		spike_hist_tot, /// entire history of spikes of the population packed into one vector; 
						/// requires num_spikes_pop to unpack it
		num_spikes_pop, /// number of spikes at each time step
		num_ref_pop; /// number of refractory neurons at each time step

	bool
		stats_record; /// whether stats should be recorded (false by default)
	vector<double>
		V_mean, /// mean of membrane potential averaged over neurons at each time step
		V_std, /// std of membrane potential averaged over neurons at each time step
		I_input_mean, /// mean of I_input averaged over neurons at each time step
		I_input_std, /// std of I_input averaged over neurons at each time step
		I_AMPA_acc, /// accumulator for AMPA input currents into each neuron
		I_AMPA_time_avg, /// I_AMPA averaged over time for each neuron
		I_NMDA_acc,  /// accumulator for NMDA input currents into each neuron
		I_NMDA_time_avg, /// I_NMDA averaged over time for each neuron
		I_GABA_acc, /// accumulator for GABA input currents into each neuron
		I_GABA_time_avg, /// I_GABA averaged over time for each neuron
		IE_ratio; /// I-E ratio for each neuron

	bool
		LFP_record; /// whether LFP should be recorded (false by default)
	vector< vector<bool> >
		LFP_neurons; /// each component vector defines a LFP measure by specifying which neurons should be included
	vector< vector<double> >
		LFP; /// each component vector is a LFP time series
	
	// parameters for Generate Gaussian random external current
	vector<double> // a vector for each neuron
		I_ext_mean, /// mean of external currents (Gaussian noise) for each neuron
		I_ext_std; /// std of external currents (Gaussian noise) for each neuron
	
	// parameters for Generate Gaussian random external conductance
	vector<double> // a vector for each neuron
		g_ext_mean, /// mean of external conductance (Gaussian noise) for each neuron
		g_ext_std; /// std of external conductance (Gaussian noise) for each neuron
	double 
		V_ext; // reversal potential for external conductance (Gaussian noise)
			
	// spike-frequency adaptation
	bool 
		spike_freq_adpt; /// whether spike-frequency adaptation should be used (false by default)
	vector<double>
		g_K, /// potassium conductance that produces spike-frequency adaptation (nS)
	    I_K; /// potassium currents that produces spike-frequency adaptation (nS)
	double
		V_K, /// reversal potential for the potassium conductance 
		dg_K, 
		tau_K, 
		exp_K_step; 
		
	// Data sampling
	vector<int> 
		sample_neurons; /// indices of the neurons to be sampled
	vector<bool> 
		sample_type; /// specifies which data to be sampled;
					 /// must correspond to [V,I_leak,I_AMPA,I_GABA,I_NMDA,I_GJ,I_ext, I_K]
	vector<bool> 
		sample_time_points; /// specifies which time steps to be sampled;
	vector< vector< vector<double> > >
		sample; /// sampled data with a dimension (types of data) by (sampled neurons) by (time points)


	// perturbation
	int
		step_perturb, /// the step where the perturbation takes place (removal of one spike)
		spike_removed; /// the index of the neuron of which the spike should be removed;
	

	// random number generator
	int 
		my_seed; /// the seed for the random number generators used for this population
	typedef default_random_engine
		base_generator_type; /// a typedef is used so that base generator type can be changed
    base_generator_type 
		gen; /// the random number generator
	

	// runaway-killer parameters
	bool killer_license; /// you need a license to kill (false by default)
	bool runaway_killed; /// whether the simulation has been killed
	int step_killed; /// if killed, record when
	int Hz_steps; /// the number of steps to be averaged over to calculate population firing rate (Hz)
	double runaway_Hz; /// if the population firing rate is higher than this, kill the simulation
	int min_steps; /// minimum number of steps the simulation should run before being killed
	int min_pop_size; /// No women, no kids;

}; //class declaration must end with a semi-colon.


inline NeuroPop::NeuroPop(){}; /// default constructor

#endif
