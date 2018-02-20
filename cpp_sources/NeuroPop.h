#ifndef NEUROPOP_H
#define NEUROPOP_H
#include <vector>
#include <random> // <boost/random.hpp>
#include <string.h> // 
#include <string> //
#include <iostream> // cout/cin, ofstream: Stream class to write on files, ifstream : Stream class to read from files, istringstream is for input, ostringstream for output
#include <fstream> // fstream : Stream class to both read and write from / to files

#include "MyIO.h"
#include "ChemSyn.h"

#ifdef HDF5
#include <H5Cpp.h>
#include <hdf5_hl.h>
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif
#endif

class NeuronNetwork; //forward declaration, better than #include "NeuronNetwork.h", if you do not need to access the internal of the class
class ChemSyn;
// class ElectricalSynapses;
class SimulatorInterface;
using namespace std;

// Leaky integrate-and-fire point neuron populations
// Two types of synapses are supported:
// 1) chemical synapses, both excitatory and inhibitory, conductance-based, with transmission delay, post-synaptic conductance described by alpha-function
// 2) electrical synapses (the most common type is modelled, i.e., gap junction)
class NeuroPop {
public:
	NeuroPop(); /// default constructor
	NeuroPop(const int pop_ind, const int N_input, const double dt_input, const int step_tot); /// parameterised constructor

	void init(); // initialise neurons, called by constructor after parameter assignment
	void set_para(string para); /// set parameters if not using default ones
	void set_seed(int seed); /// manually set RNG seed

	void recv_I(vector<double>& I_add, const int pop_ind_pre, const int syn_type);
	const vector< int >& get_spikes_current();
	const vector< double >& get_V();
	const vector< int >& get_ref_step_left();
	const bool & get_runaway_killed();
	const double & get_Cm();

	void start_stats_record();
	void start_cov_record(const int time_start, const int time_end);
	
	void start_LFP_record(const vector< vector<double> >& LFP_neurons);

	void random_V(const double firing_probability); /// Generate random initial condition for V. This function is deprecated!

	/// uniform random distribution [V_rt, V_rt + (V_th - V_rt)*r_V0] and then randomly set neurons to fire according to p_fire
	void set_init_condition(const double r_V0, const double p_fire); 
	/// set initial V according to externally given vector
	void set_init_V_external(const vector<double>& external_init_V);
	
	
	void set_neuron_model(int n_mod);
	void set_ELIF_Params(double elif_delT, double elif_VT);
	void init_poisson_pop(double rate);

	void update_spikes(const int step_current); /// Find the firing neurons, record them, reset their potential and update nonref
	// Following member(s) should not be inherited
	void update_V(const int step_current); // Update potential
	void set_gaussian_I_ext(const vector<double>& mean, const vector<double>& std);
	void set_gaussian_g_ext(const vector<double>& mean, const vector<double>& std);

	void load_file_spike_input(string fname);
	void load_file_current_input(string fname);
	void get_current_from_file();


	void add_JH_Learn(double noise);
	void get_all_rhat_JHLearn(vector<ChemSyn*> &ChemSynArray, const int step_current);
	void reset_rhat();
	const vector< int > & get_spikes_noise(); 
	const double & get_noise(); 
	
	void import_restart(H5File & file, int pop_ind, string out_filename);
	void export_restart(Group & group);
	void output_results(H5File& file_HDF5);
	void add_sampling_real_time_HDF5(const vector<int>& sample_neurons_input, const vector<bool>& sample_type_input, const vector<bool>& sample_time_points_input, string samp_file_name);
	void add_sampling(const vector<int>& sample_neurons_input, const vector<bool>& sample_type_input, const vector<bool>& sample_time_points_input);
	void init_runaway_killer(const double min_ms, const double Hz, const double Hz_ms); /// kill the simulation when runaway activity of the network is detected:
	// mean number of refractory neurons over previous steps "runaway_steps" in any population exceeding "mean_num_ref"

	void add_perturbation(const int step_perturb);
	void add_spike_freq_adpt(); /// add spike-frequency adaptation	
	void set_spike_freq_adpt_para(const double dg_K_input); /// add spike-frequency adaptation
	void set_spike_freq_adpt_para_heter(const vector<double>& dg_K_heter_input, const int start_step, const int end_step); /// add heterogenous spike-frequency adaptation

private:
	void generate_I_ext();
	void record_stats(const int step_current);
	void output_sampled_data_real_time_HDF5(const int step_current);
	string dump_para(); /// dump all the parameter values used
	void sample_data(const int step_current);
	void runaway_check(const int step_current);
	void record_LFP();

protected:

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
	Cm; /// membrane capacitance (uF=1000nF)
	double
	tau_ref;  /// absolute refractory time (ms)
	double
	// Potential constants (mV)
	V_rt, /// reset potential
	V_lk, /// leaky reversal
	V_th, ///  firing threshold
	// Leak conductance
	g_lk; /// leaky conductance (nS)
	int
	neuron_model = 0; // 0 for LIF, 1 for Exp LIF

	//Exponential Leaky Integrate and Fire neuron model parameters
	struct ELIF {
		double delT = 0;
		double V_T;

	} elif;

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

	vector<int> ref_step_left; /// current number of refractory steps left for the neurons
	/// ref_left = 0 for non-refractory, ref_left > 0 for time steps left in refraction
	vector<int>
	spikes_current; /// index vector of current spiking neurons
	vector<int>
	spike_hist_tot, /// entire history of spikes of the population packed into one vector;
	/// requires num_spikes_pop to unpack it
	num_spikes_pop, /// number of spikes at each time step
	num_ref_pop; /// number of refractory neurons at each time step

	struct Stats {
		bool
			record, /// whether stats should be recorded (false by default)
			record_cov;
			vector<double>
			V_mean, /// mean of membrane potential averaged over neurons at each time step
			V_std, /// std of membrane potential averaged over neurons at each time step
			I_input_mean, /// mean of I_input averaged over neurons at each time step
			I_input_std, /// std of I_input averaged over neurons at each time step
			I_AMPA_time_avg, /// I_AMPA averaged over time for each neuron
			I_NMDA_time_avg, /// I_NMDA averaged over time for each neuron
			I_GABA_time_avg, /// I_GABA averaged over time for each neuron
			I_ext_time_avg, /// I_ext averaged over time for each neuron
			I_tot_time_mean, /// mean of total current for each neuron
			I_tot_time_var, /// variance of total current for each neuron
			V_time_mean, /// mean of membrane potential for each neuron
			V_time_mean_dumb, /// mean of membrane potential for each neuron
			V_time_var,
			IE_ratio; /// I-E ratio for each neuron
			vector< vector <double> >
			V_time_cov; /// covariance of membrane potential for each neuron
			int
			time_start_cov,
			time_end_cov;
		} stats;

		struct Lfp {
			bool
		record; /// whether LFP should be recorded (false by default)
		vector< vector<double> >
		neurons; /// each component vector defines a LFP measure by specifying which neurons should be included
		vector< vector<double> >
		data; /// each component vector is a LFP time series
	} LFP;

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
	struct Spike_freq_adpt_heter
	{
		bool on = false; /// whether heterogenous spike-frequency adaptation should be used (false by default)
		int
		start_step,
		end_step;
	} spike_freq_adpt_heter;
	
	vector<double>
	g_K, /// potassium conductance that produces spike-frequency adaptation (nS)
	I_K, /// potassium currents that produces spike-frequency adaptation (nS)
	dg_K_heter; /// heterogeneous spike-freq-adap
	double
	V_K, /// reversal potential for the potassium conductance
	dg_K,
	tau_K,
	exp_K_step;

	// Data sampling
	struct Sample {
		int file_type; /// 0 for no sample, 1 for text-based, 2 for hdf5-based
		string
		file_name; /// the file name for sampled time series
		H5File *
		file_HDF5;
		DataSet
		V_dataset, I_leak_dataset, I_AMPA_dataset, I_GABA_dataset, I_NMDA_dataset, I_GJ_dataset, I_ext_dataset, I_K_dataset, rhatE_dataset, rhatI_dataset;;
		int
		ctr = 0; /// counter that counts how many steps have been sampled
		vector<int>
		neurons; /// indices of the neurons to be sampled
		vector<bool>
		type; /// specifies which data to be sampled;
		/// must correspond to [V,I_leak,I_AMPA,I_GABA,I_NMDA,I_GJ,I_ext, I_K]
		vector<bool>
		time_points; /// specifies which time steps to be sampled;
		int
		N_steps, /// number of sample steps
		N_neurons; /// number of sample neurons
		vector< vector< vector<double> > >
		data; /// sampled data with a dimension (types of data) by (sampled neurons) by (time points)
	} sample;


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
	struct Runaway_killer {
		bool license; /// you need a license to kill (false by default)
		bool runaway_killed; /// whether the simulation has been killed
		int step_killed; /// if killed, record when
		int Hz_steps; /// the number of steps to be averaged over to calculate population firing rate (Hz)
		double runaway_Hz; /// if the population firing rate is higher than this, kill the simulation
		int min_steps; /// minimum number of steps the simulation should run before being killed
		int min_pop_size; /// No women, no kids;
	} killer;


	// external spike file
	struct Spike_File {
		int on = 0;
		unsigned int spike_ind = 0;
		vector< vector < int > > spikes;
		string file_name;

	} spike_file;

	struct Poisson_Pop{
		int on=0;
		vector <int>  next_spike_time; // timestep of next spike for each neuron
		double rate; //in spikes per timestep
	}poisson_pop;


	// external current file
	struct Current_file {
		int on = 0; /// flag to add the external current file
		int start_step = 0; /// specific step to start adding external current file
		int end_step = 0; /// specific step to remove external current file
		unsigned int current_ind = 0;
		vector< vector < int > > neurons;
		vector< vector < double > > current;
		string file_name;
		unsigned int steps_per_frame; /// simulation timesteps per frame i.e. 1000/fps
		unsigned int framestep = 0; /// framestep is tracking how many timesteps have passed for the current frame
		double mean_curr; /// strength scale of external current from file

	} current_file;

	struct JH_Learn_Pop{
		bool on=false; //indicates if this learning is to be used
		vector<double> rhatI,rhatE;		
		double noise;
		vector<int> noise_spikes;
	} jh_learn_pop;
}; //class declaration must end with a semi-colon.

void Welford_online(const vector<double>& new_data, vector<double>& M, vector<double>& S, const int K, const bool is_end); /// online mean and var calculation: Welford's method (1962, Technometrixcs)
void Welford_online(const vector<double>& new_data, vector<double>& M, vector< vector <double> >& Cov, const int K, const bool is_end); 
void Welford_online(const vector<double>& new_data, vector<double>& M, const int K);
void Welford_online(const vector<double>& data, double& M, double& S);

inline NeuroPop::NeuroPop() {}; /// default constructor

#endif
