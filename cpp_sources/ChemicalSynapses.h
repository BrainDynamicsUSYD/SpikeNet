#ifndef CHEMICALSYNAPSES_H
#define CHEMICALSYNAPSES_H
#include <vector>
#include <functional> // pass function as parameter
#include <string> 
#include <iostream> 

#include "Neurons.h"
class Neurons; // #include ".h" for accessing its members, forward declaration for "syntax error: identifier xx" (why both are needed??)

using namespace std;

// Excitatory and inhibitory chemical synapses with transmission delay
class ChemicalSynapses{
public:
	ChemicalSynapses(); // default constructor
	ChemicalSynapses(double dt, int step_tot, char delim, char indicator); // parameterised constructor
	friend class NeuronNetwork; // Let NeuronNetwork access its private members
	friend class SimulatorInterface;

	void init(int synapses_type, int pop_ind_pre, int pop_ind_post, int N_pre, int N_post, vector<int> &C_i, vector<int> &C_j, vector<double> &K_ij, vector<double> &D_ij); // initialise chemical synapses by reading already prepared connections

	void init(int synapses_type, int pop_ind_post, int N_pre, double K_ext, int Num_ext, vector<double> &rate_ext_t, int ia, int ib); // initialise chemical synapses for simulating external Poissonian neuron population
	// [ia,ib] specifies the neuron index range in post population to receive the external stimulation 

	char delim;
	char indicator;
		

	void init();
	void set_para(string para_str);
	string dump_para(); // dump all the parameter values used
	void output_results(ofstream& output_file);

	void recv_pop_data(vector<Neurons> &NeuronPopArray);
	void update(int step_current);
	void send_pop_data(vector<Neurons> &NeuronPopArray);

	void add_short_term_depression(int STD_on_step);
	void add_inh_STDP(int inh_STDP_on_step);
	void add_sampling(vector<int> sample_neurons, vector<bool> sample_time_points); 
	void sample_data(int step_current);


	void write2file(ofstream& output_file, vector< vector<int> >& v);
	void write2file(ofstream& output_file, vector< vector<double> >& v);
	void write2file(ofstream& output_file, vector<int>& v);
	void write2file(ofstream& output_file, vector<double>& v);
	
	void start_stats_record();
	void record_stats(); //
	
	
protected:
	// constants
	double
		dt; // (msec) simulation time step
	int
		step_tot,
		pop_ind_pre,
		pop_ind_post,
		N_pre, // pre-synaptic population size
		N_post;
	int
		synapses_type; // 0 = AMPA, 1 = GABA, 2 = NMDA
	double
		V_ex, // Excitatory reversal
		V_in; // Inhibitory reversal
	int
		max_delay_steps;
		
	// A copy of data from pre-synaptic population
	vector<double> // This is problematic!!!
		*V_post; // from post-synaptic population
	vector<int>
		*spikes_pre,
		*spikes_post; // current spikes from pre-synaptic population


	// currents into post-synaptic population
	vector<double>
		I; 
	//
	bool
		stats_record;
	vector<double>
		I_mean,
		I_std;


	// Data sampling
	vector<int> 
		sample_neurons; // post-synaptic neuron indices
	vector<bool> 
		sample_time_points; // logical vector as long as time vector
	vector< vector<double> >
		sample; //  sampled neurons x time points

	
	


	// Build-in paramters for time-evolution of post-synaptic conductance change
	// 1-variable "s(t)" kinetic synapses model
	double 
		Dt_trans_AMPA, // msec, duration of transmitter release pulse (square-shape) activated by spike
		Dt_trans_GABA,
		Dt_trans_NMDA; // Note that here "transmitter" actually means "the effect of transmitter on gating variable"

	double	
		tau_decay_AMPA,
		tau_decay_GABA,
		tau_decay_NMDA;
	double
		tau_decay; // msec, decay time
	int
		steps_trans; // tranmitter duration in simulation steps
	vector<double>
		K_trans; // 1.0/transmitter_steps!
	double
		exp_step; // exp(-dt/tau)

	// voltage-dependent part B(V) (look-up table):
	double
		miuMg_NMDA, // mM^-1, concentration of [Mg2+] is around 1 mM
		gamma_NMDA, // mV^-1
		B_V_min, // < V_in = -80
		B_V_max, // > V_th = -55
		B_dV;
	vector<double>
		B; // B = 1 / (1 + miuMg_NMDA*exp(-gamma_NMDA*V))


	// Inhibitory-to-excitatory coupling STDP plasticity as in
	// ref: Inhibitory plasticity balances excitation and inhibition in sensory pathways and memory networks
	bool
		inh_STDP;
	vector<double>
		x_trace_pre,
		x_trace_post;
	double 
		tau_STDP,
		exp_step_STDP, // exp(-dt/tau_STPD)
		eta_STDP, // learning rate
		rho_0_STDP,
		alpha_STDP; // depression factor
	int
		inh_STDP_on_step;
	vector< vector<int> >
		j_2_i, // j_2_i[j_post] gives all the i_pre's (indices of pre-synaptic neurons)
		j_2_syn_ind; // j_2_syn_ind[j_post] gives all the syn_ind's so that K[i_pre][syn_ind] is a synapse onto j_post


	// connection matrices and bookkeeping for 1-variable kinetic synapse model
	double // short-term depression constants
		p_ves, // ves for vesicle
		tau_ves,
		exp_ves;
	bool
		STD; //  
	int
		STD_on_step; // the step where STD should turned on

	vector<double>
		s, // pre-synaptic dynamics
		gs_sum, // post-synaptic dynamics
		f_ves; // short-term depression: the fraction of available vesicles
	vector<int>
		trans_left; // 
	vector< vector<double> >
		d_gs_sum_buffer; // d_gs_sum_buffer[time index][post-synaptic neuron index], ring buffer
	int
		buffer_steps;
	
	vector< vector<int> >
		C, // connection index
		   // each entry in the C matrix is the index of a POST-SYNAPTIC neuron (pre to post)
		D; // connection delay (in simulation time steps)
	vector< vector<double> >
		K; // connection strength, measuring the strength of synaptic conductance between two cells

	vector< vector<double> > 
		tmp_data; // temporary data container for debugging


	// Simulating external Poisson population (noise)
	double 
		K_ext; // identical connection strength for external pre-synaptic neurons (chemical synapses)
	int 
		Num_ext, // number of external pre-synaptic neurons (chemical synapses) per post-synaptic neuron
		ia, 
		ib; // range of neuron index in post population [ia,ib]
		
	vector<double> 
		rate_ext_t; // identical rate of firing for external pre-synaptic neurons (chemical synapses)

	// Random number generator
	int 
		my_seed;
	typedef mt19937 
		base_generator_type; // A typedef is used so that base generator type can be changed
	base_generator_type 
		gen;

};

inline ChemicalSynapses::ChemicalSynapses(){};

#endif
