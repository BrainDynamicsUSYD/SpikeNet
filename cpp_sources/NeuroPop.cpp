#include <stdlib.h>
#include <iomanip>
#include <cmath>
#include <algorithm>
// for transform
#include <stdio.h> // for printf
#include <time.h>       /* time */
#include <sstream>  // stringstream is input and output
#include "NeuroPop.h"
#include <functional> // for bind(), plus

NeuroPop::NeuroPop(const int pop_ind_input, const int N_input, const double dt_input, const int step_tot_input)
{
	pop_ind = pop_ind_input;
	N = N_input;
	dt = dt_input;
	step_tot = step_tot_input; // this parameter is designed to be self-adapting (step_killed), so should be any other stuff that relies on it!!

	// Using consistant units: msec+mV+nF+miuS+nA
	// Initialise default parameters
	// Model of <<The Asynchronous State in Cortical Circuits>>
	Cm = 0.25; // nF
	tau_ref = 2.0; // absolute refractory time (ms), 3.0
	// Potential constants (mV)
	V_rt = -60.0;   // Reset, usually same as leak reversal, use -75.0 to model relative refractory period??
	V_lk = -70.0;   // Leak reversal, -70.0
	V_th = -50.0;   // Threshold // -55.0
	V_ext = 0.0;    // reversal potential for external currents
	// Leak conductance
	g_lk = 0.0167;   // (uS=miuSiemens), time constants=Cm/gL=15 ms!
	// spike-frequency adaptation parameter
	V_K = -85.0; // mV
	dg_K = 0.01; // (uS=miuSiemens)
	tau_K = 80; // ms
	// Initialise defualt parameters
	init();

}

void NeuroPop::init()
{
	// Initialise arrarys storing instantaneous neuron states
	V.assign(N, V_lk); // All zeros
	I_input.assign(N, 0.0);
	I_leak.assign(N, 0.0);
	I_AMPA.assign(N, 0.0);
	I_GABA.assign(N, 0.0);
	I_NMDA.assign(N, 0.0);
	I_GJ.assign(N, 0.0);
	I_ext.assign(N, 0.0);
	I_K.assign(N, 0.0);
	I_ext_mean.assign(N, 0.0);
	ref_step_left.assign(N, 0);
	ref_steps = (int)round(tau_ref / dt);
	// heterogenous spike-freq-adap
	dg_K_heter.assign(N,0.0);

	spike_hist_tot.reserve(step_tot * 50); // reserve!
	num_ref_pop.reserve(step_tot); // reserve and push_back so that it won't be affected by adapting step_tot
	num_spikes_pop.reserve(step_tot); // reserve and push_back so that it won't be affected by adapting step_tot

	// Random seed (random engine should be feed with DIFFERENT seed at every implementation)
	random_device rd; // random number from operating system for seed
	my_seed = rd(); // record seed
	gen.seed(my_seed);
	//my_seed = 321;
	//cout << "My_seed is: " << my_seed << endl;

	// Runaway killer is initially a sleeper agent
	killer.license = false;
	killer.runaway_killed = false;
	killer.step_killed = -1;

	//
	stats.record = false;
	stats.record_cov = false;
	LFP.record = false;
	spike_freq_adpt = false;	

	// perturbation
	step_perturb = -1;
	spike_removed = -1;

}

void NeuroPop::set_seed(int seed_input) {
	my_seed = seed_input; // This will overwrite the auto generated seed.
	gen.seed(my_seed);
}

const vector< int > & NeuroPop::get_spikes_noise() 
{
	return jh_learn_pop.noise_spikes;
}

const double & NeuroPop::get_noise() 
{
	return jh_learn_pop.noise;
}

const vector< int > & NeuroPop::get_spikes_current()
{
	return spikes_current;
}

const vector< double > & NeuroPop::get_V()
{
	return V;
}


const vector< int > & NeuroPop::get_ref_step_left() {
	return ref_step_left;
}

const bool & NeuroPop::get_runaway_killed()
{
	return killer.runaway_killed;
}


const double & NeuroPop::get_Cm()
{
	return Cm;
}

void NeuroPop::set_neuron_model(int n_mod) {
	neuron_model = n_mod;
}

void NeuroPop::set_ELIF_Params(double elif_delT, double elif_VT) {
	elif.delT = elif_delT;
	elif.V_T = elif_VT;
}

void NeuroPop::init_poisson_pop(double rate) {
	poisson_pop.on=1;
	// note that the input to the simulator is assumed to be a rate in units of spikes per second (Hz)
	// thus we convert to spikes per timestep (the simulator timestep is in units of ms)
	poisson_pop.rate=rate*dt/1000.0; //spikes per timestep

	// In this type of population, each neuron should generate Poisson spike trains with rate lambda
	// To do this we use the fact that in a Poisson process, the time intervals between pairs of subsequent events 
	// are distributed exponentially according to lambda*exp(-lambda*t)
	// So, each time a spike is generated, we then generate a time interval until the next spike and wait until that time is reached before producing a spike 
	// Note: Poisson events, by definition, are generated independently of previous events

	exponential_distribution<double> exp_dist(poisson_pop.rate);	// random numbers distributed according to lambda*exp(-lambda*X)	
	auto interval = bind(exp_dist, ref(gen));
	poisson_pop.next_spike_time.resize(N,0);
	for (int i = 0; i < N; ++i) {
		// generate a time interval until the first spike
		poisson_pop.next_spike_time[i]+=floor(interval());
	}
}

void NeuroPop::recv_I(vector<double>& I, const int pop_ind_pre, const int syn_type)
{
	if (pop_ind_pre < 0) { // if noisy external currents, always send to I_ext regardless of the synapse type
		for (int j = 0; j < N; ++j) {
			I_ext[j] += I[j];
		}
	}
	// AMPA
	else if (syn_type == 0) {
		for (int j = 0; j < N; ++j) {
			I_AMPA[j] += I[j];
		}
	}
	//GABA
	else if (syn_type == 1) {
		for (int j = 0; j < N; ++j) {
			I_GABA[j] += I[j];
		}
	}
	//NMDA
	else if (syn_type == 2) {
		for (int j = 0; j < N; ++j) {
			I_NMDA[j] += I[j];
		}
	}
}

void NeuroPop::start_stats_record()
{
	stats.record = true;

	stats.V_mean.reserve(step_tot);
	stats.V_std.reserve(step_tot);

	stats.I_input_mean.reserve(step_tot);
	stats.I_input_std.reserve(step_tot);

	stats.I_AMPA_time_avg.assign(N, 0.0);
	stats.I_NMDA_time_avg.assign(N, 0.0);
	stats.I_GABA_time_avg.assign(N, 0.0);
	stats.I_ext_time_avg.assign(N, 0.0);
	stats.I_tot_time_mean.assign(N, 0.0);
	stats.I_tot_time_var.assign(N, 0.0);
	stats.V_time_mean.assign(N, 0.0);
	stats.V_time_var.assign(N, 0.0);

	stats.IE_ratio.assign(N, 0.0);
}


void NeuroPop::start_cov_record(const int time_start, const int time_end)
{
	stats.record_cov = true;
	stats.V_time_mean_dumb.assign(N, 0.0);
	stats.V_time_cov.resize(N);
	for (int i = 0; i < N; ++i) {
		stats.V_time_cov[i].assign(N, 0.0);
	}
	// over-write the following two variables
	stats.time_start_cov = time_start;
	stats.time_end_cov = time_end;
}

void NeuroPop::start_LFP_record(const vector <vector<double> >& LFP_neurons_input) {
	if (int(LFP_neurons_input[0].size()) != N) {
		cout << "start_LFP_record failed: LFP_neurons should be 1-by-N logical vector!" << endl;
	}
	else {
		LFP.record = true;
		int n_LFP = int(LFP_neurons_input.size());
		LFP.neurons.resize(n_LFP);
		LFP.data.resize(n_LFP);
		for (int ind = 0; ind < n_LFP; ++ind) {
			LFP.neurons[ind] = LFP_neurons_input[ind];
			LFP.data[ind].reserve(step_tot);
		}
	}
}

void NeuroPop::set_para(string para_str) {
	const char delim = ',';
	if (!para_str.empty()) {
		istringstream para(para_str);
		string para_name, para_value_str;
		double para_value;
		while (getline(para, para_name, delim)) {
			getline(para, para_value_str, delim); // get parameter value (assume double)
			stringstream(para_value_str) >> para_value; // from string to numerical value
			if (para_name.find("Cm") != string::npos) {Cm = para_value;}
			else if (para_name.find("tau_ref") != string::npos) {tau_ref = para_value;}
			else if (para_name.find("V_rt") != string::npos) {V_rt = para_value;}
			else if (para_name.find("V_lk") != string::npos) {V_lk = para_value;}
			else if (para_name.find("V_th") != string::npos) {V_th = para_value;}
			else if (para_name.find("g_lk") != string::npos) {g_lk = para_value;}
			else if (para_name.find("V_ext") != string::npos) {V_ext = para_value;}
			else {cout << "Unrecognized parameter: " << para_name << endl;}
		}
	}
	// re-initialise population!
	init();
}

string NeuroPop::dump_para() {
	stringstream dump;
	dump << "Cm" << delim << Cm << delim << endl;
	dump << "tau_ref" << delim << tau_ref << delim << endl;
	dump << "V_rt" << delim << V_rt << delim << endl;
	dump << "V_lk" << delim << V_lk << delim << endl;
	dump << "V_th" << delim << V_th << delim << endl;
	dump << "g_lk" << delim << g_lk << delim << endl;
	dump << "V_K" << delim << V_K << delim << endl;
	dump << "dg_K" << delim << dg_K << delim << endl;
	dump << "tau_K" << delim << tau_K << delim << endl;
	dump << "V_ext" << delim << V_ext << delim << endl;
	dump << "seed" << delim << my_seed << delim << endl;
	return dump.str();
}


void NeuroPop::random_V(const double p) {
	// Generate random initial condition for V.
	// Generate uniform random distribution
	cout << "Function NeuroPop::random_V(double p) is deprecated!" << endl;
	if (p < 1.0) {
		uniform_real_distribution<double> uniform_dis(0.0, 1.0);
		auto ZeroOne = bind(uniform_dis, ref(gen));

		for (int i = 0; i < N; ++i) {
			// Generate random number.
			//ZeroOne = uniform_dis(gen); // range is [0,1]
			V[i] = V_rt + (V_th - V_rt) / (1.0 - p) * ZeroOne();
		}
	}
	else {cout << "Initial firing rate cannot be 100%!" << endl;}
}

void NeuroPop::set_init_condition(const double r_V0, const double p_fire) {
	// Set V to be uniformly distributed between [V_rt, V_rt + (V_th - V_rt)*r_V0]
	// And then randomly set some of them above firing threshold according to p_fire
	uniform_real_distribution<double> uniform_dis(0.0, 1.0);
	auto ZeroOne = bind(uniform_dis, ref(gen));
	for (int i = 0; i < N; ++i) {
		if (ZeroOne() < p_fire) {
			V[i] = V_th + 1.0; // above threshold for firing
		}
		else {V[i] = V_rt + (V_th - V_rt) * r_V0 * ZeroOne();}
	}
}

void NeuroPop::set_init_V_external(const vector<double>& external_init_V) {
	if (*max_element(V.begin(),V.end()) != V_lk ) {
		cout << "\t\t Warning: random initial V are overwritted by external intial V!" << endl;
	}
	for (int i = 0; i < N; ++i) {
		V[i] = external_init_V[i];
	}
}


void NeuroPop::update_spikes(const int step_current) {
	// Reset currents to be zeros, because they need to be re-calculated at every step
	// no need to reset I_leak
	fill(I_AMPA.begin(), I_AMPA.end(), 0);
	fill(I_GABA.begin(), I_GABA.end(), 0);
	fill(I_NMDA.begin(), I_NMDA.end(), 0);
	fill(I_GJ.begin(), I_GJ.end(), 0);
	fill(I_ext.begin(), I_ext.end(), 0);// fast way to do it, dump in 16 bytes at a time until it gets close to the end

	// Find the firing neurons, record them, reset their potential and set them to be refractory
	spikes_current.clear(); // empty vector
	int spike_counter = 0;
	if(spike_file.on==1) {
		spikes_current.clear(); // empty vector
		int neu=-1;
		for (unsigned int i = 0; i < spike_file.spikes[spike_file.spike_ind].size(); ++i){
			neu= spike_file.spikes[spike_file.spike_ind][i];
			if (ref_step_left[neu] == 0){
				spikes_current.push_back(neu); // record firing neurons
				ref_step_left[neu] = ref_steps; // steps left for being refractory
				spike_counter += 1;
			}
		}
		spike_file.spike_ind++;
		if(spike_file.spike_ind>spike_file.spikes.size()){
			spike_file.spike_ind=0; //go back to start of list
		}
	}
	else if(poisson_pop.on==1){
		// In this type of population, each neuron should generate Poisson spike trains with rate lambda
		// To do this we use the fact that in a Poisson process, the time intervals between pairs of subsequent events 
		// are distributed exponentially according to lambda*exp(-lambda*t)
		// So, each time a spike is generated, we then generate a time interval until the next spike and wait until that time is reached before producing a spike 
		// Note: Poisson events, by definition, are generated independently of previous events

		exponential_distribution<double> exp_dist(poisson_pop.rate);	// random numbers distributed according to lambda*exp(-lambda*X)	
		auto interval = bind(exp_dist, ref(gen));
		for (int i = 0; i < N; ++i) {
			if(poisson_pop.next_spike_time[i]<=step_current){
				spikes_current.push_back(i); // record firing neurons
				spike_counter += 1;
				//now generate a new spike interval
				poisson_pop.next_spike_time[i]+=floor(interval());
			}
		}
	}
	else{
		for (int i = 0; i < N; ++i) {
			if (ref_step_left[i] == 0 && V[i] >= V_th) {
				spikes_current.push_back(i); // record firing neurons
				V[i] = V_rt; // reset potential
				ref_step_left[i] = ref_steps; // steps left for being refractory
				spike_counter += 1;
			}
		}
	}

	if(jh_learn_pop.on){
		jh_learn_pop.noise_spikes.clear();
		jh_learn_pop.noise_spikes.resize(spikes_current.size(),0); //initialise to no noise
		if(jh_learn_pop.noise>0.0){
			//noise spikes
			uniform_real_distribution<double> uniform_dis(0.0, 1.0);
			auto ZeroOne = bind(uniform_dis,ref(gen));
			for(unsigned int i=0;i<spikes_current.size();i++){
				if (ZeroOne() < jh_learn_pop.noise){
					jh_learn_pop.noise_spikes[i]=1; //noise this spike
				}
			}
		}	
	}

	// perturbation: remove the last spike
	if (step_current == step_perturb) {
		if (spikes_current.size() > 0) {
			spike_removed = spikes_current.back(); // record the spike that is removed
			spikes_current.pop_back(); // remove the last element as the perturbation
			spike_counter -= 1;
			// output
			cout << endl << "Perturbation: the spike of neuron #" << spike_removed << " removed at time step " << step_perturb << endl;
		}
		else { step_perturb += 1; } // if there is no spike to be removed at this step, try the next step
	}

	// record number of spikes
	num_spikes_pop.push_back(spike_counter);

	// Collect total spike history data
	if (spike_counter > 0) {
		// .end(): ONE ELEMENT PAST to the last location
		// see http://www.cs.northwestern.edu/~riesbeck/programming/c++/stl-iterators.html
		copy(spikes_current.begin(), spikes_current.end(), back_inserter(spike_hist_tot));
	}

	// Refraction count-down and record number of refractory neurons
	int num_ref_temp = 0;
	for (int i = 0; i < N; ++i) {
		if (ref_step_left[i] > 0) {
			ref_step_left[i] -= 1;
			num_ref_temp += 1;
		}
	}
	num_ref_pop.push_back(num_ref_temp);


	// runaway check
	runaway_check(step_current);

}

void NeuroPop::generate_I_ext() {

	// Gaussian white external currents
	if (I_ext_mean.size() != 0) {
		if (I_ext_std.size() != 0) {
			double one_on_sqrt_dt = 1.0 / sqrt(dt); // Here sqrt_dt is on the denominator because I_ext will be multiplied by dt later.
			// Gaussian random generator
			normal_distribution<double> nrm_dist(0.0, 1.0);
			auto gaus = bind(nrm_dist, ref(gen));

			for (int i = 0; i < N; ++i) {
				I_ext[i] += I_ext_mean[i] + gaus() * I_ext_std[i] * one_on_sqrt_dt; // be careful about the sqrt(dt) term (Wiener Process)
			}
		}
		else {
			for (int i = 0; i < N; ++i) {
				I_ext[i] += I_ext_mean[i];
			}
		}
	}

	// Gaussian white external conductance
	if (g_ext_mean.size() != 0) {
		if (g_ext_std.size() != 0) {
			double one_on_sqrt_dt = 1.0 / sqrt(dt); // Here sqrt_dt is on the denominator because I_ext will be multiplied by dt later.
			// Gaussian random generator
			normal_distribution<double> nrm_dist(0.0, 1.0);
			auto gaus = bind(nrm_dist, ref(gen));
			for (int i = 0; i < N; ++i) {
				I_ext[i] += -(g_ext_mean[i] + gaus() * g_ext_std[i] * one_on_sqrt_dt) * (V[i] - V_ext); // be careful about the sqrt(dt) term (Wiener Process)
			}
		}
		else {
			for (int i = 0; i < N; ++i) {
				I_ext[i] += -g_ext_mean[i] * (V[i] - V_ext);
			}
		}
	}
}

void NeuroPop::get_current_from_file() {

	for (int i = 0; i < N; i++) {
		I_ext[i] += current_file.mean_curr * current_file.current[current_file.current_ind][i];
	}
	current_file.framestep++;
	if (current_file.framestep >= current_file.steps_per_frame) {
		current_file.framestep = 0;
		current_file.current_ind++;
		if (current_file.current_ind >= current_file.current.size()) {
			current_file.current_ind = 0; //go back to start of list
		}
	}
}

void NeuroPop::update_V(const int step_current) {
	// This function updates menbrane potentials for non-refractory neurons

	// Generate external currents
	generate_I_ext();
	// add external current file in [start_step, end_step)
	if (current_file.on && (step_current >= current_file.start_step) && (step_current < current_file.end_step)) {
		get_current_from_file();
	}


	// potassium conductance for spike-frequency adaptation
	if (spike_freq_adpt == true) {
		if (spike_freq_adpt_heter.on && (step_current >= spike_freq_adpt_heter.start_step) && (step_current < spike_freq_adpt_heter.end_step)){
		// heterogenous potassium conductance for spike-frequency adaptation
			for (unsigned int ind = 0; ind < spikes_current.size(); ++ind) {
				g_K[ spikes_current[ind] ] += dg_K_heter[ spikes_current[ind] ];
			}
		}
		// homogenous potassium conductance for spike-frequency adaptation
		else{
			for (unsigned int ind = 0; ind < spikes_current.size(); ++ind) {
				g_K[ spikes_current[ind] ] += dg_K;
			}
		}		
		
		for (int i = 0; i < N; ++i) {
			g_K[i] *= exp_K_step;
			I_K[i] = -g_K[i] * (V[i] - V_K);
		}
	}

	


	// Collect Currents from all pre-synapses (for MPI job)!!!!!!!!!!!!!!!



	// Data sampling, which must be done here!
	// sample_data(step_current);  // this is deprecated due to poor memory performance
	output_sampled_data_real_time_HDF5(step_current);

	// update menbrane potentials
	double Vdot;
	for (int i = 0; i < N; ++i) {
		I_input[i] = I_AMPA[i] + I_GABA[i] + I_NMDA[i] + I_GJ[i] + I_ext[i] + I_K[i];
		if (ref_step_left[i] == 0) { // Only update the non-refractory neurons
			// leaky current
			if (neuron_model == 0) {
				//LIF
				I_leak[i] = -g_lk * (V[i] - V_lk);
			}
			else if (neuron_model == 1) {
				//exponential LIF
				I_leak[i] = -g_lk * (V[i] - V_lk) + g_lk * elif.delT * exp( ( V[i] - elif.V_T ) / elif.delT );
			}
			// using simple Euler method
			Vdot = (I_leak[i] + I_input[i]) / Cm;
			V[i] += Vdot * dt;
			// Note that delta-function coupling is very different from the above conductance-based model!
		}
	}

	// record mean and std of membrane potentials
	record_stats(step_current);
	record_LFP();

}


void NeuroPop::add_sampling(const vector<int>& sample_neurons_input, const vector<bool>& sample_type_input, const vector<bool>& sample_time_points_input) {
	sample.neurons = sample_neurons_input;
	sample.type = sample_type_input;
	sample.time_points = sample_time_points_input;


	// initialise
	sample.N_steps = 0;
	for (int tt = 0; tt < step_tot; ++tt) {
		if (sample.time_points[tt]) {sample.N_steps += 1;}
	}
	sample.N_neurons = sample.neurons.size();

	int sample_type_tot = sample.type.size(); // 8 different data types


	sample.data.resize(sample_type_tot);
	for (int c = 0; c < sample_type_tot; ++c) {
		if (sample.type[c]) {
			sample.data[c].resize(sample.N_neurons);
			for (int i = 0; i < sample.N_neurons; ++i) {
				sample.data[c][i].reserve(sample.N_neurons); // reserve and push_back so that it won't be affected by adapting step_tot
			}
		}
	}
}



void NeuroPop::sample_data(const int step_current) {

	if (!sample.neurons.empty()) {
		if (sample.time_points[step_current]) { // push_back is amazing
			for (int i = 0; i < sample.N_neurons; ++i) { // performance issue when sampling many neurons?
				int ind_temp = sample.neurons[i];
				if (sample.type[0]) {sample.data[0][i].push_back( V[ind_temp] );}
				if (sample.type[1]) {sample.data[1][i].push_back( I_leak[ind_temp] );}
				if (sample.type[2]) {sample.data[2][i].push_back( I_AMPA[ind_temp] );}
				if (sample.type[3]) {sample.data[3][i].push_back( I_GABA[ind_temp] );}
				if (sample.type[4]) {sample.data[4][i].push_back( I_NMDA[ind_temp] );}
				if (sample.type[5]) {sample.data[5][i].push_back( I_GJ[ind_temp] );}
				if (sample.type[6]) {sample.data[6][i].push_back( I_ext[ind_temp] );}
				if (sample.type[7]) {sample.data[7][i].push_back( I_K[ind_temp] );}
			}
		}
	}

}


void NeuroPop::set_gaussian_I_ext(const vector<double>& mean, const vector<double>& std) {
	I_ext_mean = mean;
	I_ext_std = std;

	double max_std = *max_element(I_ext_std.begin(), I_ext_std.end());
	if (max_std == 0.0) {
		I_ext_std.resize(0);
	}
}

void NeuroPop::set_gaussian_g_ext(const vector<double>& mean, const vector<double>& std) {
	g_ext_mean = mean;
	g_ext_std = std;

	double max_std = *max_element(g_ext_std.begin(), g_ext_std.end());
	if (max_std == 0.0) {
		g_ext_std.resize(0);
	}
}


void NeuroPop::add_perturbation(const int step_perturb_input) {
	step_perturb = step_perturb_input;
}

void NeuroPop::add_spike_freq_adpt() {
	spike_freq_adpt = true;
	exp_K_step = exp( -dt / tau_K );
	g_K.assign(N, 0.0);
}

void NeuroPop::set_spike_freq_adpt_para(const double dg_K_input) {
	dg_K = dg_K_input; // default value 0.01 (uS=miuSiemens)
}

void NeuroPop::set_spike_freq_adpt_para_heter(const vector<double>& dg_K_heter_input, const int start_step, const int end_step) {
	dg_K_heter = dg_K_heter_input; // default value 0.01 (uS=miuSiemens)
	spike_freq_adpt_heter.on = true;
	spike_freq_adpt_heter.start_step = start_step;
	spike_freq_adpt_heter.end_step = end_step;
}


void NeuroPop::init_runaway_killer(const double min_ms, const double Hz, const double Hz_ms)
{
	killer.min_pop_size = 100;
	if (N > killer.min_pop_size) {
		killer.license = true;
		killer.min_steps = int(round(min_ms / dt));
		killer.runaway_Hz = Hz;
		killer.Hz_steps = int(round(Hz_ms / dt));
	}
}

void NeuroPop::runaway_check(const int step_current)
{
	if (killer.license == true && killer.runaway_killed == false && step_current > killer.min_steps && step_current > killer.Hz_steps) {
		// find mean value of num_ref over the last runaway_steps
		vector<int>::const_iterator first, last;
		// first element to be accumulated
		first = num_spikes_pop.begin() + (step_current - killer.Hz_steps + 1);
		// one element pass the last element to be accumulated
		last = num_spikes_pop.begin() + (step_current + 1);
		double mean_Hz = accumulate(first, last, 0.0) / (killer.Hz_steps * dt * 0.001 * N); // 0.001 for converting from ms to sec.
		//be careful!! accumulate range is : [first,last)
		if (mean_Hz >= killer.runaway_Hz) {
			killer.runaway_killed = true;
			killer.step_killed = step_current;
			cout << "warning: runaway killed at " << step_current*dt << " (ms) in population" << pop_ind << flush;
			cout << "\t with firing rate at " << mean_Hz << " Hz." << flush;
		}
	}
}


void NeuroPop::add_JH_Learn(double noise){
	jh_learn_pop.on=true;
	jh_learn_pop.rhatE.resize(N,0);
	jh_learn_pop.rhatI.resize(N,0);


	jh_learn_pop.noise=noise;
}

void NeuroPop::reset_rhat(){
	if(jh_learn_pop.on){
		fill(jh_learn_pop.rhatE.begin(), jh_learn_pop.rhatE.end(), 0.0);
		fill(jh_learn_pop.rhatI.begin(), jh_learn_pop.rhatI.end(), 0.0);
	}
}

void NeuroPop::get_all_rhat_JHLearn(vector<ChemSyn*> &ChemSynArray,const int step_current){
	if(jh_learn_pop.on){
		if (sample.type[8]){
			if (sample.time_points[step_current]){ 
				if (!sample.neurons.empty() && sample.file_type == 2){			
					reset_rhat();
					for (unsigned int syn_ind = 0; syn_ind < ChemSynArray.size(); ++syn_ind){
						if(ChemSynArray[syn_ind]->get_pop_ind_pre()==pop_ind){
							vector <double> temp;
							temp=ChemSynArray[syn_ind]->get_all_rhat_JH_Learn(sample.neurons);
							if(ChemSynArray[syn_ind]->get_post_neu_type()==0){ //AMPA excitatory
								for(unsigned int i=0;i<sample.neurons.size();i++){
									jh_learn_pop.rhatE[sample.neurons[i]]+=temp[sample.neurons[i]];
								}
							}
							else if(ChemSynArray[syn_ind]->get_post_neu_type()==1){ //GABA inhibitory
								for(unsigned int i=0;i<sample.neurons.size();i++){
									jh_learn_pop.rhatI[sample.neurons[i]]+=temp[sample.neurons[i]];
								}
							}
						}
					}
				}
			}
		}
	}
}


void NeuroPop::record_LFP() {
	if (LFP.record) {
		for (unsigned int ind = 0; ind < LFP.neurons.size(); ++ind) {
			LFP.data[ind].push_back(0.0);
			for (int i = 0; i < N; ++i) {
				if (LFP.neurons[ind][i] > 0) {
					LFP.data[ind].back() += LFP.neurons[ind][i] * (abs(I_AMPA[i]) + abs(I_GABA[i]));
				}
			}
		}
	}
}


void NeuroPop::record_stats(const int step_current) {
	if (stats.record) {
		// V and I_input averaged over the population
		double mean_tmp_V, var_tmp_V, mean_tmp_I, var_tmp_I;
		Welford_online(V, mean_tmp_V, var_tmp_V);
		Welford_online(I_input, mean_tmp_I, var_tmp_I);
		stats.V_mean.push_back(mean_tmp_V);
		stats.V_std.push_back(sqrt(var_tmp_V));
		stats.I_input_mean.push_back(mean_tmp_I);
		stats.I_input_std.push_back(sqrt(var_tmp_I));

		// online mean and var calculation: Welford's method (1962, Technometrixcs)
		bool is_end = step_current == (step_tot - 1);
		int k = step_current;
		Welford_online(I_input, stats.I_tot_time_mean, stats.I_tot_time_var, k, is_end);
		Welford_online(I_GABA, stats.I_GABA_time_avg, k);
		Welford_online(I_NMDA, stats.I_NMDA_time_avg, k);
		Welford_online(I_AMPA, stats.I_AMPA_time_avg, k);
		Welford_online(I_ext, stats.I_ext_time_avg, k);
		Welford_online(V, stats.V_time_mean, stats.V_time_var,  k, is_end);

		// get time average for each neuron
		if (step_current == step_tot - 1) { // at the end of the last time step
			for (int i = 0; i < N; ++i) {
				// be careful here, for IE_ratio, I_ext is assumed to be always excitatory and I_GJ is not considered
				// also, the only source of I_ext is generated internally
				stats.IE_ratio[i] = stats.I_GABA_time_avg[i] / (stats.I_AMPA_time_avg[i] + stats.I_NMDA_time_avg[i] + stats.I_ext_time_avg[i]);
			}
		}
	}

	if (stats.record_cov) {
		if (step_current >= stats.time_start_cov && step_current <= stats.time_end_cov) {
			bool is_end = step_current == (stats.time_end_cov - 1);
			int k = step_current - stats.time_start_cov;
			Welford_online(V, stats.V_time_mean_dumb, stats.V_time_cov,  k, is_end);
		}
	}
}

void Welford_online(const vector<double>& data, double& M, double& S) {
	M = 0.0;
	S = 0.0;
	double M_old, x;
	for (unsigned int i = 0; i < data.size(); ++i) {
		M_old = M;
		x = data[i];
		M += (x - M_old) / double(i + 1.0);
		S += (x - M_old) * (x - M);
	}
	S = S / (double(data.size()) - 1.0); // N-1, sample variance
}

void Welford_online(const vector<double>& new_data, vector<double>& M, vector<double>& S, const int K, const bool is_end) {
	// Note that K follows C++ index, K = 0, 1, ....
	// online mean and var calculation: Welford's method (1962, Technometrixcs)
	double M_old, x;
	for (unsigned int i = 0; i < M.size(); ++i) {
		M_old = M[i];
		x = new_data[i];
		M[i] += (x - M_old) / double(K + 1.0);
		S[i] += (x - M_old) * (x - M[i]);
		if (is_end == true) {
			S[i] = S[i] / double(K); // N-1, sample variance
		}
	}

}

void Welford_online(const vector<double>& new_data, vector<double>& M, vector< vector <double> >& Cov, const int K, const bool is_end) {
	double M_old, x;
	int N_tmp = int(M.size());
	vector<double> M_add;
	M_add.resize(N_tmp);

	for (int i = 0; i < N_tmp; ++i) {
		M_old = M[i];
		x = new_data[i];
		M_add[i] = (x - M_old) / double(K + 1.0);
		M[i] += M_add[i];
	}
	// covariance
	/*
	for (int i = 0; i < N_tmp; ++i){
		for (int j = i; j < N_tmp; ++j){
			Cov[i][j] -= Cov[i][j] / double(K + 1.0);
			Cov[i][j] += double(K) * M_add[i] * M_add[j];
			if (is_end == true){
				Cov[i][j] = Cov[i][j] * double(K + 1.0) / double(K);
				Cov[j][i] = Cov[i][j];
			}
		}
	}
	*/
	vector<double> Cov_add;
	Cov_add.resize(N_tmp);
	for (int i = 0; i < N_tmp; ++i) {
		// Cov[i][j] -= Cov[i][j] / double(K + 1.0);
		transform(Cov[i].begin(), Cov[i].end(), Cov[i].begin(), [K](double x) {return x - (x / double(K + 1.0));});
		// Cov_add = double(K) * M_add[i];
		fill(Cov_add.begin(), Cov_add.end(), double(K) * M_add[i]);
		// Cov_add = Cov_add * M_add[j];
		transform(Cov_add.begin(), Cov_add.end(), M_add.begin(), Cov_add.begin(), multiplies<double>());
		// Cov[i][j] += Cov_add;
		transform(Cov[i].begin(), Cov[i].end(), Cov_add.begin(), Cov[i].begin(), plus<double>());
		if (is_end == true) {
			// Cov[i][j] = Cov[i][j] * double(K + 1.0) / double(K);
			transform(Cov[i].begin(), Cov[i].end(), Cov[i].begin(), [K](double x) {return x * double(K + 1.0) / double(K);});
		}
	}

}


void Welford_online(const vector<double>& new_data, vector<double>& M, const int K) {
	// Note that K follows C++ index, K = 0, 1, ....
	// online mean and var calculation: Welford's method (1962, Technometrixcs)
	double M_old, x;
	for (unsigned int i = 0; i < M.size(); ++i) {
		M_old = M[i];
		x = new_data[i];
		M[i] += (x - M_old) / double(K + 1.0);
	}
}


void NeuroPop::load_file_current_input(string fname) {
	cout << "\t\tLoading currents from file... " << fname;
	H5File file(  fname, H5F_ACC_RDONLY );
	vector< int> neurons, t, I;
	double fps;
	read_matrix_HDF5(file, string("neurons"), current_file.neurons);
	read_matrix_HDF5(file, string("current"), current_file.current);
	fps = read_scalar_HDF5<double>(file, string("frame_rate"));
	current_file.mean_curr = read_scalar_HDF5<double>(file, string("mean_curr"));
	current_file.steps_per_frame = 1000 / (double(fps) * dt);
	current_file.file_name = fname;
	current_file.current_ind = 0;
	current_file.on = 1;
	current_file.start_step = read_scalar_HDF5<int>(file, string("start_step"));
	current_file.end_step = read_scalar_HDF5<int>(file, string("end_step"));
	cout << "done.\n";
}

void NeuroPop::load_file_spike_input(string fname) {
	cout<<"\t\tLoading spikes from file..."<<fname;
	H5File file(  fname, H5F_ACC_RDONLY );
	vector< int> x,y,t,pol;
	int max_x=1;
	int max_y=1;
	read_vector_HDF5(file, string("x"), x);
	read_vector_HDF5(file, string("y"), y);
	read_vector_HDF5(file, string("t"), t);
	read_vector_HDF5(file, string("pol"), pol);
	max_x=read_scalar_HDF5<int>(file,string("max_x"));
	max_y=read_scalar_HDF5<int>(file,string("max_y"));

	vector<int> sp_neu;
	sp_neu.resize(max_x*max_y,0); //  used to mark if already a spike for that neuron to limit to 1 spike per timestep
	// Now convert list of events into sets of spikes for each timestep dt
	unsigned long int  i=0;
	int tmax;
	tmax=t[1]+int(1000*dt); // convert from units of microseconds to milliseconds
	vector< int> tmp_spikes;
	while(i<t.size()){
		while((t[i]<tmax)&(i<t.size())){
			if(sp_neu[x[i]+max_x*y[i]]==0){ //check if this neuron has already spiked this timestep
				tmp_spikes.push_back(x[i]+max_x*y[i]);
				sp_neu[x[i]+max_x*y[i]]=1; //mark this neuron as having spiked
			}
			i++;
		}
		spike_file.spikes.push_back(tmp_spikes);
		tmp_spikes.clear();
		tmax+=int(1000*dt);

		sp_neu.assign(sp_neu.size(), 0); // reset to 0 so all neurons are marked as NOT having spiked for next timestep
	}
	spike_file.file_name=fname;
	spike_file.spike_ind=0;
	spike_file.on=1;
	cout<<"done.\n";
}


void NeuroPop::import_restart(H5File& file, int pop_ind, string out_filename) {
	string str;
	string pop_n = "/pops/pop" + to_string(pop_ind) + "/";

	neuron_model = read_scalar_HDF5<int>(file, pop_n + string("neuron_model"));
	if (neuron_model == 1) {
		elif.delT = read_scalar_HDF5<double>(file, pop_n + string("/ELIF/delT"));
		elif.V_T = read_scalar_HDF5<double>(file, pop_n + string("/ELIF/V_T"));
	}

	pop_ind = read_scalar_HDF5<double>(file, pop_n + string("pop_ind"));
	N = read_scalar_HDF5<double>(file, pop_n + string("N"));
	dt = read_scalar_HDF5<double>(file, pop_n + string("dt"));
	step_tot = read_scalar_HDF5<double>(file, pop_n + string("step_tot"));
	tau_ref = read_scalar_HDF5<double>(file, pop_n + string("tau_ref"));
	Cm = read_scalar_HDF5<double>(file, pop_n + string("Cm"));
	V_rt = read_scalar_HDF5<double>(file, pop_n + string("V_rt"));
	V_lk = read_scalar_HDF5<double>(file, pop_n + string("V_lk"));
	V_th = read_scalar_HDF5<double>(file, pop_n + string("V_th"));
	g_lk = read_scalar_HDF5<double>(file, pop_n + string("g_lk"));

	read_vector_HDF5(file, pop_n + "V", V);
	read_vector_HDF5(file, pop_n + "I_leak", I_leak);
	read_vector_HDF5(file, pop_n + "I_input", I_input);
	read_vector_HDF5(file, pop_n + "I_AMPA", I_AMPA);
	read_vector_HDF5(file, pop_n + "I_GABA", I_GABA);
	read_vector_HDF5(file, pop_n + "I_NMDA", I_NMDA);
	read_vector_HDF5(file, pop_n + "I_GJ", I_GJ);
	read_vector_HDF5(file, pop_n + "I_ext", I_ext);

	read_vector_HDF5(file, pop_n + string("I_ext_mean"), I_ext_mean);

	if (dataset_exist_HDF5(file, pop_n + string("I_ext_std"))) {

		read_vector_HDF5(file, pop_n + string("I_ext_std"), I_ext_std);
	}

	read_vector_HDF5(file, pop_n + string("g_ext_mean"), g_ext_mean);

	if (dataset_exist_HDF5(file, pop_n + string("g_ext_std"))) {

		read_vector_HDF5(file, pop_n + "g_ext_std", g_ext_std);
	}

	ref_steps = read_scalar_HDF5<int>(file, pop_n + string("ref_steps"));

	read_vector_HDF5(file, pop_n + "ref_step_left", ref_step_left);
	if (dataset_exist_HDF5(file, pop_n + string("spikes_current"))) {
		read_vector_HDF5(file, pop_n + "spikes_current", spikes_current);
	}
	// read_vector_HDF5(file,pop_n+"spike_hist_tot",spike_hist_tot);
	// read_vector_HDF5(file,pop_n+"num_spikes_pop",num_spikes_pop);
	// read_vector_HDF5(file,pop_n+"num_ref_pop",num_ref_pop);

	str = pop_n + "/Stats/";
	if (group_exist_HDF5(file, str)) {
		stats.record = read_scalar_HDF5<bool>(file, str + "record");
		start_stats_record();
		// read_vector_HDF5(file,str+"V_mean",stats.V_mean);
		// read_vector_HDF5(file,str+"V_std",stats.V_std);
		// read_vector_HDF5(file,str+"I_input_mean",stats.I_input_mean);
		// read_vector_HDF5(file,str+"I_input_std",stats.I_input_std);
		// read_vector_HDF5(file,str+"I_AMPA_time_avg",stats.I_AMPA_time_avg);
		// read_vector_HDF5(file,str+"I_NMDA_time_avg",stats.I_NMDA_time_avg);
		// read_vector_HDF5(file,str+"I_GABA_time_avg",stats.I_GABA_time_avg);
		// read_vector_HDF5(file,str+"I_ext_time_avg",stats.I_ext_time_avg);
		// read_vector_HDF5(file,str+"I_tot_time_var",stats.I_tot_time_var);
		// read_vector_HDF5(file,str+"I_tot_time_mean",stats.I_tot_time_mean);
		// read_vector_HDF5(file,str+"IE_ratio",stats.IE_ratio);
	}

	str = pop_n + "/LFP/";
	if (group_exist_HDF5(file, str)) {
		LFP.record = read_scalar_HDF5<bool>(file, str + "record");
		read_matrix_HDF5(file, str + "neurons", LFP.neurons);
		start_LFP_record(LFP.neurons);
		// read_matrix_HDF5(file,str+"data",LFP.data);
	}



	V_ext = read_scalar_HDF5<double>(file, pop_n + string("V_ext"));
	spike_freq_adpt = read_scalar_HDF5<bool>(file, pop_n + string("spike_freq_adpt"));	
	if (dataset_exist_HDF5(file, pop_n + string("g_K"))) {
		read_vector_HDF5(file, pop_n + string("g_K"), g_K);
	}
	read_vector_HDF5(file, pop_n + "I_K", I_K);

	V_K = read_scalar_HDF5<double>(file, pop_n + string("V_K"));
	dg_K = read_scalar_HDF5<double>(file, pop_n + string("dg_K"));
	// heterogenous spike-freq-adap
	read_vector_HDF5(file, pop_n + string("dg_K_heter"),dg_K_heter);
	tau_K = read_scalar_HDF5<double>(file, pop_n + string("tau_K"));
	exp_K_step = read_scalar_HDF5<double>(file, pop_n + string("exp_K_step"));

	str = pop_n + "/sample/";
	if (group_exist_HDF5(file, str)) {
		sample.file_type = read_scalar_HDF5<int>(file, str + string("file_type"));
		// read_string_HDF5(file, str+string("file_name"), sample.file_name);
		// sample.ctr=read_scalar_HDF5<int>(file,str+string("ctr"));
		read_vector_HDF5(file, str + "neurons", sample.neurons);
		read_vector_HDF5(file, str + "type", sample.type);
		read_vector_HDF5(file, str + "time_points", sample.time_points);
		sample.N_steps = read_scalar_HDF5<int>(file, str + string("N_steps"));
		sample.N_neurons = read_scalar_HDF5<int>(file, str + string("N_neurons"));
		add_sampling_real_time_HDF5(sample.neurons, sample.type, sample.time_points, out_filename);
		// read_3Dmatrix_HDF5(file,str+"data",sample.data);

	}

	step_perturb = read_scalar_HDF5<int>(file, pop_n + string("step_perturb"));
	spike_removed = read_scalar_HDF5<int>(file, pop_n + string("spike_removed"));
	my_seed = read_scalar_HDF5<int>(file, pop_n + string("my_seed"));

	str = pop_n + "/killer/";
	if (group_exist_HDF5(file, str)) {
		killer.license = read_scalar_HDF5<bool>(file, str + string("license"));
		killer.runaway_killed = read_scalar_HDF5<bool>(file, str + string("runaway_killed"));
		killer.step_killed = read_scalar_HDF5<int>(file, str + string("step_killed"));
		killer.Hz_steps = read_scalar_HDF5<int>(file, str + string("Hz_steps"));
		killer.runaway_Hz = read_scalar_HDF5<double>(file, str + string("runaway_Hz"));
		killer.min_steps = read_scalar_HDF5<int>(file, str + string("min_steps"));
		killer.min_pop_size = read_scalar_HDF5<int>(file, str + string("min_pop_size"));
	}

	str = pop_n + "/jh_learn_pop/";
	if (group_exist_HDF5(file, str)) {
		jh_learn_pop.on=read_scalar_HDF5<bool>(file,str+string("on"));
		read_vector_HDF5(file,str+string("rhatE"),jh_learn_pop.rhatE);	
		read_vector_HDF5(file,str+string("rhatI"),jh_learn_pop.rhatI);
		read_vector_HDF5(file,str+string("noise_spikes"),jh_learn_pop.noise_spikes);
		jh_learn_pop.noise=read_scalar_HDF5<double>(file,str+string("noise"));	
	}
	str = pop_n + "/spike_file/";
	if (group_exist_HDF5(file, str)) {
		spike_file.on = 1;
		read_string_HDF5(file, str + string("file_name"), spike_file.file_name);
		load_file_spike_input(spike_file.file_name);
		spike_file.spike_ind = read_scalar_HDF5<unsigned int>(file, str + string("spike_ind")); //must go after load_file_spike_input
	}

	str = pop_n + "/poisson_pop/";
	if (group_exist_HDF5(file, str)) {
		poisson_pop.on = 1;
		poisson_pop.rate = read_scalar_HDF5<double>(file, str + string("rate")); 
		read_vector_HDF5(file,str+string("next_spike_time"),poisson_pop.next_spike_time);
	}

	str = pop_n + "/current_file/";
	if (group_exist_HDF5(file, str)) {
		current_file.on = 1;
		read_string_HDF5(file, str + string("file_name"), current_file.file_name);
		load_file_current_input(current_file.file_name);
		current_file.current_ind = read_scalar_HDF5<unsigned int>(file, str + string("current_ind")); //must go after load_file_current_input
		current_file.framestep = read_scalar_HDF5<unsigned int>(file, str + string("framestep")); //must go after load_file_current_input
	}

}

void NeuroPop::export_restart(Group& group) {

	string pop_n = "/pops/pop" + to_string(pop_ind) + "/";
	Group group_pop = group.createGroup(pop_n);

	write_scalar_HDF5(group_pop, neuron_model, string("neuron_model"));
	if (neuron_model == 1) {
		string str = pop_n + "/ELIF/";
		Group group_ELIF = group_pop.createGroup(str);
		write_scalar_HDF5(group_ELIF, elif.delT, string("delT"));
		write_scalar_HDF5(group_ELIF, elif.delT, string("V_T"));
	}

	write_scalar_HDF5(group_pop, pop_ind, string("pop_ind"));
	write_scalar_HDF5(group_pop, N, string("N"));
	write_scalar_HDF5(group_pop, dt, string("dt"));
	write_scalar_HDF5(group_pop, step_tot, string("step_tot"));
	write_scalar_HDF5(group_pop, tau_ref, string("tau_ref"));
	write_scalar_HDF5(group_pop, Cm, string("Cm"));
	write_scalar_HDF5(group_pop, V_rt, string("V_rt"));
	write_scalar_HDF5(group_pop, V_lk, string("V_lk"));
	write_scalar_HDF5(group_pop, V_th, string("V_th"));
	write_scalar_HDF5(group_pop, g_lk, string("g_lk"));

	write_vector_HDF5(group_pop, V, "V");
	write_vector_HDF5(group_pop, I_leak, "I_leak");
	write_vector_HDF5(group_pop, I_input, "I_input");
	write_vector_HDF5(group_pop, I_AMPA, "I_AMPA");
	write_vector_HDF5(group_pop, I_GABA, "I_GABA");
	write_vector_HDF5(group_pop, I_NMDA, "I_NMDA");
	write_vector_HDF5(group_pop, I_GJ, "I_GJ");
	write_vector_HDF5(group_pop, I_ext, "I_ext");
	write_scalar_HDF5(group_pop, ref_steps, string("ref_steps"));
	write_vector_HDF5(group_pop, ref_step_left, "ref_step_left");
	if (spikes_current.size() != 0) {
		write_vector_HDF5(group_pop, spikes_current, "spikes_current");
	}
	write_vector_HDF5(group_pop, spike_hist_tot, "spike_hist_tot");
	write_vector_HDF5(group_pop, num_spikes_pop, "num_spikes_pop");
	write_vector_HDF5(group_pop, num_ref_pop, "num_ref_pop");

	if (stats.record) {
		string str = pop_n + "/Stats/";
		Group group_stats = group_pop.createGroup(str);
		write_scalar_HDF5(group_stats, stats.record, string("record"));
		write_vector_HDF5(group_stats, stats.V_mean, "V_mean");
		write_vector_HDF5(group_stats, stats.V_std, "V_std");
		write_vector_HDF5(group_stats, stats.I_input_mean, "I_input_mean");
		write_vector_HDF5(group_stats, stats.I_input_std, "I_input_std");
		write_vector_HDF5(group_stats, stats.I_AMPA_time_avg, "I_AMPA_time_avg");
		write_vector_HDF5(group_stats, stats.I_NMDA_time_avg, "I_NMDA_time_avg");
		write_vector_HDF5(group_stats, stats.I_GABA_time_avg, "I_GABA_time_avg");
		write_vector_HDF5(group_stats, stats.I_ext_time_avg, "I_ext_time_avg");
		write_vector_HDF5(group_stats, stats.I_tot_time_mean, "I_tot_time_mean");
		write_vector_HDF5(group_stats, stats.I_tot_time_var, "I_tot_time_var");
		write_vector_HDF5(group_stats, stats.IE_ratio, "IE_ratio");
	}

	if (LFP.record) {
		string str = pop_n + "/LFP/";
		Group group_LFP = group_pop.createGroup(str);
		write_scalar_HDF5(group_LFP, LFP.record, string("record"));
		write_matrix_HDF5(group_LFP, LFP.neurons, "neurons");
		write_matrix_HDF5(group_LFP, LFP.data, "data");
	}

	write_vector_HDF5(group_pop, I_ext_mean, string("I_ext_mean"));
	if (I_ext_std.size() != 0) {
		write_vector_HDF5(group_pop, I_ext_std, string("I_ext_std"));
	}
	write_vector_HDF5(group_pop, g_ext_mean, string("g_ext_mean"));
	if (g_ext_std.size() != 0) {
		write_vector_HDF5(group_pop, g_ext_std, "g_ext_std");
	}

	write_scalar_HDF5(group_pop, V_ext, string("V_ext"));

	write_scalar_HDF5(group_pop, spike_freq_adpt, string("spike_freq_adpt"));	
	write_vector_HDF5(group_pop, g_K, string("g_K"));
	write_vector_HDF5(group_pop, I_K, "I_K");

	write_scalar_HDF5(group_pop, V_K, string("V_K"));
	write_scalar_HDF5(group_pop, dg_K, string("dg_K"));
	write_vector_HDF5(group_pop, dg_K_heter, string("dg_K_heter"));
	write_scalar_HDF5(group_pop, tau_K, string("tau_K"));
	write_scalar_HDF5(group_pop, exp_K_step, string("exp_K_step"));

	if (sample.file_type == 2) {
		string str = pop_n + "/sample/";
		Group group_sample = group_pop.createGroup(str);
		write_scalar_HDF5(group_sample, sample.file_type, string("file_type"));
		write_string_HDF5(group_sample, sample.file_name, string("file_name"));
		write_scalar_HDF5(group_sample, sample.ctr, string("ctr"));
		write_vector_HDF5(group_sample, sample.neurons, "neurons");
		write_vector_HDF5(group_sample, sample.type, "type");
		write_vector_HDF5(group_sample, sample.time_points, "time_points");

		write_scalar_HDF5(group_sample, sample.N_steps, string("N_steps"));
		write_scalar_HDF5(group_sample, sample.N_neurons, string("N_neurons"));
		// write_3Dmatrix_HDF5(group_sample,data,"data");
	}

	write_scalar_HDF5(group_pop, step_perturb, string("step_perturb"));
	write_scalar_HDF5(group_pop, spike_removed, string("spike_removed"));
	write_scalar_HDF5(group_pop, my_seed, string("my_seed"));

	if (killer.license) {
		string str = pop_n + "/killer/";
		Group group_killer = group_pop.createGroup(str);
		write_scalar_HDF5(group_killer, killer.license, string("license"));
		write_scalar_HDF5(group_killer, killer.runaway_killed, string("runaway_killed"));
		write_scalar_HDF5(group_killer, killer.step_killed, string("step_killed"));
		write_scalar_HDF5(group_killer, killer.Hz_steps, string("Hz_steps"));
		write_scalar_HDF5(group_killer, killer.runaway_Hz, string("runaway_Hz"));
		write_scalar_HDF5(group_killer, killer.min_steps, string("min_steps"));
		write_scalar_HDF5(group_killer, killer.min_pop_size, string("min_pop_size"));
	}

	if (jh_learn_pop.on) {
		string str = pop_n+"/jh_learn_pop/";
		Group group_jh_learn_pop = group_pop.createGroup(str);
		write_scalar_HDF5(group_jh_learn_pop,jh_learn_pop.on,string("on"));
		write_vector_HDF5(group_jh_learn_pop,jh_learn_pop.rhatE,string("rhatE"));	
		write_vector_HDF5(group_jh_learn_pop,jh_learn_pop.rhatI,string("rhatI"));	
		write_vector_HDF5(group_jh_learn_pop,jh_learn_pop.noise_spikes,string("noise_spikes"));	
		write_scalar_HDF5(group_jh_learn_pop,jh_learn_pop.noise,string("noise"));
	}
	if (spike_file.on) {
		string str = pop_n + "/spike_file/";
		Group group_spike_file = group_pop.createGroup(str);
		write_scalar_HDF5(group_spike_file, spike_file.spike_ind, string("spike_ind"));
		write_string_HDF5(group_spike_file, spike_file.file_name, string("file_name"));
	}

	if (poisson_pop.on) {
		string str = pop_n + "/poisson_pop/";
		Group group_spike_file = group_pop.createGroup(str);
		write_scalar_HDF5(group_spike_file, poisson_pop.rate, string("rate"));
		write_vector_HDF5(group_spike_file, poisson_pop.next_spike_time, string("next_spike_time"));
	}
	if (current_file.on) {
		string str = pop_n + "/current_file/";
		Group group_current_file = group_pop.createGroup(str);
		write_scalar_HDF5(group_current_file, current_file.current_ind, string("current_ind"));
		write_scalar_HDF5(group_current_file, current_file.framestep, string("framestep"));
		write_string_HDF5(group_current_file, current_file.file_name, string("file_name"));
	}

}

void NeuroPop::output_results(H5File& file) {

	// new group
	string group_name = "/pop_result_";
	group_name.append(to_string(pop_ind));
	Group group_pop = file.createGroup(group_name);

	write_vector_HDF5(group_pop, spike_hist_tot, string("spike_hist_tot"));
	write_vector_HDF5(group_pop, num_spikes_pop, string("num_spikes_pop"));
	write_vector_HDF5(group_pop, num_ref_pop, string("num_ref_pop"));

	write_string_HDF5(group_pop, dump_para(), string("pop_para"));

	if (stats.record) {
		write_vector_HDF5(group_pop, stats.V_mean, string("stats_V_mean"));
		write_vector_HDF5(group_pop, stats.V_std, string("stats_V_std"));
		write_vector_HDF5(group_pop, stats.I_input_mean, string("stats_I_input_mean"));
		write_vector_HDF5(group_pop, stats.I_input_std, string("stats_I_input_std"));
		write_vector_HDF5(group_pop, stats.I_AMPA_time_avg, string("stats_I_AMPA_time_avg"));
		write_vector_HDF5(group_pop, stats.I_NMDA_time_avg, string("stats_I_NMDA_time_avg"));
		write_vector_HDF5(group_pop, stats.I_GABA_time_avg, string("stats_I_GABA_time_avg"));
		write_vector_HDF5(group_pop, stats.I_ext_time_avg, string("stats_I_ext_time_avg"));
		write_vector_HDF5(group_pop, stats.I_tot_time_mean, string("stats_I_tot_time_mean"));
		write_vector_HDF5(group_pop, stats.I_tot_time_var, string("stats_I_tot_time_var"));
		write_vector_HDF5(group_pop, stats.V_time_mean, string("stats_V_time_mean"));
		write_vector_HDF5(group_pop, stats.V_time_var, string("stats_V_time_var"));
		write_vector_HDF5(group_pop, stats.IE_ratio, string("stats_IE_ratio"));
	}
	if (stats.record_cov) {
		write_matrix_HDF5(group_pop, stats.V_time_cov, string("stats_V_time_cov"));
	}

	if (LFP.record) {
		write_matrix_HDF5(group_pop, LFP.data, string("LFP_data"));
	}

}


void NeuroPop::add_sampling_real_time_HDF5(const vector<int>& sample_neurons_input, const vector<bool>& sample_type_input,  const vector<bool>& sample_time_points_input, string sample_file_name_input) {
	sample.file_type = 2;
	sample.neurons = sample_neurons_input;
	sample.type = sample_type_input;
	sample.time_points = sample_time_points_input;

	sample.N_steps = 0;
	for (int tt = 0; tt < step_tot; ++tt) {
		if (sample.time_points[tt]) {sample.N_steps += 1;}
	}
	sample.N_neurons = sample.neurons.size();

	sample.file_name = sample_file_name_input;
	sample.file_name.append("_");
	sample.file_name.append(to_string(pop_ind));

	sample.file_name.append("_neurosamp.h5");
	// Create a new file using default properties.
	// cout << "Creating HDF5 output file..." << sample.file_name;
	sample.file_HDF5 = new H5File( sample.file_name.c_str(), H5F_ACC_TRUNC );
	// dataset dimensions
	hsize_t dims[2];
	dims[1] = sample.N_neurons;
	dims[0] = sample.N_steps;
	// file dataspace
	DataSpace fspace(2, dims);
	// create datasets
	if (sample.type[0]) {
		sample.V_dataset = sample.file_HDF5->createDataSet("V", PredType::NATIVE_DOUBLE, fspace);
	}
	if (sample.type[1]) {
		sample.I_leak_dataset = sample.file_HDF5->createDataSet("I_leak", PredType::NATIVE_DOUBLE, fspace);
	}
	if (sample.type[2]) {
		sample.I_AMPA_dataset = sample.file_HDF5->createDataSet("I_AMPA", PredType::NATIVE_DOUBLE, fspace);
	}
	if (sample.type[3]) {
		sample.I_GABA_dataset = sample.file_HDF5->createDataSet("I_GABA", PredType::NATIVE_DOUBLE, fspace);
	}
	if (sample.type[4]) {
		sample.I_NMDA_dataset = sample.file_HDF5->createDataSet("I_NMDA", PredType::NATIVE_DOUBLE, fspace);
	}
	if (sample.type[5]) {
		sample.I_GJ_dataset = sample.file_HDF5->createDataSet("I_GJ", PredType::NATIVE_DOUBLE, fspace);
	}
	if (sample.type[6]) {
		sample.I_ext_dataset = sample.file_HDF5->createDataSet("I_ext", PredType::NATIVE_DOUBLE, fspace);
	}
	if (sample.type[7]) {
		sample.I_K_dataset = sample.file_HDF5->createDataSet("I_K", PredType::NATIVE_DOUBLE, fspace);
	}
	if(sample.type[8]){
		sample.rhatE_dataset = sample.file_HDF5->createDataSet("rhatE", PredType::NATIVE_DOUBLE, fspace);
		sample.rhatI_dataset = sample.file_HDF5->createDataSet("rhatI", PredType::NATIVE_DOUBLE, fspace);
	}
}

void NeuroPop::output_sampled_data_real_time_HDF5(const int step_current) {
	if (!sample.neurons.empty() && sample.file_type == 2) {
		if (sample.time_points[step_current]) { // push_back is amazing

			// use double *temp_data=new double[sample.N_neurons]; and later delete[] tmp_data?
			vector<double> temp_data;
			temp_data.resize(sample.N_neurons);

			int ind_temp;
			if (sample.type[0]) {	// "V"
				fill(temp_data.begin(), temp_data.end(), 0);
				// Collect the data to write in a temp vector
				for (int i = 0; i < sample.N_neurons; ++i) { // performance issue when sampling many neurons?
					ind_temp = sample.neurons[i];
					temp_data[i] = V[ind_temp];
				}
				append_vector_to_matrix_HDF5(sample.V_dataset, temp_data,  sample.ctr);
			}
			if (sample.type[1]) {	// "I_leak"
				fill(temp_data.begin(), temp_data.end(), 0);
				// Collect the data to write in a temp vector
				for (int i = 0; i < sample.N_neurons; ++i) { // performance issue when sampling many neurons?
					ind_temp = sample.neurons[i];
					temp_data[i] = I_leak[ind_temp];
				}
				append_vector_to_matrix_HDF5(sample.I_leak_dataset, temp_data,  sample.ctr);
			}
			if (sample.type[2]) {	// "I_AMPA"
				fill(temp_data.begin(), temp_data.end(), 0);
				// Collect the data to write in a temp vector
				for (int i = 0; i < sample.N_neurons; ++i) { // performance issue when sampling many neurons?
					ind_temp = sample.neurons[i];
					temp_data[i] = I_AMPA[ind_temp];
				}
				append_vector_to_matrix_HDF5(sample.I_AMPA_dataset, temp_data,  sample.ctr);
			}
			if (sample.type[3]) {	// "I_GABA"
				fill(temp_data.begin(), temp_data.end(), 0);
				// Collect the data to write in a temp vector
				for (int i = 0; i < sample.N_neurons; ++i) { // performance issue when sampling many neurons?
					ind_temp = sample.neurons[i];
					temp_data[i] = I_GABA[ind_temp];
				}
				append_vector_to_matrix_HDF5(sample.I_GABA_dataset, temp_data,  sample.ctr);
			}
			if (sample.type[4]) {	// "I_NMDA"
				fill(temp_data.begin(), temp_data.end(), 0);
				// Collect the data to write in a temp vector
				for (int i = 0; i < sample.N_neurons; ++i) { // performance issue when sampling many neurons?
					ind_temp = sample.neurons[i];
					temp_data[i] = I_NMDA[ind_temp];
				}
				append_vector_to_matrix_HDF5(sample.I_NMDA_dataset, temp_data,  sample.ctr);
			}
			if (sample.type[5]) {	// "I_GJ"
				fill(temp_data.begin(), temp_data.end(), 0);
				// Collect the data to write in a temp vector
				for (int i = 0; i < sample.N_neurons; ++i) { // performance issue when sampling many neurons?
					ind_temp = sample.neurons[i];
					temp_data[i] = I_GJ[ind_temp];
				}
				append_vector_to_matrix_HDF5(sample.I_GJ_dataset, temp_data,  sample.ctr);
			}
			if (sample.type[6]) {	// "I_ext"
				fill(temp_data.begin(), temp_data.end(), 0);
				// Collect the data to write in a temp vector
				for (int i = 0; i < sample.N_neurons; ++i) { // performance issue when sampling many neurons?
					ind_temp = sample.neurons[i];
					temp_data[i] = I_ext[ind_temp];
				}
				append_vector_to_matrix_HDF5(sample.I_ext_dataset, temp_data,  sample.ctr);
			}
			if (sample.type[7]) {	// "I_K"
				fill(temp_data.begin(), temp_data.end(), 0);
				// Collect the data to write in a temp vector
				for (int i = 0; i < sample.N_neurons; ++i) { // performance issue when sampling many neurons?
					ind_temp = sample.neurons[i];
					temp_data[i] = I_K[ind_temp];
				}
				append_vector_to_matrix_HDF5(sample.I_K_dataset, temp_data,  sample.ctr);
			}
			if (sample.type[8]){	// "rhat"
				//rhat
				fill(temp_data.begin(), temp_data.end(), 0);
				// Collect the data to write in a temp vector
				for (int i = 0; i < sample.N_neurons; ++i){ // performance issue when sampling many neurons?
					ind_temp = sample.neurons[i];
					temp_data[i]= jh_learn_pop.rhatE[ind_temp];
				}
				append_vector_to_matrix_HDF5(sample.rhatE_dataset, temp_data,  sample.ctr);

				fill(temp_data.begin(), temp_data.end(), 0);
				// Collect the data to write in a temp vector
				for (int i = 0; i < sample.N_neurons; ++i){ // performance issue when sampling many neurons?
					ind_temp = sample.neurons[i];
					temp_data[i]= jh_learn_pop.rhatI[ind_temp];
				}
				append_vector_to_matrix_HDF5(sample.rhatI_dataset, temp_data,  sample.ctr);

			}	
			sample.ctr++;
		}
	}
}


