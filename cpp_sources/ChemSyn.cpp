/// \file
//#include <iostream>
#include <cmath> // always cmath instead of math.h ?
#include <algorithm>
#include <sstream>
//#include <fstream>

#include "ChemSyn.h"

ChemSyn::ChemSyn(const double dt_input, const int step_tot_input){
	
	dt = dt_input;
	step_tot = step_tot_input;
	
	// Default parameters
	V_ex = 0.0;     // Excitatory reversal, 0.0
	V_in = -80.0;   // Inhibitory reversal, -80.0

	//  time-evolution of post-synaptic conductance change (msec)
	Dt_trans_AMPA = 1.0; // 0.5
	Dt_trans_GABA = 1.0; // 1.0
	Dt_trans_NMDA = 5.0; // 5.0
	tau_decay_AMPA = 5.0; // 3.0
	tau_decay_GABA = 3.0; // 7.0
	tau_decay_NMDA = 80.0; // 80.0

	// default settting
	stats.record = false;
	stats.record_cov = false;
	STD.on = false; 
	STD.on_step = -1;
	inh_STDP.on_step = -1;
	inh_STDP.on = false;
	synapse_model = 0; // default model
	
}



void ChemSyn::init(const int syn_type_input, const int i_pre, const int j_post, const int N_pre_input, const int N_post_input, const vector<int> &C_i, const vector<int> &C_j, const vector<double> &K_ij, const vector<double> &D_ij){

	// input check
	if (*max_element(C_i.begin(), C_i.end()) >= N_pre_input){
		cout << "Wrong Input!In the chemical synapses " << syn_type_input + 1 <<" connection from pop "<< i_pre + 1 << " to pop " << j_post + 1 << ",the index of pre_population synapses exceed the number of pre-population neurons." << endl;
	}
				 
	if (*max_element(C_j.begin(), C_j.end()) >= N_post_input){
		cout << "Wrong Input!In the chemical synapses " << syn_type_input + 1 <<" connection from pop "<< i_pre + 1 << " to pop " << j_post + 1 << ",the index of post_population synapses exceed the number of post-population neurons." << endl;
	}

	// read parameter
	syn_type = syn_type_input;
	pop_ind_pre = i_pre;
	pop_ind_post = j_post;
	N_pre = N_pre_input;
	N_post = N_post_input;
	double max_delay = *max_element(D_ij.begin(), D_ij.end()); // max_element returns an iterator 
	max_delay_steps = int(round(max_delay / dt)) ;
	
	
	// read in C, K, D
	// Initialise s_TALS, s_VALS
	int i_temp, j_temp;
	for (unsigned int ind = 0; ind < K_ij.size(); ++ind){
		if (K_ij[ind] >= 0.0){ // must be no less than zero! unit: miuSiemens
			if (K.empty()){ // If empty, initialise them
				C.resize(N_pre);
				K.resize(N_pre);
				D.resize(N_pre);
			}
			i_temp = C_i[ind];
			j_temp = C_j[ind];
			C[i_temp].push_back(j_temp);
			K[i_temp].push_back(K_ij[ind]);
			
			D[i_temp].push_back((int)round(D_ij[ind] / dt)); // note that D_ij is in msec
		}
		// discard all the zeros
		else{ continue; }
	}


	
	// parameter-dependent initialisation
	init();

}


void ChemSyn::init(const int syn_type_input, const int j_post, const int N_post_input, const double K_ext, const int Num_ext, const vector<double> &rate_ext_t, const vector<bool> &neurons){

	// Initialise chemical synapses for simulating external neuron population
	syn_type = syn_type_input;
	pop_ind_pre = -1; // -1 for external noisy population
	pop_ind_post = j_post;
	N_pre = 1; // just for initialization
	N_post = N_post_input;
	max_delay_steps = 0; // no delay;


	// Parameters for noise generation
	ext_noise.K_ext = K_ext;
	ext_noise.Num_ext = Num_ext;
	ext_noise.rate_ext_t = rate_ext_t;	
	ext_noise.neurons = neurons;
	
	// Random seed (random engine should be feed with DIFFERENT seed at every implementation)
	random_device rd; // random number from operating system for seed
	my_seed = rd(); // record seed
	//my_seed = 321;
	//cout << "My_seed is: " << my_seed << endl;
	
	// parameter-dependent initialisation
	init();
}


void ChemSyn::init(const int syn_type_input, const int j_post, const int N_post_input, const double K_ext, const int Num_ext, const vector<bool> &rate_ext_on, const vector<double> &rate_ext_neuron){

	// Initialise chemical synapses for simulating external neuron population
	syn_type = syn_type_input;
	pop_ind_pre = -2; // -2 for external noisy population with different time-invariant rates for each neuron
	pop_ind_post = j_post;
	N_pre = 1; // just for initialization
	N_post = N_post_input;
	max_delay_steps = 0; // no delay;


	// Parameters for noise generation
	ext_noise_t_inv.K_ext = K_ext;
	ext_noise_t_inv.Num_ext = Num_ext;
	ext_noise_t_inv.rate_ext_neuron = rate_ext_neuron; // time-invariant rate for each neuron
	ext_noise_t_inv.rate_ext_on = rate_ext_on;
	for (int j = 0; j < N_post_input; ++j){
		ext_noise_t_inv.poi_dist.push_back(new poisson_distribution<int>(ext_noise_t_inv.Num_ext * ext_noise_t_inv.rate_ext_neuron[j] * (dt / 1000.0) ) );	
	}


	// Random seed (random engine should be feed with DIFFERENT seed at every implementation)
	random_device rd; // random number from operating system for seed
	my_seed = rd(); // record seed
	gen.seed(my_seed);

	// parameter-dependent initialisation
	init();
}




void ChemSyn::set_seed(int seed_input){
	if (pop_ind_pre != -1){
		cout << "Cannot set RNG seed to the synapses with type " << syn_type + 1 <<" from pop "<< pop_ind_pre + 1 << " to pop " << pop_ind_post + 1 << " since this object does not have RNG." << endl;
	}
	else{
		my_seed = seed_input; // This will overwrite the auto generated seed.
		gen.seed(my_seed);
	}
}

void ChemSyn::init(){
	// parameter-dependent initialisation

	
	// Initialise chemical synapse parameters
	if (syn_type == 0){
		tau_decay = tau_decay_AMPA;
		tau_rise = Dt_trans_AMPA;
	}
	else if (syn_type == 1){	
		tau_decay = tau_decay_GABA;
		tau_rise = Dt_trans_GABA;
	}
	else if (syn_type == 2){
		tau_decay = tau_decay_NMDA;
		tau_rise = Dt_trans_NMDA;
		// non-linearity of NMDA
		// voltage-dependent part B(V) (look-up table):
		miuMg_NMDA = 0.33; // mM^-1, concentration of [Mg2+] is around 1 mM, 0.33
		gamma_NMDA = 0.06; // mV^-1, 0.06
		B_V_min = -80.0 - 1.0; // < V_in = -80, check if they are consistent!!
		B_V_max = -55.0 + 1.0; // > V_th = -55
		B_dV = 0.1; // 0.1
		int i_B = 0;
		double V_temp, B_temp;
		B.resize(0); // for re-initialization!!! 
		while (true){
			V_temp = B_V_min + i_B*B_dV;
			if (V_temp > B_V_max){ break; }
			B_temp = 1 / (1 + miuMg_NMDA*exp(-gamma_NMDA*V_temp));
			B.push_back(B_temp);
			i_B += 1;
		}
	}
	steps_trans = int(round(tau_rise / dt));

	// Initialise exp_step
	exp_step_decay = exp(-dt / tau_decay); // single step
	exp_step_rise = exp(-dt / tau_rise);
	I.assign(N_post, 0);
	gs_sum.assign(N_post, 0);
	
	// transmitter_strength
	K_trans.assign(N_pre, 1.0 / steps_trans); // be careful! 1 / transmitter steps gives zero (int)!!
	
	// Synapse model choice
	// model 0, the default model
	if (synapse_model == 0){ 
		// Initialize pre- and post-synaptic dynamic variables
		gsm_0.buffer_steps = max_delay_steps + steps_trans + 1; // the +1 is vital
		gsm_0.s.assign(N_pre, 0);
		gsm_0.d_gs_sum_buffer.resize(gsm_0.buffer_steps);
		for (int i = 0; i < gsm_0.buffer_steps; ++i){
			gsm_0.d_gs_sum_buffer[i].assign(N_post, 0);
		}
		gsm_0.trans_left.assign(N_pre, 0);
	}
	else if (synapse_model == 1){
		// model 1
		gsm_1.buffer_steps = max_delay_steps + 1; // the +1 is vital
		gsm_1.gs_rise_sum.assign(N_post, 0);
		gsm_1.gs_decay_sum.assign(N_post, 0);
		gsm_1.d_gs_rd_sum_buffer.resize(gsm_1.buffer_steps);
		for (int i = 0; i < gsm_1.buffer_steps; ++i){
			gsm_1.d_gs_rd_sum_buffer[i].assign(N_post, 0);
		}
		// clear model 0
		gsm_0.s.clear();
		gsm_0.d_gs_sum_buffer.clear(); // clear() clears all of its components recursively
		gsm_0.trans_left.clear();
	}
}

const int & ChemSyn::get_direction()
{
	return jh_learn_syn.direction;
}

const int & ChemSyn::get_syn_type()
{
	return syn_type;
}
const int & ChemSyn::get_pop_ind_pre()
{
	return pop_ind_pre;
}
const int & ChemSyn::get_pop_ind_post()
{
	return pop_ind_post;
}

const int & ChemSyn::get_post_neu_type()
{
	return jh_learn_syn.ntype_post;
}

void ChemSyn::set_synapse_model(const int synapse_model_input){
	if (synapse_model_input != 0){
		synapse_model = synapse_model_input;
		init(); // initialize again
	}
}

	

void ChemSyn::update(const int step_current){

	if (synapse_model == 0){ 
		// short-term depression
		update_STD(step_current);
		
		// inhibitory STDP
		update_inh_STDP(step_current);
		
		// Update transmitter dynamics
		update_gs_sum_model_0(step_current);
	}
	else if (synapse_model == 1){
		update_gs_sum_model_1(step_current);
	}
	
	// Calculate chemical currents
	calc_I();

	// sample data
	sample_data(step_current);
	
	//
	record_stats(step_current);
	
}// update



void ChemSyn::update_STD(const int step_current){
	// STD modifies K_trans
	
	if (STD.on_step == step_current){STD.on = true;}
	// short-term depression
	if (STD.on == true){
		for (unsigned int i = 0; i < spikes_pre.size(); ++i){ 
			K_trans[spikes_pre.at(i)] = 1.0 / steps_trans * STD.f_ves[i];
			STD.f_ves[spikes_pre.at(i)] *= 1.0 - STD.p_ves; // decrease at spikes
		}
		for (int i = 0; i < N_pre; ++i){
			STD.f_ves[i] = 1.0 - STD.exp_ves * (1.0 - STD.f_ves[i]); // decay to 1.0
		}
	}

}


void ChemSyn::add_inh_STDP(const int inh_STDP_on_step_input){
	if (syn_type != 1){
		cout << "Warning: initializing inhibitory STDP on non-GABA synapses!" << endl;
	}
	if (synapse_model != 0){
		cout << "Warning: inhibitory STDP is not supported for synapse models other than 0 yet!" << endl;
	}
	
	inh_STDP.on_step = inh_STDP_on_step_input;
	if (inh_STDP.on_step == 0){
		inh_STDP.on = true;
	}
	
	inh_STDP.x_trace_pre.assign(N_pre, 0.0);
	inh_STDP.x_trace_post.assign(N_post, 0.0);
	
	inh_STDP.tau = 20; // ms
	inh_STDP.exp_step = exp(-dt / inh_STDP.tau);
	inh_STDP.eta = 0.0001; // learning rate, 0.0001 is the published value but requires 60min of simulation
	inh_STDP.rho_0 = 0.010; //kHz
	inh_STDP.alpha = 2.0 * inh_STDP.rho_0 * inh_STDP.tau; // depression factor

	// j_2_i and j_2_syn_ind
	inh_STDP.j_2_i.resize(N_post);
	inh_STDP.j_2_syn_ind.resize(N_post);
	int j_post;
	for (int i_pre = 0; i_pre < N_pre; ++i_pre){ 
		for (unsigned int syn_ind = 0; syn_ind < C[i_pre].size(); ++syn_ind){
			j_post = C[i_pre][syn_ind];
			inh_STDP.j_2_i[j_post].push_back( i_pre );
			inh_STDP.j_2_syn_ind[j_post].push_back( int(syn_ind) );
		}
	}
}	

void ChemSyn::calc_I(){
	// need V from post population
	// need gs_sum from previous calculations
	
	if (syn_type == 0){ //AMPA
		for (int j = 0; j < N_post; ++j){
			I[j] = -gs_sum[j] * (V_post.at(j) - V_ex);
		}
	}	
	else if (syn_type == 1){ //GABA
		for (int j = 0; j < N_post; ++j){
			I[j] = -gs_sum[j] * (V_post.at(j) - V_in);
			// For inhibition, every equation is in the same form as excitation. 
			// Only "V_in" encodes its inhibitory nature.
		}
	}
	else if (syn_type == 2){ //NMDA
		for (int j = 0; j < N_post; ++j){
			I[j] = -gs_sum[j] * B[(int)round((V_post.at(j) - B_V_min) / B_dV)] * (V_post.at(j) - V_ex);
		}
	}

	
	
}

void ChemSyn::update_gs_sum_model_0(const int step_current){
	// See Gu, Yifan, Gong, Pulin, 2016, The dynamics of memory retrieval in hierarchical networks: a modeling study
	// this function updates gs_sum
	if (pop_ind_pre >= 0){
		for (unsigned int i = 0; i < spikes_pre.size(); ++i){ // add spikes (transmitter release)
			gsm_0.trans_left[spikes_pre.at(i)] += steps_trans;
		}
		for (int i_pre = 0; i_pre < N_pre; ++i_pre){
			if (gsm_0.trans_left[i_pre] > 0){
				// conduction delay: put into d_gs_sum_buffer
				int j_post, delay_step, t_ring;
				for (unsigned int syn_ind = 0; syn_ind < C[i_pre].size(); ++syn_ind){ //loop through all the post-synapses
					j_post = C[i_pre][syn_ind]; // index of the post-synaptic neuron
					delay_step = D[i_pre][syn_ind]; // delay in steps for this post-synaptic neuron
					t_ring = int( (step_current + delay_step) % gsm_0.buffer_steps ); // index in the gs_buffer
					gsm_0.d_gs_sum_buffer[t_ring][j_post] += K_trans[i_pre] * (1.0 - gsm_0.s[i_pre]) * K[i_pre][syn_ind];
				}
				gsm_0.trans_left[i_pre] -= 1;
				gsm_0.s[i_pre] += K_trans[i_pre] * (1.0 - gsm_0.s[i_pre]);
			}
		}
		// decay pre-synaptic dynamics
		for (int i = 0; i < N_pre; ++i){ gsm_0.s[i] *= exp_step_decay; };
	}
	else if (pop_ind_pre == -1){ // if external noisy population
		// Contribution of external spikes, assuming square pulse transmitter release
		// Generate current random number generator, note that rate_ext_t is in Hz
		poisson_distribution<int> dist(ext_noise.Num_ext * ext_noise.rate_ext_t[step_current] * (dt / 1000.0));		
		auto ext_spikes = bind(dist, ref(gen));

		// Post-synaptic dynamics
		int t_ring;
		for (int t_trans = 0; t_trans < steps_trans; ++t_trans){
			t_ring = int( (step_current + t_trans) % gsm_0.buffer_steps );
			for (int j_post = 0; j_post < N_post; ++j_post){
				if (ext_noise.neurons[j_post]){
					gsm_0.d_gs_sum_buffer[t_ring][j_post] += K_trans[0] * ext_noise.K_ext * ext_spikes(); 
				}
			}
		}
	}
	else if (pop_ind_pre == -2){ // if external noisy population
		// Contribution of external spikes, assuming square pulse transmitter release
		// Generate current random number generator, note that rate_ext_t is in Hz

		// Post-synaptic dynamics
		if (ext_noise_t_inv.rate_ext_on[step_current]){
			int t_ring;
			for (int t_trans = 0; t_trans < steps_trans; ++t_trans){
				t_ring = int( (step_current + t_trans) % gsm_0.buffer_steps );
				for (int j_post = 0; j_post < N_post; ++j_post){
					gsm_0.d_gs_sum_buffer[t_ring][j_post] += K_trans[0] * ext_noise_t_inv.K_ext * ext_noise_t_inv.poi_dist[j_post]->operator()(gen); 
				}
			}
		}
	}
	
	
	// update post-synaptic dynamics
	int t_ring = int( step_current % gsm_0.buffer_steps );
	for (int j_post = 0; j_post < N_post; ++j_post){
		gs_sum[j_post] += gsm_0.d_gs_sum_buffer[t_ring][j_post];
		// should I decay gs_sum here??
		// decay gs_sum
		// numerical error of this integration scheme should be less 1.5%
		gs_sum[j_post] *= exp_step_decay;
	}
	// immediately reset the current buffer to zeros after being used!!
	fill(gsm_0.d_gs_sum_buffer[t_ring].begin(), gsm_0.d_gs_sum_buffer[t_ring].end(), 0.0);
}

void ChemSyn::update_gs_sum_model_1(const int step_current){
	// See Keane, A., Gong, P., 2015, Propagating Waves Can Explain Irregular Neural Dynamics
	// this function updates gs_sum
	for (unsigned int ind = 0; ind < spikes_pre.size(); ++ind){ // loop through all the spikes
		int i_pre = spikes_pre.at(ind);
		int j_post, delay_step, t_ring;
		for (unsigned int syn_ind = 0; syn_ind < C[i_pre].size(); ++syn_ind){ //loop through all the post-synapses
			j_post = C[i_pre][syn_ind]; // index of the post-synaptic neuron
			delay_step = D[i_pre][syn_ind]; // delay in steps for this post-synaptic neuron
			t_ring = int( (step_current + delay_step) % gsm_1.buffer_steps ); // index in the gs_buffer
			gsm_1.d_gs_rd_sum_buffer[t_ring][j_post] += K[i_pre][syn_ind];  // the peak value is linear to the initial impulse. 
		}
	}
	int t_ring = int( step_current % gsm_1.buffer_steps );
	for (int j_post = 0; j_post < N_post; ++j_post){ // Check the error of the following numerical scheme!
		gsm_1.gs_rise_sum[j_post] *= exp_step_rise;
		gsm_1.gs_decay_sum[j_post] *= exp_step_decay;
		gs_sum[j_post] = (gsm_1.gs_decay_sum[j_post] - gsm_1.gs_rise_sum[j_post]) / (tau_decay - tau_rise); 
		gsm_1.gs_rise_sum[j_post] += gsm_1.d_gs_rd_sum_buffer[t_ring][j_post];
		gsm_1.gs_decay_sum[j_post] += gsm_1.d_gs_rd_sum_buffer[t_ring][j_post];
	}
	// immediately reset the current buffer to zeros after being used!!
	fill(gsm_1.d_gs_rd_sum_buffer[t_ring].begin(), gsm_1.d_gs_rd_sum_buffer[t_ring].end(), 0.0);
}


void ChemSyn::add_sampling(const vector<int> & sample_neurons_input, const vector<bool> & sample_time_points_input){
	sample.neurons = sample_neurons_input;
	sample.time_points = sample_time_points_input;
	
	// initialise
	int sample_time_points_tot = 0;// count non zero elements in sample_time_points
	for (unsigned int i = 0; i < sample.time_points.size(); ++i){
		if (sample.time_points[i]){
			sample_time_points_tot += 1;
		}
	}
	int sample_neurons_tot = sample.neurons.size();// count non zero elements in sample_time_points

	sample.data.resize(sample_neurons_tot);
	for (int i = 0; i < sample_neurons_tot; ++i){
		sample.data[i].reserve(sample_time_points_tot); // reserve and push_back so that it won't be affected by adapting step_tot
	}

}

void ChemSyn::add_short_term_depression(const int STD_on_step_input){
	if (syn_type != 0){
		cout << "Warning: initializing STD on non-AMPA synapses!" << endl;
	}
	if (synapse_model != 0){
		cout << "Warning: STD is not supported for synapse models other than 0 yet!" << endl;
	}
	
	STD.on_step = STD_on_step_input;
	if (STD.on_step == 0){
		STD.on = true;
	}
	// short term depression
	STD.p_ves =  0.4; // see X. Wang, 1999, The Journal of Neuroscience
	STD.tau_ves =  700; // ms
	STD.f_ves.assign(N_pre, 1.0);
	STD.exp_ves = exp(-dt / STD.tau_ves);
}




void ChemSyn::update_inh_STDP(const int step_current){
	// inh_STDP modifies K
	
	if (inh_STDP.on_step == step_current){inh_STDP.on = true;}
	if (inh_STDP.on == true){
		// update x_trace
		for (unsigned int i = 0; i < spikes_pre.size(); ++i){
			inh_STDP.x_trace_pre[spikes_pre.at(i)] += 1.0;
		}
		for (unsigned int j = 0; j < spikes_post.size(); ++j){
			inh_STDP.x_trace_post[spikes_post.at(j)] += 1.0;
		}
		// update K
		int i_pre, j_post;
		for (unsigned int ind_spike = 0; ind_spike < spikes_pre.size(); ++ind_spike){
			i_pre = spikes_pre.at(ind_spike);
			for (unsigned int syn_ind = 0; syn_ind < C[i_pre].size(); ++syn_ind){
				j_post = C[i_pre][syn_ind];
				K[i_pre][syn_ind] += inh_STDP.eta * ( inh_STDP.x_trace_post[j_post] - inh_STDP.alpha );
			}
		}
		int syn_ind;
		for (unsigned int ind_spike = 0; ind_spike < spikes_post.size(); ++ind_spike){
			j_post = spikes_post.at(ind_spike);
			for (unsigned int ind = 0; ind < inh_STDP.j_2_i[j_post].size(); ++ind){
				// j_2_i and j_2_syn_ind together serve as the "inverse function" of j_post = C[i_pre][syn_ind]
				i_pre = inh_STDP.j_2_i[j_post][ind];
				syn_ind = inh_STDP.j_2_syn_ind[j_post][ind];
				K[i_pre][syn_ind] += inh_STDP.eta * inh_STDP.x_trace_pre[i_pre];
			}
		}
		// for testing
		//if (tmp_data.size() == 0){
		//	tmp_data.resize(2);
		//}
		//tmp_data[0].push_back(K[0][0]);
		//tmp_data[1].push_back(K[100][0]);
		for (int i = 0; i < N_pre; ++i){ inh_STDP.x_trace_pre[i] *= inh_STDP.exp_step; }
		for (int j = 0; j < N_post; ++j){ inh_STDP.x_trace_post[j] *= inh_STDP.exp_step; }
	}
}


void ChemSyn::sample_data(const int step_current){
	if (!sample.neurons.empty()){
		if (sample.time_points[step_current]){ // push_back is amazing
			for (unsigned int i = 0; i < sample.neurons.size(); ++i){ // performance issue when sampling many neurons?
				int ind_temp = sample.neurons[i];
				sample.data[i].push_back( I[ind_temp] );
			}
		}
	}
}



void ChemSyn::set_para(string para_str){
	const char delim = ',';
	if (!para_str.empty()){
		istringstream para(para_str);
		string para_name, para_value_str; 
		double para_value;
		while (getline(para, para_name, delim)){
			getline(para, para_value_str, delim); // get parameter value (assume double)
			stringstream(para_value_str) >> para_value; // from string to numerical value
			if (para_name.find("V_ex") != string::npos){V_ex = para_value;}
			else if (para_name.find("V_in") != string::npos){V_in = para_value;}
			else if (para_name.find("Dt_trans_AMPA") != string::npos){Dt_trans_AMPA = para_value;}
			else if (para_name.find("Dt_trans_GABA") != string::npos){Dt_trans_GABA = para_value;}
			else if (para_name.find("Dt_trans_NMDA") != string::npos){Dt_trans_NMDA = para_value;}
			else if (para_name.find("tau_decay_AMPA") != string::npos){tau_decay_AMPA = para_value;}
			else if (para_name.find("tau_decay_GABA") != string::npos){tau_decay_GABA = para_value;}
			else if (para_name.find("tau_decay_NMDA") != string::npos){tau_decay_NMDA = para_value;}
			else {cout << "Unrecognized parameter: " << para_name << endl;}
		}
	}
	// re-initialise it!
	init();
}


string ChemSyn::dump_para(){
	const char delim = ',';
	
	stringstream dump;

	dump << "pop_ind_pre" << delim << pop_ind_pre << delim << endl;
	dump << "pop_ind_post" << delim << pop_ind_post << delim << endl;
	dump << "syn_type" << delim << syn_type << delim << endl;

	dump << "V_ex" << delim << V_ex << delim << endl;
	dump << "V_in" << delim << V_in << delim << endl;

	dump << "seed" << delim << my_seed << delim << endl;

	dump << "synapse_model" << delim << synapse_model << endl;
	
	if (syn_type == 0){
		dump << "Dt_trans_AMPA" << delim << Dt_trans_AMPA << delim << endl;
		dump << "tau_decay_AMPA" << delim << tau_decay_AMPA << delim << endl;
	}
	else if (syn_type == 1){
		dump << "Dt_trans_GABA" << delim << Dt_trans_GABA << delim << endl;
		dump << "tau_decay_GABA" << delim << tau_decay_GABA << delim << endl;
	}
	else if (syn_type == 2){
		dump << "Dt_trans_NMDA" << delim << Dt_trans_NMDA << delim << endl;
		dump << "tau_decay_NMDA" << delim << tau_decay_NMDA << delim << endl;
	}

		
	return dump.str();
}

void ChemSyn::start_stats_record(){
	stats.record = true;
	stats.I_mean.reserve(step_tot);
	stats.I_std.reserve(step_tot);
	if (synapse_model == 0){
		stats.s_time_mean.assign(N_pre, 0.0);
		stats.s_time_var.assign(N_pre, 0.0);
		stats.I_time_mean.assign(N_pre, 0.0);
		stats.I_time_var.assign(N_pre,0.0);
	}
}

void ChemSyn::start_cov_record(const int time_start, const int time_end){
	if (synapse_model == 0){
		stats.record_cov = true;
		stats.time_start_cov = time_start;
		stats.time_end_cov = time_end;
		stats.s_time_mean_dumb.assign(N_pre, 0.0);
		stats.s_time_cov.resize(N_pre);
		for (int i = 0;i < N_pre; ++i){
			stats.s_time_cov[i].assign(N_pre,0.0);
		}
	}
}

void ChemSyn::add_JH_Learning(vector<NeuroPop*> &NeuronPopArray,int isteps, double iscale,double lrate, double lrateall,int intau, double innoise,int type_pre,int type_post, int direction){
	jh_learn_syn.on=true; //indicates this learning is to be used
	jh_learn_syn.direction=direction;

	jh_learn_syn.noise_pre=NeuronPopArray[pop_ind_pre]->get_noise();
	jh_learn_syn.noise_post=NeuronPopArray[pop_ind_post]->get_noise();

	jh_learn_syn.inf_steps=isteps;
	double Vinit=-70;
	jh_learn_syn.post_V_hist.resize(isteps+1, vector<double> (N_post,Vinit));  //+++TODO UPDATE init value to something else???
	jh_learn_syn.post_R_hist.resize(isteps+1, vector<int> (N_post,0)); 
	jh_learn_syn.t_ind=isteps+1; // will start by incrementing to 0
	jh_learn_syn.post_hist_len=isteps +2; //+1 for ring buffer empty space 
	jh_learn_syn.post_t_hist.resize(isteps+2,vector<int>(N_post,0)); 
	jh_learn_syn.post_noise_hist.resize(isteps+2,vector<int>(N_post,0)); 
	jh_learn_syn.ind_post_new.resize(N_post,isteps+1);
	jh_learn_syn.ind_post_old.resize(N_post,0);
	jh_learn_syn.pre_hist_len=isteps+2; //+1 for storing 0 to -inf_Steps, +1 for ring buffer empty space
	jh_learn_syn.pre_t_hist.resize(isteps+2,vector<int>(N_pre,0)); 
	jh_learn_syn.pre_noise_hist.resize(isteps+2,vector<int>(N_pre,0)); 
	jh_learn_syn.ind_pre_new.resize(N_pre,isteps+1);
	jh_learn_syn.ind_pre_old.resize(N_pre,0);
	jh_learn_syn.Vint.resize(N_post,0);
	jh_learn_syn.Vint_ctr.resize(N_post,intau);
	if(direction==1){
		//TODO
	}
	jh_learn_syn.rhat.resize(N_pre,0);
	jh_learn_syn.inf_scale=iscale;
	jh_learn_syn.learn_rate=lrate;
	jh_learn_syn.learn_rate_all=lrateall;
	jh_learn_syn.tau=intau;
	jh_learn_syn.C=NeuronPopArray[pop_ind_post]->get_Cm();
	jh_learn_syn.noise=innoise;
	jh_learn_syn.rhat.resize(N_pre,0);
	jh_learn_syn.ntype_pre=type_pre;
	jh_learn_syn.ntype_post=type_post;

	jh_learn_syn.j_2_i.resize(N_post);
	jh_learn_syn.j_2_syn_ind.resize(N_post);
	int j_post;
	for (int i_pre = 0; i_pre < N_pre; i_pre++){ 
		for (unsigned int syn_ind = 0; syn_ind < C[i_pre].size(); syn_ind++){
			j_post = C[i_pre][syn_ind];
			jh_learn_syn.j_2_i[j_post].push_back( i_pre );
			jh_learn_syn.j_2_syn_ind[j_post].push_back( int(syn_ind) );
		}
	}
}


void ChemSyn::record_V_post_JH_Learn(vector<NeuroPop*> &NeuronPopArray){
	if(jh_learn_syn.on){
		jh_learn_syn.t_ind=(jh_learn_syn.t_ind+1)%jh_learn_syn.post_V_hist.size();
		jh_learn_syn.post_V_hist[jh_learn_syn.t_ind]=V_post;
		jh_learn_syn.post_R_hist[jh_learn_syn.t_ind]=NeuronPopArray[pop_ind_post]->get_ref_step_left();
		if(jh_learn_syn.direction==1){
			jh_learn_syn.post_V_hist_b4_th[jh_learn_syn.t_ind]=V_post;
			for(unsigned int i=0;i<spikes_post.size();i++){
				jh_learn_syn.post_V_hist_b4_th[jh_learn_syn.t_ind][spikes_post[i]]=jh_learn_syn.V_th;
			}
		}
	}
}

void ChemSyn::update_post_spike_hist_JH_Learn(){
	if(jh_learn_syn.on){
		int j_post,i;
		for(j_post=0;j_post<N_post;j_post++){
			for(i=jh_learn_syn.ind_post_old[j_post];i!=(jh_learn_syn.ind_post_new[j_post]+1)%jh_learn_syn.post_hist_len;i=(i+1)%jh_learn_syn.post_hist_len){
				jh_learn_syn.post_t_hist[i][j_post]++;
				if(jh_learn_syn.post_t_hist[i][j_post]>=jh_learn_syn.inf_steps){
					//if spike is too old increment index past it
					jh_learn_syn.ind_post_old[j_post]=(jh_learn_syn.ind_post_old[j_post]+1)%jh_learn_syn.post_hist_len;
				}
			}
		}
	}
}


void ChemSyn::update_Vint_JH_Learn(){
	if(jh_learn_syn.on){
		int j_post;
		for(j_post=0;j_post<N_post;j_post++){
			jh_learn_syn.Vint[j_post]*=exp(-1/jh_learn_syn.tau); //dt=1 since tau is in units of timesteps
			jh_learn_syn.Vint_ctr[j_post]++;
			if(syn_type==0){
				if(V_post[j_post]<V_ex){
					jh_learn_syn.Vint[j_post]-=(V_post[j_post]-V_ex);
				}
			}
			else{
				if(V_post[j_post]>V_in){
					jh_learn_syn.Vint[j_post]+=(V_post[j_post]-V_in);
				}
			}
		}
	}
}


void ChemSyn::new_post_spikes_JH_Learn(){
	if(jh_learn_syn.on){
		int  j_post;
		unsigned int i;
		//update new post spikes
		for(i=0; i<spikes_post.size(); i++){
			j_post=spikes_post[i];
			jh_learn_syn.ind_post_new[j_post]=(jh_learn_syn.ind_post_new[j_post]+1)%jh_learn_syn.post_hist_len;
			jh_learn_syn.post_t_hist[jh_learn_syn.ind_post_new[j_post]][j_post]=0;
			jh_learn_syn.post_noise_hist[jh_learn_syn.ind_post_new[j_post]][j_post]=jh_learn_syn.spikes_post_noise[i];
		}
	}
}

void ChemSyn::wchange_non_Hebbian_outgoing(vector<NeuroPop*> &NeuronPopArray){
	if(jh_learn_syn.on){
		int  j_post, i_pre;
		unsigned int syn_ind,ind, i;
		double Kscale=0;
		//update new post spikes
		for(i=0; i<spikes_post.size(); i++){
			j_post=spikes_post[i];
			if(jh_learn_syn.spikes_post_noise[i]==0){
				//update weights with Vint
				Kscale=jh_learn_syn.learn_rate*(1.0/jh_learn_syn.C)*dt*jh_learn_syn.Vint[j_post]/(1.0-jh_learn_syn.noise_post);
				for( ind=0; ind< jh_learn_syn.j_2_i[j_post].size(); ind++){
					i_pre=jh_learn_syn.j_2_i[j_post][ind];
					syn_ind=jh_learn_syn.j_2_syn_ind[j_post][ind];
					K[i_pre][syn_ind]-=Kscale*K[i_pre][syn_ind];
									

					if(K[i_pre][syn_ind]<0){
						K[i_pre][syn_ind]=0;
					}
				}
			}
			jh_learn_syn.Vint[j_post]=0; //reset Vint
			jh_learn_syn.Vint_ctr[j_post]=0; //reset Vint
		}	
		if(jh_learn_syn.learn_rate_all>0){
			if(jh_learn_syn.spikes_post_noise[i]==0){
				//also modify weights on nonspiking neurons
				for(j_post=0; j_post<N_post; j_post++){
					jh_learn_syn.ind_post_new[j_post]=(jh_learn_syn.ind_post_new[j_post]+1)%jh_learn_syn.post_hist_len;
					jh_learn_syn.post_t_hist[jh_learn_syn.ind_post_new[j_post]][j_post]=0;
					
					//update weights with Vint
					for( ind=0; ind< jh_learn_syn.j_2_i[j_post].size(); ind++){
						i_pre=jh_learn_syn.j_2_i[j_post][ind];
						syn_ind=jh_learn_syn.j_2_syn_ind[j_post][ind];
						
						K[i_pre][syn_ind]-=jh_learn_syn.learn_rate_all;

						if(K[i_pre][syn_ind]<0){
							K[i_pre][syn_ind]=0;
						}
					}
				}
			}
		}
	}
}

void ChemSyn::new_pre_spikes_JH_Learn(){
	if(jh_learn_syn.on){
		int  i_pre;
		unsigned int i;
		//update new post spikes
		for(i=0; i<spikes_pre.size(); i++){
			i_pre=spikes_pre[i];
			jh_learn_syn.ind_pre_new[i_pre]=(jh_learn_syn.ind_pre_new[i_pre]+1)%jh_learn_syn.pre_hist_len;
			jh_learn_syn.pre_t_hist[jh_learn_syn.ind_pre_new[i_pre]][i_pre]=0;
			jh_learn_syn.pre_noise_hist[jh_learn_syn.ind_pre_new[i_pre]][i_pre]=jh_learn_syn.spikes_pre_noise[i];
		}
	}
}

const vector<double> & ChemSyn::get_all_rhat_JH_Learn(vector<int> neu_pre_samp){
	if(jh_learn_syn.on){
		int i_pre,j_post,inf_t_ind,age;
		unsigned int syn_ind;
		double expdecay;
		double post_h=0;
		for(unsigned int i=0;i<neu_pre_samp.size();i++){
			i_pre=neu_pre_samp[i];
			// Calculate rhat value
			jh_learn_syn.rhat[i_pre]=0;

			for (syn_ind = 0; syn_ind < C[i_pre].size(); ++syn_ind){
				j_post = C[i_pre][syn_ind];
				if(jh_learn_syn.ind_post_old[j_post]!=(jh_learn_syn.ind_post_new[j_post]+1)%jh_learn_syn.post_hist_len){
					//there is a j_post spike in history
									
					age=jh_learn_syn.post_t_hist[jh_learn_syn.ind_post_old[j_post]][j_post];
					inf_t_ind=(jh_learn_syn.t_ind+jh_learn_syn.post_V_hist.size()-jh_learn_syn.inf_steps)%(jh_learn_syn.post_V_hist.size()); 
					expdecay=exp(-(jh_learn_syn.inf_steps-1-age)/jh_learn_syn.tau);
					if(syn_type==0){
						post_h=-(jh_learn_syn.post_V_hist[inf_t_ind][j_post]-V_ex)*expdecay;
					}
					else
					{
						post_h=(jh_learn_syn.post_V_hist[inf_t_ind][j_post]-V_in)*expdecay;
					}
					if(post_h>0){
						jh_learn_syn.rhat[i_pre]+=post_h*K[i_pre][syn_ind];
					}
				}
			}
			jh_learn_syn.rhat[i_pre]*=jh_learn_syn.inf_scale*1.0/jh_learn_syn.C;	
		}
	}
	return jh_learn_syn.rhat;
}

void ChemSyn::old_pre_spikes(){
	if(jh_learn_syn.on){
		int i_pre;
		jh_learn_syn.old_pre.clear();
		for(i_pre=0;i_pre<N_pre;i_pre++){
			for(int i=jh_learn_syn.ind_pre_old[i_pre]; i!=(jh_learn_syn.ind_pre_new[i_pre]+1)%jh_learn_syn.pre_hist_len;i=(i+1)%jh_learn_syn.pre_hist_len){
				//increment time of history pre spikes
				jh_learn_syn.pre_t_hist[i][i_pre]++;
			}

			if(	(((jh_learn_syn.ind_pre_new[i_pre]+1)%jh_learn_syn.pre_hist_len!=jh_learn_syn.ind_pre_old[i_pre]) &&
				(jh_learn_syn.pre_t_hist[jh_learn_syn.ind_pre_old[i_pre]][i_pre]>=jh_learn_syn.inf_steps))	){
				// the oldest pre spike is too old 
				if(jh_learn_syn.pre_noise_hist[jh_learn_syn.ind_pre_old[i_pre]][i_pre]==0){
					//spike was not dropped out
					jh_learn_syn.old_pre.push_back(i_pre);
				}
				jh_learn_syn.ind_pre_old[i_pre]=(jh_learn_syn.ind_pre_old[i_pre]+1)%jh_learn_syn.pre_hist_len;
				
			}
		}
	}
}

void ChemSyn::get_rhat_spiking(){
	if(jh_learn_syn.on){
		int i_pre,j_post,inf_t_ind,age;
		unsigned int syn_ind,i;
		double expdecay;
		double post_h=0;
		for(i=0;i<jh_learn_syn.old_pre.size();i++){
			i_pre=jh_learn_syn.old_pre[i];
			jh_learn_syn.rhat[i_pre]=0;
			for (syn_ind = 0; syn_ind < C[i_pre].size(); ++syn_ind){
				j_post = C[i_pre][syn_ind];
				if(jh_learn_syn.ind_post_old[j_post]!=(jh_learn_syn.ind_post_new[j_post]+1)%jh_learn_syn.post_hist_len){
					//there is a j_post spike in history
					// if(jh_learn_syn.post_noise_hist[jh_learn_syn.ind_post_old[j_post]][j_post]==0){
						//spike was not dropped out
						age=jh_learn_syn.post_t_hist[jh_learn_syn.ind_post_old[j_post]][j_post];
						inf_t_ind=(jh_learn_syn.t_ind+jh_learn_syn.post_V_hist.size()-jh_learn_syn.inf_steps)%(jh_learn_syn.post_V_hist.size()); 
						expdecay=exp(-(jh_learn_syn.inf_steps-1-age)/jh_learn_syn.tau);
						if(syn_type==0){
							post_h=-(jh_learn_syn.post_V_hist[inf_t_ind][j_post]-V_ex)*expdecay;
						}
						else
						{
							post_h=(jh_learn_syn.post_V_hist[inf_t_ind][j_post]-V_in)*expdecay;
						}
						if(post_h>0){
							jh_learn_syn.rhat[i_pre]+=post_h*K[i_pre][syn_ind];
						}
					// }
				}
			}
			jh_learn_syn.rhat[i_pre]*=jh_learn_syn.inf_scale*(1/jh_learn_syn.C);
		}
	}
}

void ChemSyn::wchange_Hebbian_outgoing(){
	if(jh_learn_syn.on){
		int i_pre,j_post,inf_t_ind,age;
		unsigned int syn_ind,i;
		double expdecay;
		double post_h=0;
		
		for(i=0;i<jh_learn_syn.old_pre.size();i++){
			i_pre=jh_learn_syn.old_pre[i];	

			// If rhat is less than noise then make it noise level
			if(jh_learn_syn.rhat[i_pre]<jh_learn_syn.noise){
				jh_learn_syn.rhat[i_pre]=jh_learn_syn.noise;
			}
			for (syn_ind = 0; syn_ind < C[i_pre].size(); ++syn_ind){
				j_post = C[i_pre][syn_ind];
				if(jh_learn_syn.ind_post_old[j_post]!=(jh_learn_syn.ind_post_new[j_post]+1)%jh_learn_syn.post_hist_len){
					//there is a j_post spike in history
					if(jh_learn_syn.post_noise_hist[jh_learn_syn.ind_post_old[j_post]][j_post]==0){
						//spike was not dropped out
						age=jh_learn_syn.post_t_hist[jh_learn_syn.ind_post_old[j_post]][j_post];
						inf_t_ind=(jh_learn_syn.t_ind+jh_learn_syn.post_V_hist.size()-jh_learn_syn.inf_steps)%(jh_learn_syn.post_V_hist.size()); 
						expdecay=exp(-(jh_learn_syn.inf_steps-1-age)/jh_learn_syn.tau); 
						if(syn_type==0){
							post_h=-(jh_learn_syn.post_V_hist[inf_t_ind][j_post]-V_ex)*expdecay/((1.0-jh_learn_syn.noise_pre)*(1.0-jh_learn_syn.noise_post));
						}
						else
						{
							post_h=(jh_learn_syn.post_V_hist[inf_t_ind][j_post]-V_in)*expdecay/((1.0-jh_learn_syn.noise_pre)*(1.0-jh_learn_syn.noise_post));
						}
						if(post_h>0){
							K[i_pre][syn_ind]+=jh_learn_syn.learn_rate*(1.0/jh_learn_syn.C)*post_h/jh_learn_syn.rhat[i_pre]*(K[i_pre][syn_ind]);
							if(K[i_pre][syn_ind]<0){
								K[i_pre][syn_ind]=0;
							}
						}
					}
				}
				if(jh_learn_syn.learn_rate_all>0){
					//update weights to non spiking neurons as well
					if(jh_learn_syn.post_noise_hist[jh_learn_syn.ind_post_old[j_post]][j_post]==0){
						//spike was not dropped out
						K[i_pre][syn_ind]+=jh_learn_syn.learn_rate_all/jh_learn_syn.rhat[i_pre];

						if(K[i_pre][syn_ind]<0){
							K[i_pre][syn_ind]=0;
						}
					}				
				}	
			}
		}
	}
}


void ChemSyn::recv_pop_data(vector<NeuroPop*> &NeuronPopArray){
// get current spikes from pre-pop
	if (pop_ind_pre >= 0){
		spikes_pre = NeuronPopArray[pop_ind_pre]->get_spikes_current(); // This might be problematic!!!
		spikes_post = NeuronPopArray[pop_ind_post]->get_spikes_current();
		
		if(jh_learn_syn.on){

			if(jh_learn_syn.direction==0){
				jh_learn_syn.spikes_pre_noise = NeuronPopArray[pop_ind_pre]->get_spikes_noise(); 
			}
			else if(jh_learn_syn.direction==1){
				//TODO 
			}
			if(jh_learn_syn.direction==0){
				jh_learn_syn.spikes_post_noise = NeuronPopArray[pop_ind_post]->get_spikes_noise();
			}
			else if(jh_learn_syn.direction==1){
				//TODO 
			}
		}
	}
	// get current V from post-pop
	V_post = NeuronPopArray[pop_ind_post]->get_V(); // This  might be  problematic!!!
}


void ChemSyn::send_pop_data(vector<NeuroPop*> &NeuronPopArray){
	
	NeuronPopArray[pop_ind_post]->recv_I(I, pop_ind_pre, syn_type);

}

void ChemSyn::record_stats(int step_current){
	if (stats.record){
		double mean_tmp_I, var_tmp_I;
		Welford_online(I, mean_tmp_I, var_tmp_I);
		// record   
		stats.I_mean.push_back(mean_tmp_I);
		stats.I_std.push_back( sqrt(var_tmp_I) );
		
		if (synapse_model == 0){
			int k = step_current;
			bool is_end = step_current == (step_tot-1);
			Welford_online(gsm_0.s, stats.s_time_mean, stats.s_time_var, k, is_end);
			Welford_online(I, stats.I_time_mean, stats.I_time_var, k, is_end);
		}
	}
	
	if (stats.record_cov){
		if (synapse_model == 0){
			if (step_current >= stats.time_start_cov && step_current <= stats.time_end_cov){
				int k = step_current-stats.time_start_cov;
				bool is_end = step_current == (stats.time_end_cov-1);
				Welford_online(gsm_0.s, stats.s_time_mean_dumb, stats.s_time_cov, k, is_end);
			}
		}
	}
}

void ChemSyn::import_restart(H5File& file, int syn_ind){

	string str;

	string syn_str = "/syns/syn" + to_string(syn_ind)+"/";

	dt=read_scalar_HDF5<double>(file, syn_str+"dt");
	step_tot=read_scalar_HDF5<int>(file,syn_str+"step_tot");
	pop_ind_pre=read_scalar_HDF5<int>(file,syn_str+"pop_ind_pre");
	pop_ind_post=read_scalar_HDF5<int>(file,syn_str+"pop_ind_post");
	N_pre=read_scalar_HDF5<int>(file, syn_str+"N_pre");
	N_post=read_scalar_HDF5<int>(file, syn_str+"N_post");
	syn_type=read_scalar_HDF5<int>(file, syn_str+"syn_type");
	V_ex=read_scalar_HDF5<double>(file, syn_str+"V_ex");
	V_in=read_scalar_HDF5<double>(file, syn_str+"V_in");
	max_delay_steps=read_scalar_HDF5<int>(file, syn_str+"max_delay_steps");
	read_vector_HDF5(file, syn_str+"V_post",V_post);
	read_vector_HDF5(file, syn_str+"spikes_pre",spikes_pre);
	read_vector_HDF5(file,syn_str+"spikes_post",spikes_post);
	read_vector_HDF5(file, syn_str+"I",I);

		

	str = syn_str+"/Stats/";
	if(group_exist_HDF5(file,str)){
		stats.record=read_scalar_HDF5<bool>(file, str+"record");
		read_vector_HDF5(file, str+"I_mean",stats.I_mean);
		read_vector_HDF5(file, str+"I_std",stats.I_std);
		start_stats_record();
	}

	str =syn_str+"/Sample/";
	if(group_exist_HDF5(file,str)){
		read_vector_HDF5(file, str+"neurons", sample.neurons);
		read_vector_HDF5(file,str+ "time_points", sample.time_points);
		add_sampling(sample.neurons, sample.time_points);
		// read_matrix_HDF5(file, str+ "data",sample.data);
	}

	Dt_trans_AMPA=read_scalar_HDF5<double>(file, syn_str+"Dt_trans_AMPA");
	Dt_trans_GABA=read_scalar_HDF5<double>(file, syn_str+"Dt_trans_GABA");
	Dt_trans_NMDA=read_scalar_HDF5<double>(file, syn_str+"Dt_trans_NMDA");
	tau_decay_AMPA=read_scalar_HDF5<double>(file, syn_str+"tau_decay_AMPA");
	tau_decay_GABA=read_scalar_HDF5<double>(file, syn_str+"tau_decay_GABA");
	tau_decay_NMDA=read_scalar_HDF5<double>(file,syn_str+ "tau_decay_NMDA");
	tau_rise=read_scalar_HDF5<double>(file, syn_str+"tau_rise");
	tau_decay=read_scalar_HDF5<double>(file,syn_str+ "tau_decay");
	steps_trans=read_scalar_HDF5<double>(file,syn_str+ "steps_trans");
	read_vector_HDF5(file, syn_str+"K_trans",K_trans);
	exp_step_decay=read_scalar_HDF5<double>(file, syn_str+"exp_step_decay");
	exp_step_rise=read_scalar_HDF5<double>(file,syn_str+ "exp_step_rise");
	miuMg_NMDA=read_scalar_HDF5<double>(file, syn_str+"miuMg_NMDA");
	gamma_NMDA=read_scalar_HDF5<double>(file, syn_str+"gamma_NMDA");
	B_V_min=read_scalar_HDF5<double>(file, syn_str+"B_V_min");
	B_V_max=read_scalar_HDF5<double>(file, syn_str+"B_V_max");
	B_dV=read_scalar_HDF5<double>(file,syn_str+ "B_dV");
	read_vector_HDF5(file, syn_str+"B",B);

	str =syn_str+"/Inh_STDP/";
	if(group_exist_HDF5(file,str)){
		inh_STDP.on=read_scalar_HDF5<bool>(file, str+ "on");
		read_vector_HDF5(file, str+ "x_trace_pre",inh_STDP.x_trace_pre);
		read_vector_HDF5(file,  str+"x_trace_post",inh_STDP.x_trace_post);
		inh_STDP.tau=read_scalar_HDF5<double>(file,  str+"tau");
		inh_STDP.exp_step=read_scalar_HDF5<double>(file, str+ "exp_step");
		inh_STDP.eta=read_scalar_HDF5<double>(file, str+ "eta");
		inh_STDP.rho_0=read_scalar_HDF5<double>(file, str+ "rho_0");
		inh_STDP.alpha=read_scalar_HDF5<double>(file,  str+"alpha");
		inh_STDP.on_step=read_scalar_HDF5<int>(file,  str+"on_step");
		read_matrix_HDF5(file,  str+"j_2_i",inh_STDP.j_2_i);
		read_matrix_HDF5(file, str+ "j_2_syn_ind",inh_STDP.j_2_syn_ind);
	}

	str =syn_str+"/Std/";
	if(group_exist_HDF5(file,str)){
		STD.on=read_scalar_HDF5<bool>(file, str+ "on");
		STD.p_ves=read_scalar_HDF5<double>(file,  str+"p_ves");
		STD.tau_ves=read_scalar_HDF5<double>(file,  str+"tau_ves");
		STD.exp_ves=read_scalar_HDF5<double>(file,  str+"exp_ves");
		STD.on_step=read_scalar_HDF5<bool>(file,  str+"on_step");
		read_vector_HDF5(file, str+ "f_ves",STD.f_ves);
	}

	synapse_model=read_scalar_HDF5<int>(file,syn_str+ "synapse_model");
	read_vector_HDF5(file, syn_str+"gs_sum",gs_sum);

	str =syn_str+"/Gsm_0/";		
	if(group_exist_HDF5(file,str)){
		gsm_0.buffer_steps=read_scalar_HDF5<int>(file,  str+"buffer_steps");
		read_vector_HDF5(file,  str+"s",gsm_0.s);
		read_vector_HDF5(file,  str+"trans_left",gsm_0.trans_left);
		read_matrix_HDF5(file,  str+"d_gs_sum_buffer",gsm_0.d_gs_sum_buffer);
	}

	str =syn_str+"/Gsm_1/";
	if(group_exist_HDF5(file,str)){
		gsm_1.buffer_steps=read_scalar_HDF5<int>(file,  str+"buffer_steps");
		read_vector_HDF5(file, str+ "gs_rise_sum",gsm_1.gs_rise_sum);
		read_vector_HDF5(file,  str+"gs_decay_sum",gsm_1.gs_decay_sum);
		read_matrix_HDF5(file,  str+"d_gs_rd_sum_buffer",gsm_1.d_gs_rd_sum_buffer);
	}

	read_matrix_HDF5(file,syn_str+ "C",C);
	read_matrix_HDF5(file, syn_str+"D",D);
	read_matrix_HDF5(file,syn_str+"K",K);

	//need to remove trailing zero entries that are padded into the matrix storage but aren't supposed to be there
	// this happens typically when there are no periodic boundaries and neurons near the edge have fewer connections
	for(unsigned int i=0;i<C.size();i++){
		int j=C[i].size()-1;
		while((C[i][j]==-1.0)&(j>=0)){
			C[i].pop_back();
			D[i].pop_back();
			K[i].pop_back();
			j--;
		}
	}

	str =syn_str+"/Ext_noise/";
	if(group_exist_HDF5(file,str)){
		ext_noise.K_ext=read_scalar_HDF5<double>(file, str+ "K_ext");
		ext_noise.Num_ext=read_scalar_HDF5<int>(file,  str+"Num_ext");
		read_vector_HDF5(file,  str+"neurons", ext_noise.neurons);
		read_vector_HDF5(file,  str+"rate_ext_t", ext_noise.rate_ext_t);
	}
	my_seed=read_scalar_HDF5<int>(file, syn_str+"my_seed");
	// +++ TODO BASE_GENERATOR_TYP
	
	// JH Learning
	str = syn_str+"/JH_Learn/";
	if(group_exist_HDF5(file,str)){		
		jh_learn_syn.on=read_scalar_HDF5<bool>(file,str+ "on");
		
		read_matrix_HDF5(file,str+ "post_t_hist",jh_learn_syn.post_t_hist);
		read_matrix_HDF5(file,str+ "post_noise_hist",jh_learn_syn.post_noise_hist);
		read_vector_HDF5(file,str+ "ind_post_new",jh_learn_syn.ind_post_new);
		read_vector_HDF5(file,str+ "ind_post_old",jh_learn_syn.ind_post_old);
		read_matrix_HDF5(file,str+ "pre_t_hist",jh_learn_syn.pre_t_hist);
		read_matrix_HDF5(file,str+ "pre_noise_hist",jh_learn_syn.pre_noise_hist);
		read_vector_HDF5(file,str+ "ind_pre_new",jh_learn_syn.ind_pre_new);
		read_vector_HDF5(file,str+ "ind_pre_old",jh_learn_syn.ind_pre_old);
		read_vector_HDF5(file,str+ "old_pre",jh_learn_syn.old_pre);
		jh_learn_syn.post_hist_len=read_scalar_HDF5<int>(file,str+ "post_hist_len");
		jh_learn_syn.pre_hist_len=read_scalar_HDF5<int>(file,str+ "pre_hist_len");
		jh_learn_syn.ntype_pre=read_scalar_HDF5<int>(file,str+ "ntype_pre");
		jh_learn_syn.ntype_post=read_scalar_HDF5<int>(file,str+ "ntype_post");
		read_vector_HDF5(file,str+ "Vint",jh_learn_syn.Vint);
		read_vector_HDF5(file,str+ "Vint_ctr",jh_learn_syn.Vint_ctr);
		read_vector_HDF5(file,str+ "rhat",jh_learn_syn.rhat);
		read_matrix_HDF5(file,str+ "post_V_hist",jh_learn_syn.post_V_hist);
		read_matrix_HDF5(file,str+ "post_R_hist",jh_learn_syn.post_R_hist);
		jh_learn_syn.t_ind=read_scalar_HDF5<int>(file,str+ "t_ind");
		jh_learn_syn.inf_steps=read_scalar_HDF5<int>(file,str+ "inf_steps");
		jh_learn_syn.inf_scale=read_scalar_HDF5<double>(file,str+ "inf_scale");
		jh_learn_syn.learn_rate=read_scalar_HDF5<double>(file,str+ "learn_rate");
		jh_learn_syn.learn_rate_all=read_scalar_HDF5<double>(file,str+ "learn_rate_all");
		jh_learn_syn.tau=read_scalar_HDF5<double>(file,str+ "tau");
		jh_learn_syn.C=read_scalar_HDF5<double>(file, str+"C");
		jh_learn_syn.noise=read_scalar_HDF5<double>(file,str+ "noise");
		jh_learn_syn.noise_post=read_scalar_HDF5<double>(file,str+ "noise_post");
		jh_learn_syn.noise_pre=read_scalar_HDF5<double>(file,str+ "noise_pre");
		
		read_matrix_HDF5(file,str+"j_2_i",jh_learn_syn.j_2_i);
		//need to remove trailing zero entries that are padded into the matrix storage but aren't supposed to be there
		// this happens typically when there are no periodic boundaries and neurons near the edge have fewer connections
		for(unsigned int i=0;i<jh_learn_syn.j_2_i.size();i++){
			int j=jh_learn_syn.j_2_i[i].size()-1;
			while((jh_learn_syn.j_2_i[i][j]==-1.0)&(j>=0)){
				jh_learn_syn.j_2_i[i].pop_back();
				j--;
			}
		}
		read_matrix_HDF5(file,str+"j_2_syn_ind",jh_learn_syn.j_2_syn_ind);
		//need to remove trailing zero entries that are padded into the matrix storage but aren't supposed to be there
		// this happens typically when there are no periodic boundaries and neurons near the edge have fewer connections
		for(unsigned int i=0;i<jh_learn_syn.j_2_syn_ind.size();i++){
			int j=jh_learn_syn.j_2_syn_ind[i].size()-1;
			while((jh_learn_syn.j_2_syn_ind[i][j]==-1.0)&(j>=0)){
				jh_learn_syn.j_2_syn_ind[i].pop_back();
				j--;
			}
		}

		read_vector_HDF5(file,str+ "spikes_pre_noise",jh_learn_syn.spikes_pre_noise);
		read_vector_HDF5(file,str+ "spikes_post_noise",jh_learn_syn.spikes_post_noise);

		jh_learn_syn.V_rt=read_scalar_HDF5<double>(file,str+ "V_rt");
		jh_learn_syn.V_th=read_scalar_HDF5<double>(file,str+ "V_th");
		jh_learn_syn.direction=read_scalar_HDF5<double>(file,str+ "direction");

		if(jh_learn_syn.direction==1){
			//TODO
		}	
	}
}

void ChemSyn::export_restart(Group& group, int syn_ind){
	string syn_str = "/syns/syn" + to_string(syn_ind);
	Group group_syn = group.createGroup(syn_str);

	write_scalar_HDF5(group_syn,dt, "dt");
	write_scalar_HDF5(group_syn,step_tot, "step_tot");
	write_scalar_HDF5(group_syn,pop_ind_pre,"pop_ind_pre");
	write_scalar_HDF5(group_syn,pop_ind_post,"pop_ind_post");
	write_scalar_HDF5(group_syn,N_pre, "N_pre");
	write_scalar_HDF5(group_syn,N_post, "N_post");
	write_scalar_HDF5(group_syn,syn_type, "syn_type");
	write_scalar_HDF5(group_syn,V_ex, "V_ex");
	write_scalar_HDF5(group_syn,V_in, "V_in");
	write_scalar_HDF5(group_syn,max_delay_steps, "max_delay_steps");
	write_vector_HDF5(group_syn,V_post, "V_post");
	write_vector_HDF5(group_syn,spikes_pre, "spikes_pre");
	write_vector_HDF5(group_syn,spikes_post, "spikes_post");
	write_vector_HDF5(group_syn,I, "I");

		
	if(stats.record){
		string str = syn_str+"/Stats";
		Group group_stats = group_syn.createGroup(str);
		write_scalar_HDF5(group_stats,stats.record, "record");
		// write_vector_HDF5(group_stats,stats.I_mean, "I_mean");
		// write_vector_HDF5(group_stats,stats.I_std, "I_std");
	}

	if(!sample.time_points.empty()){
		string str =syn_str+"/Sample";
		Group group_sample = group_syn.createGroup(str);
		write_vector_HDF5(group_sample, sample.neurons, "neurons");
		write_vector_HDF5(group_sample, sample.time_points, "time_points");
		// write_matrix_HDF5(group_sample, sample.data, "data");
	}

	write_scalar_HDF5(group_syn,Dt_trans_AMPA, "Dt_trans_AMPA");
	write_scalar_HDF5(group_syn,Dt_trans_GABA, "Dt_trans_GABA");
	write_scalar_HDF5(group_syn,Dt_trans_NMDA, "Dt_trans_NMDA");
	write_scalar_HDF5(group_syn,tau_decay_AMPA, "tau_decay_AMPA");
	write_scalar_HDF5(group_syn,tau_decay_GABA, "tau_decay_GABA");
	write_scalar_HDF5(group_syn,tau_decay_NMDA, "tau_decay_NMDA");
	write_scalar_HDF5(group_syn,tau_rise, "tau_rise");
	write_scalar_HDF5(group_syn,tau_decay, "tau_decay");
	write_scalar_HDF5(group_syn,steps_trans, "steps_trans");
	write_vector_HDF5(group_syn,K_trans, "K_trans");
	write_scalar_HDF5(group_syn,exp_step_decay, "exp_step_decay");
	write_scalar_HDF5(group_syn,exp_step_rise, "exp_step_rise");
	write_scalar_HDF5(group_syn,miuMg_NMDA, "miuMg_NMDA");
	write_scalar_HDF5(group_syn,gamma_NMDA, "gamma_NMDA");
	write_scalar_HDF5(group_syn,B_V_min, "B_V_min");
	write_scalar_HDF5(group_syn,B_V_max, "B_V_max");
	write_scalar_HDF5(group_syn,B_dV, "B_dV");
	write_vector_HDF5(group_syn,B, "B");

	if(inh_STDP.on){
		string str =syn_str+"/Inh_STDP";
		Group group_Inh_STDP = group_syn.createGroup(str);
		write_scalar_HDF5(group_Inh_STDP,inh_STDP.on, "on");
		write_vector_HDF5(group_Inh_STDP,inh_STDP.x_trace_pre, "x_trace_pre");
		write_vector_HDF5(group_Inh_STDP,inh_STDP.x_trace_post, "x_trace_post");
		write_scalar_HDF5(group_Inh_STDP,inh_STDP.tau, "tau");
		write_scalar_HDF5(group_Inh_STDP,inh_STDP.exp_step, "exp_step");
		write_scalar_HDF5(group_Inh_STDP,inh_STDP.eta, "eta");
		write_scalar_HDF5(group_Inh_STDP,inh_STDP.rho_0, "rho_0");
		write_scalar_HDF5(group_Inh_STDP,inh_STDP.alpha, "alpha");
		write_scalar_HDF5(group_Inh_STDP,inh_STDP.on_step, "on_step");
		write_matrix_HDF5(group_Inh_STDP,inh_STDP.j_2_i, "j_2_i");
		write_matrix_HDF5(group_Inh_STDP,inh_STDP.j_2_syn_ind, "j_2_syn_ind");
	}

	if(STD.on){
		string str =syn_str+"/Std";
		Group group_STD = group_syn.createGroup(str);
		write_scalar_HDF5(group_STD,STD.on, "on");
		write_scalar_HDF5(group_STD,STD.p_ves, "p_ves");
		write_scalar_HDF5(group_STD,STD.tau_ves, "tau_ves");
		write_scalar_HDF5(group_STD,STD.exp_ves, "exp_ves");
		write_scalar_HDF5(group_STD,STD.on_step, "on_step");
		write_vector_HDF5(group_STD,STD.f_ves, "f_ves");
	}

	write_scalar_HDF5(group_syn,synapse_model, "synapse_model");
	write_vector_HDF5(group_syn,gs_sum, "gs_sum");
		
	if(synapse_model==0){
		string str =syn_str+"/Gsm_0";
		Group group_Gsm_0 = group_syn.createGroup(str);
		write_scalar_HDF5(group_Gsm_0,gsm_0.buffer_steps, "buffer_steps");
		write_vector_HDF5(group_Gsm_0,gsm_0.s, "s");
		write_vector_HDF5(group_Gsm_0,gsm_0.trans_left, "trans_left");
		write_matrix_HDF5(group_Gsm_0,gsm_0.d_gs_sum_buffer, "d_gs_sum_buffer");
	}
	
	if(synapse_model==1){
		string str =syn_str+"/Gsm_1";
		Group group_Gsm_1 = group_syn.createGroup(str);
		write_scalar_HDF5(group_Gsm_1,gsm_1.buffer_steps, "buffer_steps");
		write_vector_HDF5(group_Gsm_1,gsm_1.gs_rise_sum, "gs_rise_sum");
		write_vector_HDF5(group_Gsm_1,gsm_1.gs_decay_sum, "gs_decay_sum");
		write_matrix_HDF5(group_Gsm_1,gsm_1.d_gs_rd_sum_buffer, "d_gs_rd_sum_buffer");
	}
	
	if (pop_ind_pre >= 0){
		// every neuron doesn't need to have the same number of connections
		// so when stored as a matrix we need to pad out missing connections
		// we choose to pad by -1, as we assume entries are always positive,
		// so they can then be identified as padding later
		
		//find max number of entries
		unsigned int max=0;
		for(unsigned int i=0;i<C.size();i++){
			if(C[i].size()>max){
				max=C[i].size();
			}
		}
		for(unsigned int i=0;i<C.size();i++){
			C[i].resize(max,-1.0);
			D[i].resize(max,-1.0);
			K[i].resize(max,-1.0);
		}
		write_matrix_HDF5(group_syn, C, "C");
		write_matrix_HDF5(group_syn, D, "D");
		write_matrix_HDF5(group_syn, K, "K");
	}

	if(!ext_noise.neurons.empty()){
		string str = syn_str+"/Ext_noise";
		Group group_Ext_noise = group_syn.createGroup(str);
		write_scalar_HDF5(group_Ext_noise, ext_noise.K_ext, "K_ext");
		write_scalar_HDF5(group_Ext_noise, ext_noise.Num_ext, "Num_ext");
		write_vector_HDF5(group_Ext_noise, ext_noise.neurons, "neurons");
		
		write_vector_HDF5(group_Ext_noise, ext_noise.rate_ext_t, "rate_ext_t");
	}
	write_scalar_HDF5(group_syn, my_seed, "my_seed");
	// +++ TODO BASE_GENERATOR_TYP

	// JH Learning
	if (jh_learn_syn.on){
		string str = syn_str+"/JH_Learn/";
		Group group_JH_Learn = group_syn.createGroup(str);
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.on, "on");
		
		write_matrix_HDF5(group_JH_Learn,jh_learn_syn.post_t_hist, "post_t_hist");
		write_matrix_HDF5(group_JH_Learn,jh_learn_syn.post_noise_hist, "post_noise_hist");
		write_vector_HDF5(group_JH_Learn,jh_learn_syn.ind_post_new, "ind_post_new");
		write_vector_HDF5(group_JH_Learn,jh_learn_syn.ind_post_old, "ind_post_old");
		write_matrix_HDF5(group_JH_Learn,jh_learn_syn.pre_t_hist, "pre_t_hist");
		write_matrix_HDF5(group_JH_Learn,jh_learn_syn.pre_noise_hist, "pre_noise_hist");
		write_vector_HDF5(group_JH_Learn,jh_learn_syn.ind_pre_new, "ind_pre_new");
		write_vector_HDF5(group_JH_Learn,jh_learn_syn.ind_pre_old, "ind_pre_old");
		write_vector_HDF5(group_JH_Learn,jh_learn_syn.old_pre, "old_pre");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.post_hist_len, "post_hist_len");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.pre_hist_len, "pre_hist_len");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.ntype_pre, "ntype_pre");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.ntype_post, "ntype_post");
		write_vector_HDF5(group_JH_Learn,jh_learn_syn.Vint, "Vint");
		write_vector_HDF5(group_JH_Learn,jh_learn_syn.Vint_ctr, "Vint_ctr");
		write_vector_HDF5(group_JH_Learn,jh_learn_syn.rhat, "rhat");
		write_matrix_HDF5(group_JH_Learn,jh_learn_syn.post_V_hist, "post_V_hist");
		write_matrix_HDF5(group_JH_Learn,jh_learn_syn.post_R_hist, "post_R_hist");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.t_ind, "t_ind");

		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.inf_steps, "inf_steps");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.inf_scale, "inf_scale");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.learn_rate, "learn_rate");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.learn_rate_all, "learn_rate_all");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.tau, "tau");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.C, "C");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.noise, "noise");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.noise_pre, "noise_pre");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.noise_post, "noise_post");
		// every neuron doesn't need to have the same number of connections
		// so when stored as a matrix we need to pad out missing connections
		// we choose to pad by -1, as we assume entries are always positive,
		// so they can then be identified as padding later
			
		//find max number of entries
		unsigned int max=0;
		for(unsigned int i=0;i<jh_learn_syn.j_2_i.size();i++){
			if(jh_learn_syn.j_2_i[i].size()>max){
				max=jh_learn_syn.j_2_i[i].size();
			}
		}
		for(unsigned int i=0;i<jh_learn_syn.j_2_i.size();i++){
			jh_learn_syn.j_2_i[i].resize(max,-1.0);
		}
		write_matrix_HDF5(group_JH_Learn,jh_learn_syn.j_2_i,"j_2_i");
		// every neuron doesn't need to have the same number of connections
		// so when stored as a matrix we need to pad out missing connections
		// we choose to pad by -1, as we assume entries are always positive,
		// so they can then be identified as padding later
			
		//find max number of entries
		max=0;
		for(unsigned int i=0;i<jh_learn_syn.j_2_syn_ind.size();i++){
			if(jh_learn_syn.j_2_syn_ind[i].size()>max){
				max=jh_learn_syn.j_2_syn_ind[i].size();
			}
		}
		for(unsigned int i=0;i<jh_learn_syn.j_2_syn_ind.size();i++){
			jh_learn_syn.j_2_syn_ind[i].resize(max,-1.0);
		}
		write_matrix_HDF5(group_JH_Learn,jh_learn_syn.j_2_syn_ind ,"j_2_syn_ind");

		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.V_th, "V_th");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.V_rt, "V_rt");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.direction, "direction");

		write_vector_HDF5(group_JH_Learn,jh_learn_syn.spikes_post_noise, "spikes_post_noise");
		write_vector_HDF5(group_JH_Learn,jh_learn_syn.spikes_pre_noise, "spikes_pre_noise");
		
		if(jh_learn_syn.direction==1){
			//TODO
		}
	}
}
void ChemSyn::output_results(H5File& file, int syn_ind){
	// new group
	stringstream group_name;
	group_name << "/syn_result_"  << syn_ind;
	Group group_syn = file.createGroup(group_name.str());
	
	write_string_HDF5(group_syn, dump_para(), string("syn_para"));
		
	if (!sample.neurons.empty()){
		write_matrix_HDF5(group_syn, sample.data, string("sample_data"));
	}
	
	if (stats.record){
		write_vector_HDF5(group_syn, stats.I_mean, string("stats_I_mean"));
		write_vector_HDF5(group_syn, stats.I_std, string("stats_I_std"));
		if (synapse_model == 0){
			write_vector_HDF5(group_syn, stats.s_time_mean, string("stats_s_time_mean"));
			write_vector_HDF5(group_syn, stats.s_time_var, string("stats_s_time_var"));
			write_vector_HDF5(group_syn, stats.I_time_mean, string("stats_I_time_mean"));
			write_vector_HDF5(group_syn, stats.I_time_var, string("stats_I_time_var"));
		}
	}
	if (stats.record_cov){
		if (synapse_model == 0){
			write_matrix_HDF5(group_syn, stats.s_time_cov, string("stats_s_time_cov"));
		}
	}
	
}

