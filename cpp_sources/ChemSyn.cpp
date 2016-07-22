/// \file
//#include <iostream>
#include <cmath> // always cmath instead of math.h ?
#include <algorithm>
#include <sstream>
//#include <fstream>

#include "ChemSyn.h"

ChemSyn::ChemSyn(const double dt_input, const int step_tot_input, const char delim_input, const char indicator_input){
	
	dt = dt_input;
	step_tot = step_tot_input;
	delim = delim_input;
	indicator = indicator_input;
	
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
	STD.on = false; 
	STD.on_step = -1;
	inh_STDP.on_step = -1;
	inh_STDP.on = false;
	synapse_model = 0; // default model
	
}



void ChemSyn::init(const int syn_type_input, const int i_pre, const int j_post, const int N_pre_input, const int N_post_input, const vector<int> &C_i, const vector<int> &C_j, const vector<double> &K_ij, const vector<double> &D_ij){

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
	record_stats();
	
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

	// decay gs_sum
	// numerical error of this integration scheme should be less 1.5%
	for (int j = 0; j < N_post; ++j){ gs_sum[j] *= exp_step_decay; };
	
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
		gen.seed(my_seed + step_current);// reseed random engine!!!
		poisson_distribution<int> dist(ext_noise.Num_ext * ext_noise.rate_ext_t[step_current] * (dt / 1000.0));		
		auto ext_spikes = bind(dist, gen);

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
	// update post-synaptic dynamics
	int t_ring = int( step_current % gsm_0.buffer_steps );
	for (int j_post = 0; j_post < N_post; ++j_post){
		gs_sum[j_post] += gsm_0.d_gs_sum_buffer[t_ring][j_post];
		// should I decay gs_sum here??
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
}

void ChemSyn::output_results(ofstream& output_file){
	// SYND001 # synapse parameters
	// count number of variables
	stringstream dump_count;
	string para_str = dump_para();
	dump_count << para_str;
	string str_temp;
	int var_number = 0;
	while(getline(dump_count,str_temp)){++var_number;} // count number of variables
	output_file << indicator << " SYND001" << endl;
	output_file << var_number << delim << endl;
	output_file << para_str;
	
	
	// SYND002 # sampled synapse data
	if (!sample.neurons.empty()){
		output_file << indicator << " SYND002" << endl;
		output_file << pop_ind_pre << delim << pop_ind_post << delim << syn_type << delim << sample.neurons.size() << delim << endl;
		write2file(output_file, sample.data); // 2D matrix
	}

	
	// SYND003 # currents mean and std
	if (stats.record){
		output_file << indicator << " SYND003" << endl;
		output_file << pop_ind_pre << delim << pop_ind_post << delim << syn_type << delim << endl;
		write2file(output_file, stats.I_mean);
		write2file(output_file, stats.I_std);
	}
	
	// tmp data
	if (tmp_data.size() != 0){
		output_file << indicator << " SYND004" << endl;
		output_file << pop_ind_pre << delim << pop_ind_post << delim << syn_type << delim << tmp_data.size() << delim << endl;
		write2file(output_file, tmp_data);	
	}
}



void ChemSyn::recv_pop_data(vector<NeuroPop*> &NeuronPopArray){
	// get current spikes from pre-pop
	if (pop_ind_pre >= 0){
		spikes_pre = NeuronPopArray[pop_ind_pre]->get_spikes_current(); // This might be problematic!!!
		spikes_post = NeuronPopArray[pop_ind_post]->get_spikes_current();
	}
	// get current V from post-pop
	V_post = NeuronPopArray[pop_ind_post]->get_V(); // This  might be  problematic!!!
}


void ChemSyn::send_pop_data(vector<NeuroPop*> &NeuronPopArray){
	
	NeuronPopArray[pop_ind_post]->recv_I(I, pop_ind_pre, syn_type);

}

void ChemSyn::record_stats(){
	if (stats.record){
		// get mean
		double sum_mean = 0.0;
		for (unsigned int i = 0; i < I.size(); ++i){
			sum_mean += I[i];
		}
		double mean_tmp = sum_mean / double(I.size());
	
		// get std
		double sum_std = 0.0;
		for (unsigned int i = 0; i < I.size(); ++i){
			sum_std += (I[i]-mean_tmp)*(I[i]-mean_tmp);
		}
		double std_tmp = sqrt( sum_std / double(I.size()));
	
		// record   
		stats.I_mean.push_back(mean_tmp);
		stats.I_std.push_back(std_tmp);
	}
}


// Use function templates when you want to perform the same action on types that can be different.
// Use function overloading when you want to apply different operations depending on the type.
// In this case, just save yourself the trouble and use overloading.
void ChemSyn::write2file(ofstream& output_file, const vector< vector<int> >& v){
	if (!v.empty()){
		for (unsigned int i = 0; i < v.size(); ++i){
			//for (double f : v[i]){ output_file << f << delim; } // range-based "for" in C++11
			for (unsigned int j = 0; j < v[i].size(); ++j){
				output_file << v[i][j] << delim;
			}
			output_file << endl;
		}
	}
	else {output_file << " " << endl;}
}



void ChemSyn::write2file(ofstream& output_file, const vector< vector<double> >& v){
	if (!v.empty()){
		for (unsigned int i = 0; i < v.size(); ++i){
			//for (double f : v[i]){ output_file << f << delim; } // range-based "for" in C++11
			for (unsigned int j = 0; j < v[i].size(); ++j){
				output_file << v[i][j] << delim;
			}
			output_file << endl;
		}
	}
	else {output_file << " " << endl;}
}


void ChemSyn::write2file(ofstream& output_file, const vector<int>& v){
	if (!v.empty()){
		//for (int f : v){ output_file << f << delim; } // range-based "for" in C++11
		for (unsigned int i = 0; i < v.size(); ++i){
			output_file << v[i] << delim;
		}
		output_file << endl;
	}
	else {output_file << " " << endl;}
}

void ChemSyn::write2file(ofstream& output_file, const vector<double>& v){
	if (!v.empty()){
		//for (int f : v){ output_file << f << delim; } // range-based "for" in C++11
		for (unsigned int i = 0; i < v.size(); ++i){
			output_file << v[i] << delim;
		}
		output_file << endl;
	}
	else {output_file << " " << endl;}
}


#ifdef HDF5
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
	}
}


void ChemSyn::write_scalar_HDF5(Group & group, int s, const string & v_name){
	vector<int> v_tmp;
	v_tmp.push_back(s);
	write_vector_HDF5(group, v_tmp, v_name);
}

void ChemSyn::write_scalar_HDF5(Group & group, double s, const string & v_name){
	vector<double> v_tmp;
	v_tmp.push_back(s);
	write_vector_HDF5(group, v_tmp, v_name);
}



void ChemSyn::write_vector_HDF5(Group & group, const vector<int> & v, const string & v_name){
	hsize_t dims[1]; 
	dims[0] = v.size();
	DataSpace fspace(1, dims); 
	DataSet v_dataset = group.createDataSet(v_name, PredType::NATIVE_INT32, fspace);
	v_dataset.write( v.data(), PredType::NATIVE_INT32, fspace, fspace );	
}

void ChemSyn::write_vector_HDF5(Group & group, const vector<double> & v, const string &   v_name){
	hsize_t dims[1]; 
	dims[0] = v.size();
	DataSpace fspace(1, dims); 
	DataSet v_dataset = group.createDataSet(v_name, PredType::NATIVE_DOUBLE, fspace);
	v_dataset.write( v.data(), PredType::NATIVE_INT32, fspace, fspace );	
}

void ChemSyn::write_string_HDF5(Group & group, const string & s, const string &  s_name){
   // HDF5 only understands vector of char* :-(
   vector<const char*> arr_c_str;
   arr_c_str.push_back(s.c_str());

   hsize_t str_dimsf[1] {arr_c_str.size()};
   DataSpace dataspace(1, str_dimsf);

   // Variable length string
   StrType datatype(PredType::C_S1, H5T_VARIABLE); 
   DataSet str_dataset = group.createDataSet(s_name, datatype, dataspace);

   str_dataset.write(arr_c_str.data(), datatype);
}

void ChemSyn::write_matrix_HDF5(Group & group, const vector< vector<double> > & m, const string &  m_name){
	hsize_t dims[2]; 
	dims[0] = m.size();
	dims[1] = m[0].size();
	DataSpace fspace(2, dims); 
	DataSet m_dataset = group.createDataSet(m_name, PredType::NATIVE_DOUBLE, fspace);
	for (int i = 0; i < int(m.size()); i++){
		append_vector_to_matrix_HDF5(m_dataset, m[i],  i);
	}
}


void ChemSyn::append_vector_to_matrix_HDF5(DataSet & dataset_tmp, const vector<double> & v, const int colNum){
   	hsize_t offset[2];
    offset[1] = 0;
    offset[0] = colNum;
	
    hsize_t fdims[2];            // new data dimensions 
	fdims[0] = 1;
	fdims[1] = v.size();
	DataSpace mspace( 2, fdims );

	DataSpace fspace = dataset_tmp.getSpace();
    fspace.selectHyperslab( H5S_SELECT_SET, fdims, offset );
    dataset_tmp.write(  v.data(), PredType::NATIVE_DOUBLE, mspace, fspace );	
}

#endif
