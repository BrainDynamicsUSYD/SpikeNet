/// \file
//#include <iostream>
#include <cmath> // always cmath instead of math.h ?
#include <algorithm>
#include <sstream>
//#include <fstream>

#include "ChemSyn.h"

ChemSyn::ChemSyn(double dt_input, int step_tot_input, char delim_input, char indicator_input){
	
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
	stats_record = false;
	STD = false; 
	STD_on_step = -1;
	inh_STDP_on_step = -1;
	inh_STDP = false;
	synapse_model = 0; // default model
	
}



void ChemSyn::init(int syn_type_input, int i_pre, int j_post, int N_pre_input, int N_post_input, vector<int> &C_i, vector<int> &C_j, vector<double> &K_ij, vector<double> &D_ij){

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


void ChemSyn::init(int syn_type_input, int j_post, int N_post_input, double K_ext_input, int Num_ext_input, vector<double> &rate_ext_t_input, int ia_input, int ib_input){

	// Initialise chemical synapses for simulating external neuron population
	syn_type = syn_type_input;
	pop_ind_pre = -1; // -1 for external noisy population
	pop_ind_post = j_post;
	N_pre = 1; // just for initialization
	N_post = N_post_input;
	max_delay_steps = 0; // no delay;
	ia = ia_input; // start of neuron index range in post population
	ib = ib_input; // end of neuron index range in post population


	// Parameters for noise generation
	K_ext = K_ext_input;
	Num_ext = Num_ext_input;
	rate_ext_t = rate_ext_t_input;	

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
	
	
	// Synapse model choice
	// model 0, the default model
	if (synapse_model == 0){ 
		// Initialize pre- and post-synaptic dynamic variables
		buffer_steps = max_delay_steps + steps_trans + 1; // the +1 is vital
		s.assign(N_pre, 0);
		d_gs_sum_buffer.resize(buffer_steps);
		for (int i = 0; i < buffer_steps; ++i){
			d_gs_sum_buffer[i].assign(N_post, 0);
		}
		// transmitter_strength
		K_trans.assign(N_pre, 1.0 / steps_trans); // be careful! 1 / transmitter steps gives zero (int)!!
		trans_left.assign(N_pre, 0);
	}
	else if (synapse_model == 1){
		// model 1
		buffer_steps = max_delay_steps + 1; // the +1 is vital
		gs_rise_sum.assign(N_post, 0);
		gs_decay_sum.assign(N_post, 0);
		d_gs_rd_sum_buffer.resize(buffer_steps);
		for (int i = 0; i < buffer_steps; ++i){
			d_gs_rd_sum_buffer[i].assign(N_post, 0);
		}
		// clear model 0
		s.clear();
		d_gs_sum_buffer.clear(); // clear() clears all of its components recursively
		K_trans.clear();
		trans_left.clear();
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

void ChemSyn::set_synapse_model(int synapse_model_input){
	if (synapse_model_input != 0){
		synapse_model = synapse_model_input;
		init(); // initialize again
	}
}

	

void ChemSyn::update(int step_current){

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



void ChemSyn::update_STD(int step_current){
	// STD modifies K_trans
	
	if (STD_on_step == step_current){STD = true;}
	// short-term depression
	if (STD == true){
		for (unsigned int i = 0; i < spikes_pre.size(); ++i){ 
			K_trans[spikes_pre.at(i)] = 1.0 / steps_trans * f_ves[i];
			f_ves[spikes_pre.at(i)] *= 1.0 - p_ves; // decrease at spikes
		}
		for (int i = 0; i < N_pre; ++i){
			f_ves[i] = 1.0 - exp_ves * (1.0 - f_ves[i]); // decay to 1.0
		}
	}

}


void ChemSyn::add_inh_STDP(int inh_STDP_on_step_input){
	if (syn_type != 1){
		cout << "Warning: initializing inhibitory STDP on non-GABA synapses!" << endl;
	}
	if (synapse_model != 0){
		cout << "Warning: inhibitory STDP is not supported for synapse models other than 0 yet!" << endl;
	}
	
	inh_STDP_on_step = inh_STDP_on_step_input;
	if (inh_STDP_on_step == 0){
		inh_STDP = true;
	}
	
	x_trace_pre.assign(N_pre, 0.0);
	x_trace_post.assign(N_post, 0.0);
	
	tau_STDP = 20; // ms
	exp_step_STDP = exp(-dt / tau_STDP);
	eta_STDP = 0.0001; // learning rate, 0.0001 is the published value but requires 60min of simulation
	rho_0_STDP = 0.010; //kHz
	alpha_STDP = 2.0 * rho_0_STDP * tau_STDP; // depression factor

	// j_2_i and j_2_syn_ind
	j_2_i.resize(N_post);
	j_2_syn_ind.resize(N_post);
	int j_post;
	for (int i_pre = 0; i_pre < N_pre; ++i_pre){ 
		for (unsigned int syn_ind = 0; syn_ind < C[i_pre].size(); ++syn_ind){
			j_post = C[i_pre][syn_ind];
			j_2_i[j_post].push_back( i_pre );
			j_2_syn_ind[j_post].push_back( int(syn_ind) );
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

void ChemSyn::update_gs_sum_model_0(int step_current){
	// See Gu, Yifan, Gong, Pulin, 2016, The dynamics of memory retrieval in hierarchical networks: a modeling study
	// this function updates gs_sum
	if (pop_ind_pre >= 0){
		for (unsigned int i = 0; i < spikes_pre.size(); ++i){ // add spikes (transmitter release)
			trans_left[spikes_pre.at(i)] += steps_trans;
		}
		for (int i_pre = 0; i_pre < N_pre; ++i_pre){
			if (trans_left[i_pre] > 0){
				// conduction delay: put into d_gs_sum_buffer
				int j_post, delay_step, t_ring;
				for (unsigned int syn_ind = 0; syn_ind < C[i_pre].size(); ++syn_ind){ //loop through all the post-synapses
					j_post = C[i_pre][syn_ind]; // index of the post-synaptic neuron
					delay_step = D[i_pre][syn_ind]; // delay in steps for this post-synaptic neuron
					t_ring = int( (step_current + delay_step) % buffer_steps ); // index in the gs_buffer
					d_gs_sum_buffer[t_ring][j_post] += K_trans[i_pre] * (1.0 - s[i_pre]) * K[i_pre][syn_ind];
				}
				trans_left[i_pre] -= 1;
				s[i_pre] += K_trans[i_pre] * (1.0 - s[i_pre]);
			}
		}
		// decay pre-synaptic dynamics
		for (int i = 0; i < N_pre; ++i){ s[i] *= exp_step_decay; };
	}
	else if (pop_ind_pre == -1){ // if external noisy population
		// Contribution of external spikes, assuming square pulse transmitter release
		// Generate current random number generator, note that rate_ext_t is in Hz
		gen.seed(my_seed + step_current);// reseed random engine!!!
		poisson_distribution<int> dist(Num_ext * rate_ext_t[step_current] * (dt / 1000.0));		
		auto ext_spikes = bind(dist, gen);

		// Post-synaptic dynamics
		int t_ring;
		for (int t_trans = 0; t_trans < steps_trans; ++t_trans){
			t_ring = int( (step_current + t_trans) % buffer_steps );
			for (int j = ia; j <= ib; ++j){
				d_gs_sum_buffer[t_ring][j] += K_trans[0] * K_ext * ext_spikes(); 
			}
		}
	}
	// update post-synaptic dynamics
	int t_ring = int( step_current % buffer_steps );
	for (int j_post = 0; j_post < N_post; ++j_post){
		gs_sum[j_post] += d_gs_sum_buffer[t_ring][j_post];
		// should I decay gs_sum here??
	}
	// immediately reset the current buffer to zeros after being used!!
	fill(d_gs_sum_buffer[t_ring].begin(), d_gs_sum_buffer[t_ring].end(), 0.0);
}

void ChemSyn::update_gs_sum_model_1(int step_current){
	// See Keane, A., Gong, P., 2015, Propagating Waves Can Explain Irregular Neural Dynamics
	// this function updates gs_sum
	for (unsigned int ind = 0; ind < spikes_pre.size(); ++ind){ // loop through all the spikes
		int i_pre = spikes_pre.at(ind);
		int j_post, delay_step, t_ring;
		for (unsigned int syn_ind = 0; syn_ind < C[i_pre].size(); ++syn_ind){ //loop through all the post-synapses
			j_post = C[i_pre][syn_ind]; // index of the post-synaptic neuron
			delay_step = D[i_pre][syn_ind]; // delay in steps for this post-synaptic neuron
			t_ring = int( (step_current + delay_step) % buffer_steps ); // index in the gs_buffer
			d_gs_rd_sum_buffer[t_ring][j_post] += K[i_pre][syn_ind];  // the peak value is linear to the initial impulse. 
		}
	}
	int t_ring = int( step_current % buffer_steps );
	for (int j_post = 0; j_post < N_post; ++j_post){ // Check the error of the following numerical scheme!
		gs_rise_sum[j_post] *= exp_step_rise;
		gs_decay_sum[j_post] *= exp_step_decay;
		gs_sum[j_post] = (gs_decay_sum[j_post] - gs_rise_sum[j_post]) / (tau_decay - tau_rise); 
		gs_rise_sum[j_post] += d_gs_rd_sum_buffer[t_ring][j_post];
		gs_decay_sum[j_post] += d_gs_rd_sum_buffer[t_ring][j_post];
	}
	// immediately reset the current buffer to zeros after being used!!
	fill(d_gs_rd_sum_buffer[t_ring].begin(), d_gs_rd_sum_buffer[t_ring].end(), 0.0);
}


void ChemSyn::add_sampling(vector<int> sample_neurons_input, vector<bool> sample_time_points_input){
	sample_neurons = sample_neurons_input;
	sample_time_points = sample_time_points_input;
	
	// initialise
	int sample_time_points_tot = 0;// count non zero elements in sample_time_points
	for (unsigned int i = 0; i < sample_time_points.size(); ++i){
		if (sample_time_points[i]){
			sample_time_points_tot += 1;
		}
	}
	int sample_neurons_tot = sample_neurons.size();// count non zero elements in sample_time_points

	sample.resize(sample_neurons_tot);
	for (int i = 0; i < sample_neurons_tot; ++i){
		sample[i].reserve(sample_time_points_tot); // reserve and push_back so that it won't be affected by adapting step_tot
	}

}

void ChemSyn::add_short_term_depression(int STD_on_step_input){
	if (syn_type != 0){
		cout << "Warning: initializing STD on non-AMPA synapses!" << endl;
	}
	if (synapse_model != 0){
		cout << "Warning: STD is not supported for synapse models other than 0 yet!" << endl;
	}
	
	STD_on_step = STD_on_step_input;
	if (STD_on_step == 0){
		STD = true;
	}
	// short term depression
	p_ves =  0.4; // see X. Wang, 1999, The Journal of Neuroscience
	tau_ves =  700; // ms
	f_ves.assign(N_pre, 1.0);
	exp_ves = exp(-dt / tau_ves);
}




void ChemSyn::update_inh_STDP(int step_current){
	// inh_STDP modifies K
	
	if (inh_STDP_on_step == step_current){inh_STDP = true;}
	if (inh_STDP == true){
		// update x_trace
		for (unsigned int i = 0; i < spikes_pre.size(); ++i){
			x_trace_pre[spikes_pre.at(i)] += 1.0;
		}
		for (unsigned int j = 0; j < spikes_post.size(); ++j){
			x_trace_post[spikes_post.at(j)] += 1.0;
		}
		// update K
		int i_pre, j_post;
		for (unsigned int ind_spike = 0; ind_spike < spikes_pre.size(); ++ind_spike){
			i_pre = spikes_pre.at(ind_spike);
			for (unsigned int syn_ind = 0; syn_ind < C[i_pre].size(); ++syn_ind){
				j_post = C[i_pre][syn_ind];
				K[i_pre][syn_ind] += eta_STDP * ( x_trace_post[j_post] - alpha_STDP );
			}
		}
		int syn_ind;
		for (unsigned int ind_spike = 0; ind_spike < spikes_post.size(); ++ind_spike){
			j_post = spikes_post.at(ind_spike);
			for (unsigned int ind = 0; ind < j_2_i[j_post].size(); ++ind){
				// j_2_i and j_2_syn_ind together serve as the "inverse function" of j_post = C[i_pre][syn_ind]
				i_pre = j_2_i[j_post][ind];
				syn_ind = j_2_syn_ind[j_post][ind];
				K[i_pre][syn_ind] += eta_STDP * x_trace_pre[i_pre];
			}
		}
		// for testing
		//if (tmp_data.size() == 0){
		//	tmp_data.resize(2);
		//}
		//tmp_data[0].push_back(K[0][0]);
		//tmp_data[1].push_back(K[100][0]);
		for (int i = 0; i < N_pre; ++i){ x_trace_pre[i] *= exp_step_STDP; }
		for (int j = 0; j < N_post; ++j){ x_trace_post[j] *= exp_step_STDP; }
	}
}


void ChemSyn::sample_data(int step_current){
	if (!sample_neurons.empty()){
		if (sample_time_points[step_current]){ // push_back is amazing
			for (unsigned int i = 0; i < sample_neurons.size(); ++i){ // performance issue when sampling many neurons?
				int ind_temp = sample_neurons[i];
				sample[i].push_back( I[ind_temp] );
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
	stats_record = true;
	I_mean.reserve(step_tot);
	I_std.reserve(step_tot);
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
	if (!sample_neurons.empty()){
		output_file << indicator << " SYND002" << endl;
		output_file << pop_ind_pre << delim << pop_ind_post << delim << syn_type << delim << sample_neurons.size() << delim << endl;
		write2file(output_file, sample); // 2D matrix
	}

	
	// SYND003 # currents mean and std
	if (stats_record){
		output_file << indicator << " SYND003" << endl;
		output_file << pop_ind_pre << delim << pop_ind_post << delim << syn_type << delim << endl;
		write2file(output_file, I_mean);
		write2file(output_file, I_std);
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


// Use function templates when you want to perform the same action on types that can be different.
// Use function overloading when you want to apply different operations depending on the type.
// In this case, just save yourself the trouble and use overloading.
void ChemSyn::write2file(ofstream& output_file, vector< vector<int> >& v){
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



void ChemSyn::write2file(ofstream& output_file, vector< vector<double> >& v){
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


void ChemSyn::write2file(ofstream& output_file, vector<int>& v){
	if (!v.empty()){
		//for (int f : v){ output_file << f << delim; } // range-based "for" in C++11
		for (unsigned int i = 0; i < v.size(); ++i){
			output_file << v[i] << delim;
		}
		output_file << endl;
	}
	else {output_file << " " << endl;}
}

void ChemSyn::write2file(ofstream& output_file, vector<double>& v){
	if (!v.empty()){
		//for (int f : v){ output_file << f << delim; } // range-based "for" in C++11
		for (unsigned int i = 0; i < v.size(); ++i){
			output_file << v[i] << delim;
		}
		output_file << endl;
	}
	else {output_file << " " << endl;}
}



void ChemSyn::record_stats(){
	if (stats_record){
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
		I_mean.push_back(mean_tmp);
		I_std.push_back(std_tmp);
	}
}