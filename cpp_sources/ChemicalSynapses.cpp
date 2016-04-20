#include <iostream>
#include <cmath> // always cmath instead of math.h ?
#include <algorithm>
#include <sstream>
#include <fstream>

#include "ChemicalSynapses.h"

ChemicalSynapses::ChemicalSynapses(double dt_input, int step_tot_input){
	
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

	// 
	stats_record = false;
}



void ChemicalSynapses::init(int synapses_type_input, int i_pre, int j_post, int N_pre_input, int N_post_input, vector<int> &C_i, vector<int> &C_j, vector<double> &K_ij, vector<double> &D_ij){

	// read parameter
	synapses_type = synapses_type_input;
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

	// default settting
	STD = false; 
	STD_on_step = -1;
	inh_STDP = false;
	
	// parameter-dependent initialisation
	init();

}


void ChemicalSynapses::init(int synapses_type_input, int j_post, int N_post_input, double K_ext_input, int Num_ext_input, vector<double> &rate_ext_t_input, int ia_input, int ib_input){

	// Initialise chemical synapses for simulating external neuron population
	synapses_type = synapses_type_input;
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

	//
	STD = false; 
	inh_STDP = false;
	
	// parameter-dependent initialisation
	init();
}



void ChemicalSynapses::init(){
	// parameter-dependent initialisation

	
	// Initialise chemical synapse parameters
	if (synapses_type == 0){
		tau_decay = tau_decay_AMPA;
		steps_trans = int(round(Dt_trans_AMPA / dt));
	}
	else if (synapses_type == 1){	
		tau_decay = tau_decay_GABA;
		steps_trans = int(round(Dt_trans_GABA / dt));
	}
	else if (synapses_type == 2){
		tau_decay = tau_decay_NMDA;
		steps_trans = int(round(Dt_trans_NMDA / dt));
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


	// Initialize pre- and post-synaptic dynamic variables
	buffer_steps = max_delay_steps + steps_trans + 1; // the +1 is vital
	s.assign(N_pre, 0);
	gs_sum.assign(N_post, 0);
	d_gs_sum_buffer.resize(buffer_steps);
	for (int i = 0; i < buffer_steps; ++i){
		d_gs_sum_buffer[i].assign(N_post, 0);
	}
	I.assign(N_post, 0);
	
	
	// transmitter_strength
	K_trans.assign(N_pre, 1.0 / steps_trans); // be careful! 1 / transmitter steps gives zero (int)!!
	trans_left.assign(N_pre, 0);

	// Initialise exp_step
	exp_step = exp(-dt / tau_decay); // single step
	


}


void ChemicalSynapses::update(int step_current){

	
	if (pop_ind_pre >= 0){
		//
		if (STD_on_step == step_current){STD = true;}
		if (inh_STDP_on_step == step_current){inh_STDP = true;}
		
		// update pre-synaptic dynamics
		for (unsigned int i = 0; i < spikes_pre->size(); ++i){ // add spikes (transmitter release)
			trans_left[spikes_pre->at(i)] += steps_trans;
		}
		
		// short-term depression
		if (STD == true){
			for (unsigned int i = 0; i < spikes_pre->size(); ++i){ 
				K_trans[spikes_pre->at(i)] = 1.0 / steps_trans * f_ves[i];
			}
		}
		
		// inhibitory STDP
		if (inh_STDP == true){
			// update x_trace
			for (unsigned int i = 0; i < spikes_pre->size(); ++i){
				x_trace_pre[spikes_pre->at(i)] += 1.0;
			}
			for (unsigned int j = 0; j < spikes_post->size(); ++j){
				x_trace_post[spikes_post->at(j)] += 1.0;
			}
			// update K
			int i_pre, j_post;
			for (unsigned int ind_spike = 0; ind_spike < spikes_pre->size(); ++ind_spike){
				i_pre = spikes_pre->at(ind_spike);
				for (unsigned int syn_ind = 0; syn_ind < C[i_pre].size(); ++syn_ind){
					j_post = C[i_pre][syn_ind];
					K[i_pre][syn_ind] += eta_STDP * ( x_trace_post[j_post] - alpha_STDP );
				}
			}
			int syn_ind;
			for (unsigned int ind_spike = 0; ind_spike < spikes_post->size(); ++ind_spike){
				j_post = spikes_post->at(ind_spike);
				for (unsigned int ind = 0; ind < j_2_i[j_post].size(); ++ind){
					// j_2_i and j_2_syn_ind together serve as the "inverse function" of j_post = C[i_pre][syn_ind]
					i_pre = j_2_i[j_post][ind];
					syn_ind = j_2_syn_ind[j_post][ind];
					K[i_pre][syn_ind] += eta_STDP * x_trace_pre[i_pre];
				}
			}
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
	} // if pop_ind_pre >=0

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
	}
	// immediately reset the current buffer to zeros!!
	fill(d_gs_sum_buffer[t_ring].begin(), d_gs_sum_buffer[t_ring].end(), 0.0);
	

	

	// Calculate chemical currents
	// need V from post population!!
	if (synapses_type == 0){ //AMPA
		for (int j = 0; j < N_post; ++j){
			I[j] = -gs_sum[j] * (V_post->at(j) - V_ex);
		}
	}	
	else if (synapses_type == 1){ //GABA
		for (int j = 0; j < N_post; ++j){
			I[j] = -gs_sum[j] * (V_post->at(j) - V_in);
			// For inhibition, every equation is in the same form as excitation. 
			// Only "V_in" encodes its inhibitory nature.
		}
	}
	else if (synapses_type == 2){ //NMDA
		for (int j = 0; j < N_post; ++j){
			I[j] = -gs_sum[j] * B[(int)round((V_post->at(j) - B_V_min) / B_dV)] * (V_post->at(j) - V_ex);
		}
	}

	if (pop_ind_pre >= 0){
		// update short-term depression
		if (STD == true){
			for (unsigned int i = 0; i < spikes_pre->size(); ++i){
				f_ves[spikes_pre->at(i)] *= 1.0 - p_ves; // decrease
			}
			for (int i = 0; i < N_pre; ++i){
				f_ves[i] = 1.0 - exp_ves * (1.0 - f_ves[i]); // decay to 1.0
			}
		}
		// update inhibitory STDP 
		if (inh_STDP == true){
			for (int i = 0; i < N_pre; ++i){
				x_trace_pre[i] *= exp_step_STDP;
			}
			for (int j = 0; j < N_post; ++j){
				x_trace_post[j] *= exp_step_STDP;
			}
		}
		// decay pre-synaptic dynamics
		for (int i = 0; i < N_pre; ++i){ s[i] *= exp_step; };
	} 
	

	// decay post-synaptic dynamics
	for (int j = 0; j < N_post; ++j){ gs_sum[j] *= exp_step; };

	// sample data
	sample_data(step_current);
	
	//
	record_stats();
	
}// update



void ChemicalSynapses::add_sampling(vector<int> sample_neurons_input, vector<bool> sample_time_points_input){
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

void ChemicalSynapses::add_short_term_depression(int STD_on_step_input){
	if (synapses_type != 0){
		cout << "Warning: initializing STD on non-AMPA synapses!" << endl;
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

void ChemicalSynapses::add_inh_STDP(int inh_STDP_on_step_input){
	if (synapses_type != 1){
		cout << "Warning: initializing inhibitory STDP on non-GABA synapses!" << endl;
	}
	inh_STDP_on_step = inh_STDP_on_step_input;
	if (inh_STDP_on_step == 0){
		inh_STDP = true;
	}
	
	x_trace_pre.assign(N_pre, 0.0);
	x_trace_post.assign(N_post, 0.0);

	tau_STDP = 20; // ms
	exp_step_STDP = exp(-dt / tau_STDP);
	eta_STDP = 0.01; // learning rate, 0.0001 is the published value but requires 60min of simulation
	rho_0_STDP = 0.003; //kHz
	alpha_STDP = 2.0 * rho_0_STDP * tau_STDP; // depression factor
	
	// j_2_i and j_2_syn_ind
	int j_post;
	for (int i_pre = 0; i_pre < N_pre; ++i_pre){ 
		for (unsigned int syn_ind = 0; syn_ind < C[i_pre].size(); ++syn_ind){
			j_post = C[i_pre][syn_ind];
			j_2_i[j_post].push_back( i_pre );
			j_2_syn_ind[j_post].push_back( int(syn_ind) );
		}
	}
}	


void ChemicalSynapses::sample_data(int step_current){
	if (!sample_neurons.empty()){
		if (sample_time_points[step_current]){ // push_back is amazing
			for (unsigned int i = 0; i < sample_neurons.size(); ++i){ // performance issue when sampling many neurons?
				int ind_temp = sample_neurons[i];
				sample[i].push_back( I[ind_temp] );
			}
		}
	}
}



void ChemicalSynapses::set_para(string para_str, char delim){
	if (!para_str.empty()){
		istringstream para(para_str);
		string para_name, para_value_str, line_str; 
		double para_value;
		while (getline(para, line_str)){
			istringstream line_ss(line_str);
			getline(line_ss, para_name, delim); // get parameter name
			getline(line_ss, para_value_str, delim); // get parameter value (assume double)
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


string ChemicalSynapses::dump_para(char delim){
	stringstream dump;


	dump << "pop_ind_pre" << delim << pop_ind_pre << delim << endl;
	dump << "pop_ind_post" << delim << pop_ind_post << delim << endl;
	dump << "synapses_type" << delim << synapses_type << delim << endl;

	dump << "V_ex" << delim << V_ex << delim << endl;
	dump << "V_in" << delim << V_in << delim << endl;

	dump << "seed" << delim << my_seed << delim << endl;
	
	if (synapses_type == 0){
		dump << "Dt_trans_AMPA" << delim << Dt_trans_AMPA << delim << endl;
		dump << "tau_decay_AMPA" << delim << tau_decay_AMPA << delim << endl;
	}
	else if (synapses_type == 1){
		dump << "Dt_trans_GABA" << delim << Dt_trans_GABA << delim << endl;
		dump << "tau_decay_GABA" << delim << tau_decay_GABA << delim << endl;
	}
	else if (synapses_type == 2){
		dump << "Dt_trans_NMDA" << delim << Dt_trans_NMDA << delim << endl;
		dump << "tau_decay_NMDA" << delim << tau_decay_NMDA << delim << endl;
	}

	return dump.str();
}

void ChemicalSynapses::start_stats_record(){
	stats_record = true;
	I_mean.reserve(step_tot);
	I_std.reserve(step_tot);
}


void ChemicalSynapses::output_results(ofstream& output_file, char delim, char indicator){
	// SYND001 # synapse parameters
	// count number of variables
	stringstream dump_count;
	string para_str = dump_para(delim);
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
		output_file << pop_ind_pre << delim << pop_ind_post << delim << synapses_type << delim << sample_neurons.size() << delim << endl;
		write2file(output_file, delim, sample); // 2D matrix
	}

	
	// SYND003 # currents mean and std
	if (stats_record){
		output_file << indicator << " SYND003" << endl;
		output_file << pop_ind_pre << delim << pop_ind_post << delim << synapses_type << delim << endl;
		write2file(output_file, delim, I_mean);
		write2file(output_file, delim, I_std);
	}
	
	// tmp data
	if (tmp_data.size() != 0){
		output_file << indicator << " SYND004" << endl;
		output_file << pop_ind_pre << delim << pop_ind_post << delim << synapses_type << delim << tmp_data.size() << delim << endl;
		write2file(output_file, delim, tmp_data);	
	}
}



void ChemicalSynapses::recv_pop_data(vector<Neurons> &NeuronPopArray){
	// get current spikes from pre-pop
	if (pop_ind_pre >= 0){
		spikes_pre = &(NeuronPopArray[pop_ind_pre].spikes_current); // This might be problematic!!!
		if (inh_STDP == true){
			spikes_post = &(NeuronPopArray[pop_ind_post].spikes_current);
		}
	}
	// get current V from post-pop
	V_post = &(NeuronPopArray[pop_ind_post].V); // This is problematic!!!
	
}


void ChemicalSynapses::send_pop_data(vector<Neurons> &NeuronPopArray){
	
	// send currents to post-pop
	//AMPA
	if (synapses_type == 0){
		for (int j = 0; j < N_post; ++j){
			NeuronPopArray[pop_ind_post].I_AMPA[j] += I[j];
		}
	}
	//GABA
	else if (synapses_type == 1){
		for (int j = 0; j < N_post; ++j){
			NeuronPopArray[pop_ind_post].I_GABA[j] += I[j];
		}
	}
	//NMDA
	else if (synapses_type == 2){
		for (int j = 0; j < N_post; ++j){
			NeuronPopArray[pop_ind_post].I_NMDA[j] += I[j];
		}
	}
}


// Use function templates when you want to perform the same action on types that can be different.
// Use function overloading when you want to apply different operations depending on the type.
// In this case, just save yourself the trouble and use overloading.
void ChemicalSynapses::write2file(ofstream& output_file, char delim, vector< vector<int> >& v){
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



void ChemicalSynapses::write2file(ofstream& output_file, char delim, vector< vector<double> >& v){
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


void ChemicalSynapses::write2file(ofstream& output_file, char delim, vector<int>& v){
	if (!v.empty()){
		//for (int f : v){ output_file << f << delim; } // range-based "for" in C++11
		for (unsigned int i = 0; i < v.size(); ++i){
			output_file << v[i] << delim;
		}
		output_file << endl;
	}
	else {output_file << " " << endl;}
}

void ChemicalSynapses::write2file(ofstream& output_file, char delim, vector<double>& v){
	if (!v.empty()){
		//for (int f : v){ output_file << f << delim; } // range-based "for" in C++11
		for (unsigned int i = 0; i < v.size(); ++i){
			output_file << v[i] << delim;
		}
		output_file << endl;
	}
	else {output_file << " " << endl;}
}



void ChemicalSynapses::record_stats(){
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
