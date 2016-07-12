#include <stdlib.h>
#include <iomanip>
#include <cmath>
#include <algorithm> // for transform
#include <stdio.h> // for printf
#include <time.h>       /* time */
#include <sstream>  // stringstream is input and output
#include "NeuroPop.h"
#include <functional> // for bind(), plus

NeuroPop::NeuroPop(const int pop_ind_input, const int N_input, const double dt_input, const int step_tot_input, const char delim_input, const char indicator_input)
{
	pop_ind = pop_ind_input;
	N = N_input;
	dt = dt_input;
	step_tot = step_tot_input; // this parameter is designed to be self-adapting (step_killed), so should be any other stuff that relies on it!!
	delim = delim_input;
	indicator = indicator_input;
	
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

	spike_hist_tot.reserve(step_tot*50); // reserve!
	num_ref_pop.reserve(step_tot); // reserve and push_back so that it won't be affected by adapting step_tot
	num_spikes_pop.reserve(step_tot); // reserve and push_back so that it won't be affected by adapting step_tot

	// Random seed (random engine should be feed with DIFFERENT seed at every implementation)
	random_device rd; // random number from operating system for seed
	my_seed = rd(); // record seed
	//my_seed = 321;
	//cout << "My_seed is: " << my_seed << endl;

	// Runaway killer is initially a sleeper agent
	killer.license = false;
	killer.runaway_killed = false;
	killer.step_killed = -1;
	
	//
	stats.record = false;
	LFP.record = false;
	spike_freq_adpt = false;
	
	// perturbation
	step_perturb = -1;
	spike_removed = -1;

}

const vector< int > & NeuroPop::get_spikes_current() 
{
	return spikes_current;
}

const vector< double > & NeuroPop::get_V()
{
	return V;
}

const bool & NeuroPop::get_runaway_killed()
{
	return killer.runaway_killed;
};

void NeuroPop::recv_I(vector<double>& I, const int pop_ind_pre, const int syn_type)
{
	if (pop_ind_pre == -1){ // if noisy external currents, always send to I_ext regardless of the synapse type
		for (int j = 0; j < N; ++j){
			I_ext[j] += I[j];
		}
	}
	// AMPA
	else if (syn_type == 0){
		for (int j = 0; j < N; ++j){
			I_AMPA[j] += I[j];
		}
	}
	//GABA
	else if (syn_type == 1){
		for (int j = 0; j < N; ++j){
			I_GABA[j] += I[j];
		}
	}
	//NMDA
	else if (syn_type == 2){
		for (int j = 0; j < N; ++j){
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
	
	stats.I_AMPA_acc.assign(N, 0.0);
	stats.I_AMPA_time_avg.assign(N, 0.0);
	stats.I_NMDA_acc.assign(N, 0.0);
	stats.I_NMDA_time_avg.assign(N, 0.0);
	stats.I_GABA_acc.assign(N, 0.0);
	stats.I_GABA_time_avg.assign(N, 0.0);
	
	stats.IE_ratio.assign(N, 0.0);
}

void NeuroPop::start_LFP_record(vector <vector<bool> >& LFP_neurons_input){
	if (int(LFP_neurons_input[0].size()) != N){
		cout << "start_LFP_record failed: LFP_neurons should be 1-by-N logical vector!" << endl;
	}
	else{
		LFP.record = true;
		int n_LFP = int(LFP_neurons_input.size());
		LFP.neurons.resize(n_LFP);
		LFP.data.resize(n_LFP);
		for (int ind = 0; ind < n_LFP; ++ind){
			LFP.neurons[ind] = LFP_neurons_input[ind];
			LFP.data[ind].reserve(step_tot);
		}
	}	
}



void NeuroPop::set_para(string para_str){
	if (!para_str.empty()){
		istringstream para(para_str);
		string para_name, para_value_str; 
		double para_value;
		while (getline(para, para_name, delim)){
			getline(para, para_value_str, delim); // get parameter value (assume double)
			stringstream(para_value_str) >> para_value; // from string to numerical value
			if (para_name.find("Cm") != string::npos){Cm = para_value;}
			else if (para_name.find("tau_ref") != string::npos){tau_ref = para_value;}
			else if (para_name.find("V_rt") != string::npos){V_rt = para_value;}
			else if (para_name.find("V_lk") != string::npos){V_lk = para_value;}
			else if (para_name.find("V_th") != string::npos){V_th = para_value;}
			else if (para_name.find("g_lk") != string::npos){g_lk = para_value;}
			else {cout << "Unrecognized parameter: " << para_name << endl;}
		}
	}
	// re-initialise population!
	init();
}

string NeuroPop::dump_para(){
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


void NeuroPop::random_V(const double p){
	// Generate random initial condition for V.
	// Generate uniform random distribution 
	cout << "Function NeuroPop::random_V(double p) is deprecated!" << endl;
	if (p < 1.0){
		gen.seed(my_seed);// reseed random engine!
		uniform_real_distribution<double> uniform_dis(0.0, 1.0);
		auto ZeroOne = bind(uniform_dis,gen);

		for (int i = 0; i < N; ++i){
			// Generate random number.
			//ZeroOne = uniform_dis(gen); // range is [0,1]
			V[i] = V_rt + (V_th - V_rt) / (1.0 - p)*ZeroOne();
		}
	}
	else {cout << "Initial firing rate cannot be 100%!" << endl;}
}

void NeuroPop::set_init_condition(const double r_V0, const double p_fire){
	// Set V to be uniformly distributed between [V_rt, V_rt + (V_th - V_rt)*r_V0]
	// And then randomly set some of them above firing threshold according to p_fire
	gen.seed(my_seed);// reseed random engine!
	uniform_real_distribution<double> uniform_dis(0.0, 1.0);
	auto ZeroOne = bind(uniform_dis,gen);
	for (int i = 0; i < N; ++i){
		if (ZeroOne() < p_fire){
			V[i] = V_th + 1.0; // above threshold for firing
		}
		else{V[i] = V_rt + (V_th - V_rt) *r_V0 * ZeroOne();}
	}
}


void NeuroPop::update_spikes(const int step_current){
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
	for (int i = 0; i < N; ++i){
		if (ref_step_left[i] == 0 && V[i] >= V_th){
			spikes_current.push_back(i); // record firing neurons
			V[i] = V_rt; // reset potential
			ref_step_left[i] = ref_steps; // steps left for being refractory
			spike_counter += 1;
		}
	}
	
	
	// perturbation: remove the last spike
	if (step_current == step_perturb){
		if (spikes_current.size() > 0){
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
	if (spike_counter > 0){
		// .end(): ONE ELEMENT PAST to the last location
		// see http://www.cs.northwestern.edu/~riesbeck/programming/c++/stl-iterators.html
		copy(spikes_current.begin(), spikes_current.end(), back_inserter(spike_hist_tot));
	}

	// Refraction count-down and record number of refractory neurons
	int num_ref_temp = 0;
	for (int i = 0; i < N; ++i){
		if (ref_step_left[i] > 0){
			ref_step_left[i] -= 1;
			num_ref_temp += 1;
		}
	}
	num_ref_pop.push_back(num_ref_temp);


	// runaway check
	runaway_check(step_current);
	
}

void NeuroPop::generate_I_ext(const int step_current){
	
	// Gaussian white external currents
	if (I_ext_mean.size() != 0){
		if (I_ext_std.size() != 0){
			double one_on_sqrt_dt = 1.0/sqrt(dt); // Here sqrt_dt is on the denominator because I_ext will be multiplied by dt later. 
			// Gaussian random generator
			gen.seed(my_seed+step_current);// reseed random engine!
			normal_distribution<double> nrm_dist(0.0, 1.0);
			auto gaus = bind(nrm_dist,gen);

			for (int i = 0; i < N; ++i){ 
				I_ext[i] += I_ext_mean[i] + gaus() * I_ext_std[i] * one_on_sqrt_dt; // be careful about the sqrt(dt) term (Wiener Process)
			}
		}
		else{ 
			for (int i = 0; i < N; ++i){ 
				I_ext[i] += I_ext_mean[i]; 
			} 
		}
	}
	
	// Gaussian white external conductance
	if (g_ext_mean.size() != 0){
		if (g_ext_std.size() != 0){
			double one_on_sqrt_dt = 1.0/sqrt(dt); // Here sqrt_dt is on the denominator because I_ext will be multiplied by dt later. 
			// Gaussian random generator
			gen.seed(my_seed+step_current+step_tot);// reseed random engine!
			normal_distribution<double> nrm_dist(0.0, 1.0);
			auto gaus = bind(nrm_dist,gen);

			for (int i = 0; i < N; ++i){ 
				I_ext[i] += -(g_ext_mean[i] + gaus() * g_ext_std[i] * one_on_sqrt_dt) * (V[i] - V_ext); // be careful about the sqrt(dt) term (Wiener Process)
			}
		}
		else{ 
			for (int i = 0; i < N; ++i){ 
				I_ext[i] += -g_ext_mean[i] * (V[i] - V_ext);
			} 
		}
	}
}
void NeuroPop::update_V(const int step_current){
	// This function updates menbrane potentials for non-refractory neurons

	// Generate external currents
	generate_I_ext(step_current);
	
	// potassium conductance for spike-frequency adaptation
	if (spike_freq_adpt == true){
		for (unsigned int ind = 0; ind < spikes_current.size(); ++ind){
			g_K[ spikes_current[ind] ] += dg_K;
		}
		for (int i = 0; i < N; ++i){
			g_K[i] *= exp_K_step;
			I_K[i] = -g_K[i] * (V[i] - V_K);
		}
	}

	// Collect Currents from all pre-synapses (for MPI job)!!!!!!!!!!!!!!!



	// Data sampling, which must be done here!
	// sample_data(step_current);  // this is deprecated due to poor memory performance
	output_sampled_data_real_time_HDF5(step_current);
	// we want to sample the V[] before it's updated!



	// update menbrane potentials
	double Vdot;
	for (int i = 0; i < N; ++i){
		I_input[i] = I_AMPA[i] + I_GABA[i] + I_NMDA[i] + I_GJ[i] + I_ext[i] + I_K[i];
		if (ref_step_left[i] == 0){ // Only update the non-refractory neurons
			// leaky current
			I_leak[i] = -g_lk * (V[i] - V_lk); 
			// using simple Euler method
			Vdot = (I_leak[i] + I_input[i])/Cm;
			V[i] += Vdot * dt;
			// Note that delta-function coupling is very different from the above conductance-based model!
		}
	}

	// record mean and std of membrane potentials
	record_stats(step_current);
	record_LFP();
	
}


void NeuroPop::add_sampling_real_time(const vector<int>& sample_neurons_input, const vector<bool>& sample_type_input, const vector<bool>& sample_time_points_input, string sample_file_name_input){
	sample.neurons = sample_neurons_input;
	sample.type = sample_type_input;
	sample.time_points = sample_time_points_input;
	
	sample.file_name = sample_file_name_input.append(to_string(pop_ind)).append(".ygout_samp");
	sample.file.open(sample.file_name);
	
	sample.N_steps = 0;
	for (int tt = 0; tt < step_tot; ++tt){
		if (sample.time_points[tt]){sample.N_steps += 1;}
	}
	sample.N_neurons = sample.neurons.size();
}


void NeuroPop::add_sampling_real_time_HDF5(vector<int>& sample_neurons_input, vector<bool>& sample_type_input,  vector<bool>& sample_time_points_input, string sample_file_name_input){
	sample.neurons = sample_neurons_input;
	sample.type = sample_type_input;
	sample.time_points = sample_time_points_input;
	
	sample.N_steps = 0;
	for (int tt = 0; tt < step_tot; ++tt){
		if (sample.time_points[tt]){sample.N_steps += 1;}
	}
	sample.N_neurons = sample.neurons.size();
	
	sample.file_name = sample_file_name_input.append("_");
	sample.file_name.append(to_string(pop_ind));
	
	sample.file_name.append("_out.h5");
	// Create a new file using default properties. 
	cout << "Creating HDF5 output file...\n";
	sample.file_HDF5 = new H5File( sample.file_name.c_str(), H5F_ACC_TRUNC );
	// dataset dimensions
	hsize_t dims[2]; 
	dims[1] = sample.N_neurons;
	dims[0] = sample.N_steps;
	// file dataspace
	DataSpace fspace(2, dims); 
	// create datasets
	if(sample.type[0]){
		sample.V_dataset = sample.file_HDF5->createDataSet("V", PredType::NATIVE_DOUBLE, fspace);
	}
	if(sample.type[1]){
		sample.I_leak_dataset = sample.file_HDF5->createDataSet("I_leak", PredType::NATIVE_DOUBLE, fspace);
	}
	if(sample.type[2]){
		sample.I_AMPA_dataset = sample.file_HDF5->createDataSet("I_AMPA", PredType::NATIVE_DOUBLE, fspace);
	}
	if(sample.type[3]){
		sample.I_GABA_dataset = sample.file_HDF5->createDataSet("I_GABA", PredType::NATIVE_DOUBLE, fspace);
	}
	if(sample.type[4]){
		sample.I_NMDA_dataset = sample.file_HDF5->createDataSet("I_NMDA", PredType::NATIVE_DOUBLE, fspace);
	}
	if(sample.type[5]){
		sample.I_GJ_dataset = sample.file_HDF5->createDataSet("I_GJ", PredType::NATIVE_DOUBLE, fspace);
	}
	if(sample.type[6]){
		sample.I_ext_dataset = sample.file_HDF5->createDataSet("I_ext", PredType::NATIVE_DOUBLE, fspace);
	}
	if(sample.type[7]){
		sample.I_K_dataset = sample.file_HDF5->createDataSet("I_K", PredType::NATIVE_DOUBLE, fspace);
	}
}

void NeuroPop::output_sampled_data_real_time_HDF5(const int step_current){
	if (!sample.neurons.empty()){
		if (sample.time_points[step_current]){ // push_back is amazing

			// use double *temp_data=new double[sample.N_neurons]; and later delete[] tmp_data?
			vector<double> temp_data;
			temp_data.resize(sample.N_neurons);
			
			int ind_temp;	
			if (sample.type[0]){	// "V"
				fill(temp_data.begin(), temp_data.end(), 0);
				// Collect the data to write in a temp vector
				for (int i = 0; i < sample.N_neurons; ++i){ // performance issue when sampling many neurons?
					ind_temp = sample.neurons[i];
					temp_data[i] = V[ind_temp];
				}
				append_vector_to_matrix_HDF5(sample.V_dataset, temp_data,  sample.ctr);
			}
			if (sample.type[1]){	// "I_leak"
				fill(temp_data.begin(), temp_data.end(), 0);
				// Collect the data to write in a temp vector
				for (int i = 0; i < sample.N_neurons; ++i){ // performance issue when sampling many neurons?
					ind_temp = sample.neurons[i];
					temp_data[i] = I_leak[ind_temp];
				}
				append_vector_to_matrix_HDF5(sample.I_leak_dataset, temp_data,  sample.ctr);
			}	
			if (sample.type[2]){	// "I_AMPA"
				fill(temp_data.begin(), temp_data.end(), 0);
				// Collect the data to write in a temp vector
				for (int i = 0; i < sample.N_neurons; ++i){ // performance issue when sampling many neurons?
					ind_temp = sample.neurons[i];
					temp_data[i] = I_AMPA[ind_temp];
				}
				append_vector_to_matrix_HDF5(sample.I_AMPA_dataset, temp_data,  sample.ctr);
			}	
			if (sample.type[3]){	// "I_GABA"
				fill(temp_data.begin(), temp_data.end(), 0);
				// Collect the data to write in a temp vector
				for (int i = 0; i < sample.N_neurons; ++i){ // performance issue when sampling many neurons?
					ind_temp = sample.neurons[i];
					temp_data[i] = I_GABA[ind_temp];
				}
				append_vector_to_matrix_HDF5(sample.I_GABA_dataset, temp_data,  sample.ctr);	
			}	
			if (sample.type[4]){	// "I_NMDA"
				fill(temp_data.begin(), temp_data.end(), 0);
				// Collect the data to write in a temp vector
				for (int i = 0; i < sample.N_neurons; ++i){ // performance issue when sampling many neurons?
					ind_temp = sample.neurons[i];
					temp_data[i] = I_NMDA[ind_temp];
				}
				append_vector_to_matrix_HDF5(sample.I_NMDA_dataset, temp_data,  sample.ctr);	
			}	
			if (sample.type[5]){	// "I_GJ"
				fill(temp_data.begin(), temp_data.end(), 0);
				// Collect the data to write in a temp vector
				for (int i = 0; i < sample.N_neurons; ++i){ // performance issue when sampling many neurons?
					ind_temp = sample.neurons[i];
					temp_data[i] = I_GJ[ind_temp];
				}
				append_vector_to_matrix_HDF5(sample.I_GJ_dataset, temp_data,  sample.ctr);		
			}	
			if (sample.type[6]){	// "I_ext"
				fill(temp_data.begin(), temp_data.end(), 0);
				// Collect the data to write in a temp vector
				for (int i = 0; i < sample.N_neurons; ++i){ // performance issue when sampling many neurons?
					ind_temp = sample.neurons[i];
					temp_data[i] = I_ext[ind_temp];
				}
				append_vector_to_matrix_HDF5(sample.I_ext_dataset, temp_data,  sample.ctr);
			}	
			if (sample.type[7]){	// "I_K"
				fill(temp_data.begin(), temp_data.end(), 0);
				// Collect the data to write in a temp vector
				for (int i = 0; i < sample.N_neurons; ++i){ // performance issue when sampling many neurons?
					ind_temp = sample.neurons[i];
					temp_data[i]= I_K[ind_temp];
				}
				append_vector_to_matrix_HDF5(sample.I_K_dataset, temp_data,  sample.ctr);
			}	
			sample.ctr++;
		}
	}
}




void NeuroPop::output_sampled_data_real_time(const int step_current){
	
	if (!sample.neurons.empty() && step_current == 0){
		
		sample.file << indicator << " POPD006" << endl;
		sample.file << pop_ind << delim << sample.N_neurons << delim << sample.N_steps << delim << endl;
  		vector< string > data_types = { "V", "I_leak", "I_AMPA", "I_GABA", "I_NMDA", "I_GJ", "I_ext", "I_K" };
		for (unsigned int tt = 0; tt < data_types.size(); ++tt){
			if (sample.type[tt]){sample.file << data_types[tt] << delim; }
		}
		sample.file << endl;
	}
	
	if (!sample.neurons.empty()){
		if (sample.time_points[step_current]){ // push_back is amazing
			for (int i = 0; i < sample.N_neurons; ++i){ // performance issue when sampling many neurons?
				int ind_temp = sample.neurons[i];
				if (sample.type[0]){sample.file << V[ind_temp] << delim;}
				if (sample.type[1]){sample.file << I_leak[ind_temp] << delim;}
				if (sample.type[2]){sample.file << I_AMPA[ind_temp] << delim;}
				if (sample.type[3]){sample.file << I_GABA[ind_temp] << delim;}
				if (sample.type[4]){sample.file << I_NMDA[ind_temp] << delim;}
				if (sample.type[5]){sample.file << I_GJ[ind_temp] << delim;}
				if (sample.type[6]){sample.file << I_ext[ind_temp] << delim;}
				if (sample.type[7]){sample.file << I_K[ind_temp] << delim;}
				sample.file << endl;
			}
		}
	}

}



void NeuroPop::add_sampling(vector<int>& sample_neurons_input, vector<bool>& sample_type_input, vector<bool>& sample_time_points_input){
	sample.neurons = sample_neurons_input;
	sample.type = sample_type_input;
	sample.time_points = sample_time_points_input;
	
	
	// initialise
	sample.N_steps = 0;
	for (int tt = 0; tt < step_tot; ++tt){
		if (sample.time_points[tt]){sample.N_steps += 1;}
	}
	sample.N_neurons = sample.neurons.size();

	int sample_type_tot = sample.type.size(); // 8 different data types
	
	
	sample.data.resize(sample_type_tot); 
	for (int c = 0; c < sample_type_tot; ++c){
		if (sample.type[c]){
			sample.data[c].resize(sample.N_neurons);
			for (int i = 0; i < sample.N_neurons; ++i){
				sample.data[c][i].reserve(sample.N_neurons); // reserve and push_back so that it won't be affected by adapting step_tot
			}
		}
	}
}



void NeuroPop::sample_data(const int step_current){
	
	if (!sample.neurons.empty()){
		if (sample.time_points[step_current]){ // push_back is amazing
			for (int i = 0; i < sample.N_neurons; ++i){ // performance issue when sampling many neurons?
				int ind_temp = sample.neurons[i];
				if (sample.type[0]){sample.data[0][i].push_back( V[ind_temp] );}
				if (sample.type[1]){sample.data[1][i].push_back( I_leak[ind_temp] );}
				if (sample.type[2]){sample.data[2][i].push_back( I_AMPA[ind_temp] );}
				if (sample.type[3]){sample.data[3][i].push_back( I_GABA[ind_temp] );}
				if (sample.type[4]){sample.data[4][i].push_back( I_NMDA[ind_temp] );}
				if (sample.type[5]){sample.data[5][i].push_back( I_GJ[ind_temp] );}
				if (sample.type[6]){sample.data[6][i].push_back( I_ext[ind_temp] );}
				if (sample.type[7]){sample.data[7][i].push_back( I_K[ind_temp] );}
			}
		}
	}

}


void NeuroPop::set_gaussian_I_ext(vector<double>& mean, vector<double>& std){
	I_ext_mean = mean;
	I_ext_std = std;
	
	double max_std = *max_element(I_ext_std.begin(), I_ext_std.end());
	if (max_std == 0.0){
		I_ext_std.resize(0);
	}
}

void NeuroPop::set_gaussian_g_ext(vector<double>& mean, vector<double>& std){
	g_ext_mean = mean;
	g_ext_std = std;
	
	double max_std = *max_element(g_ext_std.begin(), g_ext_std.end());
	if (max_std == 0.0){
		g_ext_std.resize(0);
	}
}


void NeuroPop::add_perturbation(const int step_perturb_input){
	step_perturb = step_perturb_input;
}

void NeuroPop::add_spike_freq_adpt(){
	spike_freq_adpt = true;
	exp_K_step = exp( -dt / tau_K );
	g_K.assign(N, 0.0);
}


void NeuroPop::init_runaway_killer(const double min_ms, const double Hz, const double Hz_ms)
{
	killer.min_pop_size = 100;
	if (N > killer.min_pop_size){
		killer.license = true;
		killer.min_steps = int(round(min_ms / dt));
		killer.runaway_Hz = Hz;
		killer.Hz_steps = int(round(Hz_ms / dt));
	}
}

void NeuroPop::runaway_check(const int step_current)
{
	if (killer.license == true && killer.runaway_killed == false && step_current > killer.min_steps && step_current > killer.Hz_steps){
		// find mean value of num_ref over the last runaway_steps
		vector<int>::const_iterator first, last;
		// first element to be accumulated
		first = num_spikes_pop.begin() + (step_current - killer.Hz_steps + 1); 
		// one element pass the last element to be accumulated
		last = num_spikes_pop.begin() + (step_current + 1); 
		double mean_Hz = accumulate(first, last, 0.0) / (killer.Hz_steps * dt * 0.001 * N); // 0.001 for converting from ms to sec.
		//be careful!! accumulate range is : [first,last)
		if (mean_Hz >= killer.runaway_Hz){
			killer.runaway_killed = true;
			killer.step_killed = step_current;
			cout << "warning: runaway killed at " << step_current*dt << " (ms) in population" << pop_ind << flush;
			cout << "\t with firing rate at " << mean_Hz << " Hz."<< flush;
		}
	}
}


void NeuroPop::output_results(H5File& file){

	// new group
	string group_name = "/pop_result_";
	group_name.append(to_string(pop_ind));
	Group group_pop = file.createGroup(group_name);

	write_vector_HDF5(group_pop, spike_hist_tot, string("spike_hist_tot"));
	write_vector_HDF5(group_pop, num_spikes_pop, string("num_spikes_pop"));
	write_vector_HDF5(group_pop, num_ref_pop, string("num_ref_pop"));
	
	write_string_HDF5(group_pop, dump_para(), string("pop_para"));
		
	if (stats.record){
		write_vector_HDF5(group_pop, stats.V_mean, string("stats_V_mean"));
		write_vector_HDF5(group_pop, stats.V_std, string("stats_V_std"));
		write_vector_HDF5(group_pop, stats.I_input_mean, string("stats_I_input_mean"));
		write_vector_HDF5(group_pop, stats.I_input_std, string("stats_I_input_std"));
		write_vector_HDF5(group_pop, stats.IE_ratio, string("stats_IE_ratio"));
	}
	
	if (LFP.record){
		write_matrix_HDF5(group_pop, LFP.data, string("LFP_data"));
	}
	
}


void NeuroPop::output_results(ofstream& output_file){

	// POPD001 # spike history of neuron population
	output_file << indicator << " POPD001" << endl;
	output_file << pop_ind << delim << endl;
	write2file(output_file, spike_hist_tot);
	write2file(output_file, num_spikes_pop);
	write2file(output_file, num_ref_pop);

	// POPD002 # neuron parameters in the population
	stringstream dump_count;
	string para_str = dump_para();
	dump_count << para_str;
	string str_temp;
	int var_number = 0;
	while(getline(dump_count,str_temp)){++var_number;} // count number of variables
	output_file << indicator << " POPD002" << endl;
	output_file << pop_ind << delim << var_number << delim << endl;
	output_file << para_str;

	// POPD003 # membrane potential mean and std
	if (stats.record){
		output_file << indicator << " POPD003" << endl;
		output_file << pop_ind << delim << endl;
		write2file(output_file, stats.V_mean);
		write2file(output_file, stats.V_std);
		write2file(output_file, stats.I_input_mean);
		write2file(output_file, stats.I_input_std);
	}
	
	// POPD007 # local field potential
	if (LFP.record){
		output_file << indicator << " POPD007" << endl;
		output_file << pop_ind << delim << LFP.data.size() << delim << endl;
		write2file(output_file, LFP.data);
	}

	// POPD005 # E-I ratio for each neuron
	if (stats.record){
		output_file << indicator << " POPD005" << endl;
		output_file << pop_ind << delim << endl;
		write2file(output_file, stats.IE_ratio);
	}
	
	// SAMF001 # sampled data file name
	if (!sample.file_name.empty()){
		output_file << indicator << " SAMF001" << endl;
		output_file << sample.file_name << endl;
	}


	/* // This following output protocol is deprecated due to poor memeory performance
	// POPD004 # sampled neuron data
	if (!sample_neurons.empty()){
		output_file << indicator << " POPD004" << endl;
		output_file << pop_ind << delim << sample.neurons.size() << delim << endl;

  		vector< string > data_types = { "V", "I_leak", "I_AMPA", "I_GABA", "I_NMDA", "I_GJ", "I_ext", "I_K" };
		for (unsigned int tt = 0; tt < data_types.size(); ++tt){
			if (sample.type[tt] == true){output_file << data_types[tt] << delim;}
		}
		output_file << endl;


		for (unsigned int c = 0; c < sample.type.size(); ++c){
			if (!sample.data[c].empty()){
				NeuroPop::write2file(output_file, sample.data[c]); // 2D matrix
			}
		}
	}
	*/

}



void NeuroPop::record_LFP(){
	if (LFP.record){
		for (unsigned int ind = 0; ind < LFP.neurons.size(); ++ind){
			LFP.data[ind].push_back(0.0);
			for (int i = 0; i < N; ++i){
				if (LFP.neurons[ind][i]){
					LFP.data[ind].back() += abs(I_AMPA[i]) + abs(I_GABA[i]);
				}
			} 
		}
	}
}


void NeuroPop::record_stats(const int step_current){
	if (stats.record){
		// get mean
		double sum_mean_V = 0.0;
		double sum_mean_I = 0.0;
		for (unsigned int i = 0; i < V.size(); ++i){
			sum_mean_V += V[i];
			sum_mean_I += I_input[i];
		}
		double mean_tmp_V = sum_mean_V / double(V.size());
		double mean_tmp_I = sum_mean_I / double(V.size());
	
		// get std
		double sum_std_V = 0.0;
		double sum_std_I = 0.0;
		for (unsigned int i = 0; i < V.size(); ++i){
			sum_std_V += (V[i]-mean_tmp_V)*(V[i]-mean_tmp_V);
			sum_std_I += (I_input[i]-mean_tmp_I)*(I_input[i]-mean_tmp_I);
		}
		double std_tmp_V = sqrt( sum_std_V / double(V.size()));
		double std_tmp_I = sqrt( sum_std_I / double(V.size()));
	
		// record   
		stats.V_mean.push_back(mean_tmp_V);
		stats.V_std.push_back(std_tmp_V);
		stats.I_input_mean.push_back(mean_tmp_I);
		stats.I_input_std.push_back(std_tmp_I);
		
		// accumulate
		//for (unsigned int i = 0; i < N; ++i){ // this manual loop is slow, use transform()
		//	I_input_acc[i] += I_input[i];
		//}
		transform( stats.I_AMPA_acc.begin(), stats.I_AMPA_acc.end(), I_AMPA.begin(), stats.I_AMPA_acc.begin(), plus<double>() );
		transform( stats.I_NMDA_acc.begin(), stats.I_NMDA_acc.end(), I_NMDA.begin(), stats.I_NMDA_acc.begin(), plus<double>() );
		transform( stats.I_GABA_acc.begin(), stats.I_GABA_acc.end(), I_GABA.begin(), stats.I_GABA_acc.begin(), plus<double>() );
		// get time average for each neuron
		if (step_current == step_tot - 1){ // at the end of the last time step
			for (int i = 0; i < N; ++i){
				stats.I_AMPA_time_avg[i] = stats.I_AMPA_acc[i] / step_tot;
				stats.I_NMDA_time_avg[i] = stats.I_NMDA_acc[i] / step_tot;
				stats.I_GABA_time_avg[i] = stats.I_GABA_acc[i] / step_tot;
				// be careful here, for IE_ratio, I_ext is assumed to be always excitatory and I_GJ is not considered
				// also, the only source of I_ext is generated internally
				stats.IE_ratio[i] = stats.I_GABA_time_avg[i] / (stats.I_AMPA_time_avg[i] + stats.I_NMDA_time_avg[i] + I_ext_mean[i]);
			}
		}
	}
}



// Use function templates when you want to perform the same action on types that can be different.
// Use function overloading when you want to apply different operations depending on the type.
// In this case, just save yourself the trouble and use overloading.
void NeuroPop::write2file(ofstream& output_file,  vector< vector<int> >& v){
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



void NeuroPop::write2file(ofstream& output_file, vector< vector<double> >& v){
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


void NeuroPop::write2file(ofstream& output_file, const vector<int>& v){
	if (!v.empty()){
		//for (int f : v){ output_file << f << delim; } // range-based "for" in C++11
		for (unsigned int i = 0; i < v.size(); ++i){
			output_file << v[i] << delim;
		}
		output_file << endl;
	}
	else {output_file << " " << endl;}
}


void NeuroPop::write2file(ofstream& output_file, const vector<double>& v){
	if (!v.empty()){
		//for (int f : v){ output_file << f << delim; } // range-based "for" in C++11
		for (unsigned int i = 0; i < v.size(); ++i){
			output_file << v[i] << delim;
		}
		output_file << endl;
	}
	else {output_file << " " << endl;}
}

void NeuroPop::write_vector_HDF5(Group & group, const vector<int> & v, const string & v_name){
	hsize_t dims[1]; 
	dims[0] = v.size();
	DataSpace fspace(1, dims); 
	DataSet v_dataset = group.createDataSet(v_name, PredType::NATIVE_INT32, fspace);
	v_dataset.write( v.data(), PredType::NATIVE_INT32, fspace, fspace );	
}

void NeuroPop::write_vector_HDF5(Group & group, const vector<double> & v, const string &   v_name){
	hsize_t dims[1]; 
	dims[0] = v.size();
	DataSpace fspace(1, dims); 
	DataSet v_dataset = group.createDataSet(v_name, PredType::NATIVE_DOUBLE, fspace);
	v_dataset.write( v.data(), PredType::NATIVE_INT32, fspace, fspace );	
}

void NeuroPop::write_string_HDF5(Group & group, const string & s, const string &  s_name){
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

void NeuroPop::write_matrix_HDF5(Group & group, vector< vector<double> > & m, const string &  m_name){
	hsize_t dims[2]; 
	dims[0] = m.size();
	dims[1] = m[0].size();
	DataSpace fspace(2, dims); 
	DataSet m_dataset = group.createDataSet(m_name, PredType::NATIVE_DOUBLE, fspace);
	for (int i = 0; i < int(m.size()); i++){
		append_vector_to_matrix_HDF5(m_dataset, m[i],  i);
	}
}


void NeuroPop::append_vector_to_matrix_HDF5(DataSet & dataset_tmp, const vector<double> & v, const int colNum){
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


