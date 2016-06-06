#include <stdlib.h>
#include <iomanip>
#include <random> //<boost/random.hpp>
#include <cmath>
#include <algorithm> // for transform
#include <stdio.h> // for printf
#include <time.h>       /* time */

#include <fstream> // fstream : Stream class to both read and write from / to files
#include <sstream>  // stringstream is input and output


#include "Neurons.h"
#include <functional> // for bind(), plus
// no need to include what have been included in the header file

Neurons::Neurons(int pop_ind_input, int N_input, double dt_input, int step_tot_input, char delim_input, char indicator_input){

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
	// Leak conductance
	g_lk = 0.0167;   // (uS=miuSiemens), time constants=Cm/gL=15 ms!
	// spike-frequency adaptation parameter
	V_K = -85.0; // mV
	dg_K = 0.01; // (uS=miuSiemens)
	tau_K = 80; // ms
	// Initialise defualt parameters
	init();

}




void Neurons::init(){

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
	killer_license = false;
	runaway_killed = false;
	step_killed = -1;
	
	//
	stats_record = false;
	spike_freq_adpt = false;
	
	// perturbation
	step_perturb = -1;
	spike_removed = -1;


}

void Neurons::start_stats_record(){
	stats_record = true;
	
	V_mean.reserve(step_tot);
	V_std.reserve(step_tot);
	
	I_input_mean.reserve(step_tot);
	I_input_std.reserve(step_tot);
	
	I_AMPA_acc.assign(N, 0.0);
	I_AMPA_time_avg.assign(N, 0.0);
	I_NMDA_acc.assign(N, 0.0);
	I_NMDA_time_avg.assign(N, 0.0);
	I_GABA_acc.assign(N, 0.0);
	I_GABA_time_avg.assign(N, 0.0);
	
	IE_ratio.assign(N, 0.0);
}


void Neurons::set_para(string para_str){
	if (!para_str.empty()){
		istringstream para(para_str);
		string para_name, para_value_str, line_str; 
		double para_value;
		while (getline(para, line_str)){
			istringstream line_ss(line_str);
			getline(line_ss, para_name, delim); // get parameter name
			getline(line_ss, para_value_str, delim); // get parameter value (assume double)
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

string Neurons::dump_para(){
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
	dump << "seed" << delim << my_seed << delim << endl;
	return dump.str();
}


void Neurons::random_V(double p){
	// Generate random initial condition for V.
	// Generate uniform random distribution
	if (p < 1.0){
		gen.seed(my_seed);// reseed random engine!
		uniform_real_distribution<double> uniform_dis(0.0, 1.0);
		auto ZeroOne = bind(uniform_dis,gen);

		for (int i = 0; i < N; ++i) {
			// Generate random number.
			//ZeroOne = uniform_dis(gen); // range is [0,1]
			V[i] = V_rt + (V_th - V_rt) / (1.0 - p)*ZeroOne();
		}
	}
	else {cout << "initial firing rate cannot be 100%!" << endl;}
}


void Neurons::update_spikes(int step_current){
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


void Neurons::update_V(int step_current){
	// This function updates menbrane potentials for non-refractory neurons

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
	output_sampled_data_real_time(step_current);
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
	
}


void Neurons::add_sampling_real_time(vector<int> sample_neurons_input, vector<bool> sample_type_input, vector<bool> sample_time_points_input, string samp_file_name_input){
	sample_neurons = sample_neurons_input;
	sample_type = sample_type_input;
	sample_time_points = sample_time_points_input;
	
	samp_file_name = samp_file_name_input.append(to_string(pop_ind)).append(".ygout_samp");
	samp_file.open(samp_file_name);
			
	// initialise
	int sample_time_points_tot = 0;// count non zero elements in sample_time_points
	for (unsigned int i = 0; i < sample_time_points.size(); ++i){
		if (sample_time_points[i]){
			sample_time_points_tot += 1;
		}
	}
	
}



void Neurons::add_sampling(vector<int> sample_neurons_input, vector<bool> sample_type_input, vector<bool> sample_time_points_input){
	sample_neurons = sample_neurons_input;
	sample_type = sample_type_input;
	sample_time_points = sample_time_points_input;
	
	
	// initialise
	int sample_time_points_tot = 0;// count non zero elements in sample_time_points
	for (unsigned int i = 0; i < sample_time_points.size(); ++i){
		if (sample_time_points[i]){
			sample_time_points_tot += 1;
		}
	}
	int sample_neurons_tot = sample_neurons.size();// count non zero elements in sample_time_points
	int sample_type_tot = sample_type.size(); // 8 different data types
	
	
	sample.resize(sample_type_tot); 
	for (int c = 0; c < sample_type_tot; ++c){
		if (sample_type[c]){
			sample[c].resize(sample_neurons_tot);
			for (int i = 0; i < sample_neurons_tot; ++i){
				sample[c][i].reserve(sample_time_points_tot); // reserve and push_back so that it won't be affected by adapting step_tot
			}
		}
	}
}



void Neurons::sample_data(int step_current){
	
	if (!sample_neurons.empty()){
		if (sample_time_points[step_current]){ // push_back is amazing
			for (unsigned int i = 0; i < sample_neurons.size(); ++i){ // performance issue when sampling many neurons?
				int ind_temp = sample_neurons[i];
				if (sample_type[0]){sample[0][i].push_back( V[ind_temp] );}
				if (sample_type[1]){sample[1][i].push_back( I_leak[ind_temp] );}
				if (sample_type[2]){sample[2][i].push_back( I_AMPA[ind_temp] );}
				if (sample_type[3]){sample[3][i].push_back( I_GABA[ind_temp] );}
				if (sample_type[4]){sample[4][i].push_back( I_NMDA[ind_temp] );}
				if (sample_type[5]){sample[5][i].push_back( I_GJ[ind_temp] );}
				if (sample_type[6]){sample[6][i].push_back( I_ext[ind_temp] );}
				if (sample_type[7]){sample[7][i].push_back( I_K[ind_temp] );}
			}
		}
	}

}


void Neurons::set_gaussian_I_ext(vector<double> mean, vector<double> std){
	I_ext_mean = mean;
	I_ext_std = std;
	
	double max_std = *max_element(I_ext_std.begin(), I_ext_std.end());
	if (max_std == 0.0){
		I_ext_std.resize(0);
	}
}


void Neurons::add_perturbation(int step_perturb_input){
	step_perturb = step_perturb_input;
}

void Neurons::add_spike_freq_adpt(){
	spike_freq_adpt = true;
	exp_K_step = exp( -dt / tau_K );
	g_K.assign(N, 0.0);
}


void Neurons::init_runaway_killer(double min_ms, double Hz, double Hz_ms){

	min_pop_size = 100;
	if (N > min_pop_size){
		killer_license = true;
		min_steps = int(round(min_ms / dt));
		runaway_Hz = Hz;
		Hz_steps = int(round(Hz_ms / dt));
	}
}

void Neurons::runaway_check(int step_current){
	if (killer_license == true && runaway_killed == false && step_current > min_steps && step_current > Hz_steps){
		// find mean value of num_ref over the last runaway_steps
		vector<int>::const_iterator first, last;
		// first element to be accumulated
		first = num_spikes_pop.begin() + (step_current - Hz_steps + 1); 
		// one element pass the last element to be accumulated
		last = num_spikes_pop.begin() + (step_current + 1); 
		double mean_Hz = accumulate(first, last, 0.0) / (Hz_steps * dt * 0.001 * N); // 0.001 for converting from ms to sec.
		//be careful!! accumulate range is : [first,last)
		if (mean_Hz >= runaway_Hz){
			runaway_killed = true;
			step_killed = step_current;
			cout << "warning: runaway killed at " << step_current*dt << " (ms) in population" << pop_ind << flush;
			cout << "\t with firing rate at " << mean_Hz << " Hz."<< flush;
		}
	}
}

void Neurons::output_sampled_data_real_time(int step_current){
	
	if (!sample_neurons.empty() && step_current == 0){
		int sample_step_number = 0;
		for (int tt = 0; tt < step_tot; ++tt){
			if (sample_time_points[tt]){sample_step_number += 1;}
		}
		samp_file << indicator << " POPD006" << endl;
		samp_file << pop_ind << delim << int(sample_neurons.size())*sample_step_number << delim << endl;
  		vector< string > data_types = { "V", "I_leak", "I_AMPA", "I_GABA", "I_NMDA", "I_GJ", "I_ext", "I_K" };
		for (unsigned int tt = 0; tt < data_types.size(); ++tt){
			if (sample_type[tt]){samp_file << data_types[tt] << delim; }
		}
		samp_file << endl;
	}
	
	if (!sample_neurons.empty()){
		if (sample_time_points[step_current]){ // push_back is amazing
			for (unsigned int i = 0; i < sample_neurons.size(); ++i){ // performance issue when sampling many neurons?
				int ind_temp = sample_neurons[i];
				if (sample_type[0]){samp_file << V[ind_temp] << delim;}
				if (sample_type[1]){samp_file << I_leak[ind_temp] << delim;}
				if (sample_type[2]){samp_file << I_AMPA[ind_temp] << delim;}
				if (sample_type[3]){samp_file << I_GABA[ind_temp] << delim;}
				if (sample_type[4]){samp_file << I_NMDA[ind_temp] << delim;}
				if (sample_type[5]){samp_file << I_GJ[ind_temp] << delim;}
				if (sample_type[6]){samp_file << I_ext[ind_temp] << delim;}
				if (sample_type[7]){samp_file << I_K[ind_temp] << delim;}
				samp_file << endl;
			}
		}
	}
	
	
}

void Neurons::output_results(ofstream& output_file){

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
	if (stats_record){
		output_file << indicator << " POPD003" << endl;
		output_file << pop_ind << delim << endl;
		write2file(output_file, V_mean);
		write2file(output_file, V_std);
		write2file(output_file, I_input_mean);
		write2file(output_file, I_input_std);
	}
	
	// POPD005 # E-I ratio for each neuron
	if (stats_record){
		output_file << indicator << " POPD005" << endl;
		output_file << pop_ind << delim << endl;
		write2file(output_file, IE_ratio);
	}
	
	// SAMF001 # sampled data file name
	if (!samp_file_name.empty()){
		output_file << indicator << " SAMF001" << endl;
		output_file << samp_file_name << endl;
	}

	/* // This following output protocol is deprecated due to poor memeory performance
	// POPD004 # sampled neuron data
	if (!sample_neurons.empty()){
		output_file << indicator << " POPD004" << endl;
		output_file << pop_ind << delim << sample_neurons.size() << delim << endl;

  		vector< string > data_types = { "V", "I_leak", "I_AMPA", "I_GABA", "I_NMDA", "I_GJ", "I_ext", "I_K" };
		for (unsigned int tt = 0; tt < data_types.size(); ++tt){
			if (sample_type[tt] == true){output_file << data_types[tt] << delim;}
		}
		output_file << endl;


		for (unsigned int c = 0; c < sample_type.size(); ++c){
			if (!sample[c].empty()){
				Neurons::write2file(output_file, sample[c]); // 2D matrix
			}
		}
	}
	*/

}


// Use function templates when you want to perform the same action on types that can be different.
// Use function overloading when you want to apply different operations depending on the type.
// In this case, just save yourself the trouble and use overloading.
void Neurons::write2file(ofstream& output_file,  vector< vector<int> >& v){
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



void Neurons::write2file(ofstream& output_file, vector< vector<double> >& v){
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


void Neurons::write2file(ofstream& output_file, vector<int>& v){
	if (!v.empty()){
		//for (int f : v){ output_file << f << delim; } // range-based "for" in C++11
		for (unsigned int i = 0; i < v.size(); ++i){
			output_file << v[i] << delim;
		}
		output_file << endl;
	}
	else {output_file << " " << endl;}
}



void Neurons::write2file(ofstream& output_file, vector<double>& v){
	if (!v.empty()){
		//for (int f : v){ output_file << f << delim; } // range-based "for" in C++11
		for (unsigned int i = 0; i < v.size(); ++i){
			output_file << v[i] << delim;
		}
		output_file << endl;
	}
	else {output_file << " " << endl;}
}




void Neurons::record_stats(int step_current){
	if (stats_record){
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
		V_mean.push_back(mean_tmp_V);
		V_std.push_back(std_tmp_V);
		I_input_mean.push_back(mean_tmp_I);
		I_input_std.push_back(std_tmp_I);
		
		// accumulate
		//for (unsigned int i = 0; i < N; ++i){ // this manual loop is slow, use transform()
		//	I_input_acc[i] += I_input[i];
		//}
		transform( I_AMPA_acc.begin(), I_AMPA_acc.end(), I_AMPA.begin(), I_AMPA_acc.begin(), plus<double>() );
		transform( I_NMDA_acc.begin(), I_NMDA_acc.end(), I_NMDA.begin(), I_NMDA_acc.begin(), plus<double>() );
		transform( I_GABA_acc.begin(), I_GABA_acc.end(), I_GABA.begin(), I_GABA_acc.begin(), plus<double>() );
		// get time average for each neuron
		if (step_current == step_tot - 1){ // at the end of the last time step
			for (int i = 0; i < N; ++i){
				I_AMPA_time_avg[i] = I_AMPA_acc[i] / step_tot;
				I_NMDA_time_avg[i] = I_NMDA_acc[i] / step_tot;
				I_GABA_time_avg[i] = I_GABA_acc[i] / step_tot;
				// be careful here, for IE_ratio, I_ext is assumed to be always excitatory and I_GJ is not considered
				// also, the only source of I_ext is generated internally
				IE_ratio[i] = I_GABA_time_avg[i] / (I_AMPA_time_avg[i] + I_NMDA_time_avg[i] + I_ext_mean[i]);
			}
		}
	}
}


