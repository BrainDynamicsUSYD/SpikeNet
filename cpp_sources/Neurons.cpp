#include <stdlib.h>
#include <iomanip>
#include <random> //<boost/random.hpp>
#include <cmath>
#include <algorithm>
#include <stdio.h> // for printf
#include <time.h>       /* time */

#include <fstream> // fstream : Stream class to both read and write from / to files
#include <sstream>  // stringstream is input and output


#include "Neurons.h"
#include <functional> // for bind() 
// no need to include what have been included in the header file

Neurons::Neurons(int pop_ind_input, int N_input, double dt_input, int step_tot_input){

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
	// Leak conductance
	g_lk = 0.0167;   // (uS=miuSiemens), time constants=Cm/gL=15 ms!

	// Initialise defualt parameters
	init();

}




void Neurons::init(){
	
	// I_ext
	I_ext_mean = 0.0;
	I_ext_std = 0.0;


	// Initialise arrarys storing instantaneous neuron states
	V.assign(N, V_lk); // All zeros	
	I_leak.assign(N, 0.0);
	I_AMPA.assign(N, 0.0);
	I_GABA.assign(N, 0.0);
	I_NMDA.assign(N, 0.0);
	I_GJ.assign(N, 0.0);
	I_ext.assign(N, 0.0);
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
}



void Neurons::set_para(string para_str, char delim){
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


string Neurons::dump_para(char delim){
	stringstream dump;
	dump << "Cm" << delim << Cm << delim << endl;
	dump << "tau_ref" << delim << tau_ref << delim << endl;
	dump << "V_rt" << delim << V_rt << delim << endl;
	dump << "V_lk" << delim << V_lk << delim << endl;
	dump << "V_th" << delim << V_th << delim << endl;
	dump << "g_lk" << delim << g_lk << delim << endl;
	dump << "seed" << delim << my_seed << delim << endl;
	return dump.str();
}




void Neurons::random_V(double p){
	// Generate random initial condition for V.
	// Generate uniform random distribution
	if ( p < 1.0){
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
	if (I_ext_mean != 0.0){
		if (I_ext_std != 0.0){
			// Gaussian random generator
			gen.seed(my_seed+step_current);// reseed random engine!
			normal_distribution<double> nrm_dist(I_ext_mean, I_ext_std);
			auto gaus = bind(nrm_dist,gen);
			// Generate Gaussian white noise. White means not temporally correlated	
			for (int i = 0; i < N; ++i) { I_ext[i] = gaus(); }				
		}
		else { for (int i = 0; i < N; ++i) { I_ext[i] = I_ext_mean;} }
 	}


	// Collect Currents from all pre-synapses!!!!!!!!!!!!!!!



	// Data sampling
	sample_data(step_current); // this must be here
	// we want to sample the V[] before it's updated!


	// update menbrane potentials
	double Vdot;
	for (int i = 0; i < N; ++i){
		if (ref_step_left[i] == 0){ // Only update the non-refractory neurons
			// leaky current
			I_leak[i] = -g_lk * (V[i] - V_lk); 
			// using simple Euler method
			Vdot = (I_leak[i] + I_AMPA[i] + I_GABA[i] + I_NMDA[i] + I_GJ[i] + I_ext[i])/Cm;
			V[i] += Vdot * dt;
			// Note that delta-function coupling is very different from the above conductance-based model!
		}
	}

}




void Neurons::add_neuron_sampling(vector<int> neuron_sample_ind_input, vector<bool> neuron_sample_type_input){
	neuron_sample_ind = neuron_sample_ind_input;
	neuron_sample_type = neuron_sample_type_input;
	// initialise neuron_sample
	int sample_size = neuron_sample_ind.size();
	int data_type_tot = neuron_sample_type.size(); // 7 different data types
	neuron_sample.resize(data_type_tot); 
	for (int c = 0; c < data_type_tot; ++c){
		if (neuron_sample_type[c]){
			neuron_sample[c].resize(sample_size);
			for (int i = 0; i < sample_size; ++i){ neuron_sample[c][i].reserve(step_tot); } // reserve and push_back so that it won't be affected by adapting step_tot
		}
	}
}

void Neurons::add_pop_sampling(vector<bool> pop_sample_ind_input, vector<bool> pop_sample_type_input){
	pop_sample_ind = pop_sample_ind_input;
	pop_sample_type = pop_sample_type_input;
	pop_sample_size = 0; // accumulator

	// initialise pop_sample
	int max_pop_sample_size = 0;// count non zero elements
	for (unsigned int i = 0; i < pop_sample_ind.size(); ++i){
		if (pop_sample_ind[i]){
			max_pop_sample_size += 1;
		}
	}
	int data_type_tot = pop_sample_type.size(); // 7 different data types
	pop_sample.resize(data_type_tot); 
	for (int c = 0; c < data_type_tot; ++c){
		if (pop_sample_type[c]){
			pop_sample[c].reserve(max_pop_sample_size*N); // reserve and push_back so that it won't be affected by adapting step_tot
		}
	}
}


void Neurons::sample_data(int step_current){
	if (!neuron_sample_ind.empty()){
		for (unsigned int i = 0; i < neuron_sample_ind.size(); ++i){
			int ind_temp = neuron_sample_ind[i];
			if (neuron_sample_type[0]){neuron_sample[0][i].push_back( V[ind_temp] );}
			if (neuron_sample_type[1]){neuron_sample[1][i].push_back( I_leak[ind_temp] );}
			if (neuron_sample_type[2]){neuron_sample[2][i].push_back( I_AMPA[ind_temp] );}
			if (neuron_sample_type[3]){neuron_sample[3][i].push_back( I_GABA[ind_temp] );}
			if (neuron_sample_type[4]){neuron_sample[4][i].push_back( I_NMDA[ind_temp] );}
			if (neuron_sample_type[5]){neuron_sample[5][i].push_back( I_GJ[ind_temp] );}
			if (neuron_sample_type[6]){neuron_sample[6][i].push_back( I_ext[ind_temp] );}
		}
	}

	if (!pop_sample_ind.empty()){
		if (pop_sample_ind[step_current]){
			pop_sample_size += 1;
			if (pop_sample_type[0]){pop_sample[0].push_back(V);}
			if (pop_sample_type[1]){pop_sample[1].push_back(I_leak);}
			if (pop_sample_type[2]){pop_sample[2].push_back(I_AMPA);}
			if (pop_sample_type[3]){pop_sample[3].push_back(I_GABA);}
			if (pop_sample_type[4]){pop_sample[4].push_back(I_NMDA);}
			if (pop_sample_type[5]){pop_sample[5].push_back(I_GJ);}
			if (pop_sample_type[6]){pop_sample[6].push_back(I_ext);}
		}
	}
}



void Neurons::set_gaussian_I_ext(double mean, double std){
	I_ext_mean = mean;
	I_ext_std = std;
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



void Neurons::output_results(ofstream& output_file, char delim, char indicator){

	// POPD001 # spike history of neuron population
	output_file << indicator << " POPD001" << endl;
	output_file << pop_ind << delim << endl;
	write2file(output_file, delim, spike_hist_tot);
	write2file(output_file, delim, num_spikes_pop);
	write2file(output_file, delim, num_ref_pop);

	// POPD002 # neuron parameters in the population
	stringstream dump_count;
	string para_str = dump_para(delim);
	dump_count << para_str;
	string str_temp;
	int var_number = 0;
	while(getline(dump_count,str_temp)){++var_number;} // count number of variables
	output_file << indicator << " POPD002" << endl;
	output_file << pop_ind << delim << var_number << delim << endl;
	output_file << para_str;

	// POPD003 # sampled populational data
	if (!pop_sample_ind.empty()){
		output_file << indicator << " POPD003" << endl;
		output_file << pop_ind << delim << pop_sample_size << delim  <<  endl;

		// data_type: [V,I_leak,I_AMPA,I_GABA,I_NMDA,I_GJ,I_ext]
		if (pop_sample_type[0] == true){output_file << "V" << delim;}
		if (pop_sample_type[1] == true){output_file << "I_leak" << delim;}
		if (pop_sample_type[2] == true){output_file << "I_AMPA" << delim;}
		if (pop_sample_type[3] == true){output_file << "I_GABA" << delim;}
		if (pop_sample_type[4] == true){output_file << "I_NMDA" << delim;}
		if (pop_sample_type[5] == true){output_file << "I_GJ" << delim;}
		if (pop_sample_type[6] == true){output_file << "I_ext" << delim;}
		output_file << endl;
		// any clever way to do the above? (use map?)
		for (unsigned int c = 0; c < pop_sample_type.size(); ++c){
			if (!pop_sample[c].empty()){
				write2file(output_file, delim, pop_sample[c]); // 2D matrix
			}
		}
	}


	// POPD004 # sampled neuron data
	if (!neuron_sample_ind.empty()){
		output_file << indicator << " POPD004" << endl;
		output_file << pop_ind << delim << neuron_sample_ind.size() << delim << endl;

			
		/*/ data_type: [V,I_leak,I_AMPA,I_GABA,I_NMDA,I_GJ,I_ext]
		if (neuron_sample_type[0] == true){output_file << "V" << delim;}
		if (neuron_sample_type[1] == true){output_file << "I_leak" << delim;}
		if (neuron_sample_type[2] == true){output_file << "I_AMPA" << delim;}
		if (neuron_sample_type[3] == true){output_file << "I_GABA" << delim;}
		if (neuron_sample_type[4] == true){output_file << "I_NMDA" << delim;}
		if (neuron_sample_type[5] == true){output_file << "I_GJ" << delim;}
		if (neuron_sample_type[6] == true){output_file << "I_ext" << delim;}
		output_file << endl;
		// any clever way to do the above? (use map?) */


  		vector< string > data_types = { "V", "I_leak", "I_AMPA", "I_GABA", "I_NMDA", "I_GJ", "I_ext" };
		for (unsigned int tt = 0; tt < data_types.size(); ++tt){
			if (neuron_sample_type[tt] == true){output_file << data_types[tt] << delim;}
		}
		output_file << endl;


		for (unsigned int c = 0; c < neuron_sample_type.size(); ++c){
			if (!neuron_sample[c].empty()){
				write2file(output_file, delim, neuron_sample[c]); // 2D matrix
			}
		}
	}



}






// Use function templates when you want to perform the same action on types that can be different.
// Use function overloading when you want to apply different operations depending on the type.
// In this case, just save yourself the trouble and use overloading.
void Neurons::write2file(ofstream& output_file, char delim, vector< vector<int> >& v){
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



void Neurons::write2file(ofstream& output_file, char delim, vector< vector<double> >& v){
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


void Neurons::write2file(ofstream& output_file, char delim, vector<int>& v){
	if (!v.empty()){
		//for (int f : v){ output_file << f << delim; } // range-based "for" in C++11
		for (unsigned int i = 0; i < v.size(); ++i){
			output_file << v[i] << delim;
		}
		output_file << endl;
	}
	else {output_file << " " << endl;}
}




