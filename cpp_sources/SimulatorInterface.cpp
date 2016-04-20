#include "NeuronNetwork.h"
#include "SimulatorInterface.h"

#include <chrono> // #include <boost/chrono.hpp>

using namespace std;

SimulatorInterface::SimulatorInterface(){
	// define default output data file format
	output_suffix = ".ygout"; // filename extension
	delim = ',';
	indicator = '>';
	commentor = '#';
}

bool SimulatorInterface::import(string in_filename_input){
	// read main input file
	in_filename = in_filename_input;
	inputfile.open(in_filename, ifstream::in);
	if (inputfile){
		cout << "------------------------------------------------------------" << endl;
		cout << "Importing " << in_filename << "..." << endl;
	}
	else {
		cout << "Error: cannot open file " << in_filename << endl;
		return 0;
	}
	
	// claim variables to be assigned
	int Num_pop;
	double dt;
	int step_tot;
	vector<int> N_array;
	// temporary data storage
	vector<string> pop_para; // parameters of neuron popluation
	string syn_para; // parameters of synapses
	vector< vector<double> > ext_spike_settings; // (int pop_ind, int type_ext, double K_ext, int Num_ext, double rate_ext, int ia, int ib)
	vector< vector<double> > rate_ext_t; // vector<duoble> rate_ext(t)
	
	vector<int> ext_current_settings_pop_ind; // (int pop_ind)
	vector< vector<double> > ext_current_settings; // (vector<double> mean, vector<double> std)
	
	vector<double> init_condition_settings; // (double p_fire)
	
	vector<int> neuron_sample_pop_ind; //
	vector< vector<int> > neuron_sample_neurons; //
	vector< vector<bool> > neuron_sample_type; //
	vector< vector<bool> > neuron_sample_time_points; // the length of time vector
	
	vector< vector<int> > syn_sample_pop_ind; //
	vector< vector<int> > syn_sample_neurons; //
	vector< vector<bool> > syn_sample_time_points; // the length of time vector
	
	vector<int> neuron_stats_setting;
	vector< vector<int> > syn_stats_setting;
	vector< vector<int> > STD_setting;
	vector< vector<int> > inh_STDP_setting;
	
	string syn_filename; // name of file that defines synaptic connection
	vector< vector<double> > runaway_killer_setting; // [pop_ind, min_ms, runaway_Hz, Hz_ms]
	
	vector< vector<int> > step_perturb_setting; // [pop_ind step_perturb]

	// read data
	string line_str, entry_str; // temporary container for the entire while loop
	size_t found;
	while (getline(inputfile, line_str)){
		istringstream line_ss(line_str); // local variable: lifetime is only a signle loop (scope is within "{}")! reusing it by ".str()" may cause fatal error because the internal error flag remains the same. Everytime using "while(getline(...)){}" will change internal error flag!
		if (line_str.empty()){continue;}
		else if (line_str.front() == commentor){continue;}
		else if (line_str.front() == indicator){// read data-info line

			// read number of neurons in each population
			found = line_str.find("INIT001");
			if (found != string::npos){// if found match
				cout << "\t Reading number of neurons in each population..." << endl;
				read_next_line_as_vector(N_array);

				// Infer number of populations
				Num_pop = N_array.size();
				pop_para.resize(Num_pop); // ??????
				continue; // move to next line
			}


			// read time step length and total steps
			found = line_str.find("INIT002");
			if (found != string::npos){// if found match
				cout << "\t Reading time step length and total steps..." << endl;
				getline(inputfile, line_str);istringstream line_ss(line_str);// Read next line

				dt = read_next_entry<double>(line_ss);
				step_tot = read_next_entry<int>(line_ss);
				continue; // move to next line
			}



			// Random initial distributions for membrane potentials 
			found = line_str.find("INIT003");
			if (found != string::npos){// if found match
				cout << "\t Reading random initial distributions for membrane potentials..." << endl;
				read_next_line_as_vector(init_condition_settings); // p_fire
				continue; // move to next line
			}


			// read external current setting
			found = line_str.find("INIT004");
			if (found != string::npos){// if found match
				cout << "\t Reading external current settings..." << endl;
				getline(inputfile, line_str);istringstream line_ss(line_str);// Read next line 
				ext_current_settings_pop_ind.push_back(read_next_entry<int>(line_ss));
				ext_current_settings.resize(ext_current_settings.size()+1);
				read_next_line_as_vector(ext_current_settings.back()); // vector<double> mean
				ext_current_settings.resize(ext_current_settings.size()+1);
				read_next_line_as_vector(ext_current_settings.back()); // vector<double> std
				continue; // move to next line
			}



			// read external spike setting
			found = line_str.find("INIT005");
			if (found != string::npos){// if found match
				cout << "\t Reading external spike settings..." << endl;
				// [pop_ind,type_ext,K_ext,Num_ext,ia,ib]
				ext_spike_settings.resize(ext_spike_settings.size()+1);
				read_next_line_as_vector(ext_spike_settings.back());
				// [rate_ext_t]
				rate_ext_t.resize(rate_ext_t.size()+1);
				read_next_line_as_vector(rate_ext_t.back());
				continue; // move to next line
			}

			// read perturbation setting
			found = line_str.find("INIT007");
			if (found != string::npos){// if found match
				cout << "\t Reading perturbation settings..." << endl;
				// [pop_ind, step_perturb]
				step_perturb_setting.resize(step_perturb_setting.size()+1);
				read_next_line_as_vector(step_perturb_setting.back());
				continue; // move to next line
			}

			found = line_str.find("INIT008");
			if (found != string::npos){// if found match
				cout << "\t Reading short term depression settings..." << endl;
				// [pop_ind_pre, pop_ind_post, STD_on_step]
				STD_setting.resize(STD_setting.size()+1);
				read_next_line_as_vector(STD_setting.back());
				continue; // move to next line
			}
			
			found = line_str.find("INIT009");
			if (found != string::npos){// if found match
				cout << "\t Reading inhibitory STDP settings..." << endl;
				// [pop_ind_pre, pop_ind_post, inh_STDP_on_step]
				inh_STDP_setting.resize(inh_STDP_setting.size()+1);
				read_next_line_as_vector(inh_STDP_setting.back());
				continue; // move to next line
			}
			
			
			// read neuron population parameter setting
			found = line_str.find("PARA001");
			if (found != string::npos){// if found match
				cout << "\t Reading non-default neuron population parameters..." << endl;
				getline(inputfile, line_str);istringstream line_ss(line_str);// Read next line
				// pop_ind and number of parameters
				int pop_ind = read_next_entry<int>(line_ss); 
				int num_para = read_next_entry<int>(line_ss);
				// store parameter settings
				for (int c = 0; c < num_para; ++c){
					getline(inputfile, line_str); // read next line
					pop_para[pop_ind] += line_str; // append (+=)
				}
				continue; // move to next line
			}



			// read synapse parameter setting
			found = line_str.find("PARA002");
			if (found != string::npos){// if found match
				cout << "\t Reading non-default synapse parameters..." << endl;
				getline(inputfile, line_str);istringstream line_ss(line_str);// Read next line
				// number of parameters
				int num_para = read_next_entry<int>(line_ss);
				// store parameter settings
				for (int c = 0; c < num_para; ++c){
					getline(inputfile, line_str); // read next line
					syn_para += line_str; // append (+=)
				}
				continue; // move to next line
			}



			// read runaway killer setting
			found = line_str.find("KILL001");
			if (found != string::npos){
				cout << "\t Reading runaway killer setting..." << endl;
				// double pop_ind, min_ms, double Hz, double Hz_ms
				runaway_killer_setting.resize(runaway_killer_setting.size()+1);
				read_next_line_as_vector(runaway_killer_setting.back());
				continue;
			}



			// read neuron data sampling setting
			found = line_str.find("SAMP001");
			if (found != string::npos){// if found match
				cout << "\t Reading neuron data sampling settings..." << endl;
				getline(inputfile, line_str);istringstream line_ss(line_str);// Read next line
				// int pop_ind;
				neuron_sample_pop_ind.push_back(read_next_entry<int>(line_ss));
				// sample_type
				neuron_sample_type.resize(neuron_sample_type.size()+1);
				read_next_line_as_vector(neuron_sample_type.back());
				// sample_neurons
				neuron_sample_neurons.resize(neuron_sample_neurons.size()+1);
				read_next_line_as_vector(neuron_sample_neurons.back());
				// sample_t_ind
				neuron_sample_time_points.resize(neuron_sample_time_points.size()+1);
				read_next_line_as_vector(neuron_sample_time_points.back());
				continue; // move to next line
			}

			// read synapase data sampling setting
			found = line_str.find("SAMP002");
			if (found != string::npos){// if found match
				cout << "\t Reading synapse data sampling settings..." << endl;
				// pop indices
				syn_sample_pop_ind.resize(syn_sample_pop_ind.size()+1);
				read_next_line_as_vector(syn_sample_pop_ind.back());
				// sample_neurons
				syn_sample_neurons.resize(syn_sample_neurons.size()+1);
				read_next_line_as_vector(syn_sample_neurons.back());
				// sample_t_ind
				syn_sample_time_points.resize(syn_sample_time_points.size()+1);
				read_next_line_as_vector(syn_sample_time_points.back());
				continue; // move to next line
			}
			
			// neuron_stats_record
			found = line_str.find("SAMP003");
			if (found != string::npos){// if found match
				cout << "\t Reading neuron stats record settings..." << endl;
				getline(inputfile, line_str);istringstream line_ss(line_str);// Read next line
				// int pop_ind;
				neuron_stats_setting.push_back(read_next_entry<int>(line_ss));
				continue; // move to next line
			}
			
			
			// neuron_stats_record
			found = line_str.find("SAMP004");
			if (found != string::npos){// if found match
				cout << "\t Reading synaptic stats record settings..." << endl;
				syn_stats_setting.resize(syn_stats_setting.size()+1);
				read_next_line_as_vector(syn_stats_setting.back());
				continue; // move to next line
			}
			

			// read non-default synapse definition file name
			found = line_str.find("SYNF001");
			if (found != string::npos){
				cout << "\t Reading non-default synapse definition file name..." << endl;
				getline(inputfile, syn_filename); // Read next line
				cout << "\t\t non-default synapse definition: " << syn_filename << endl;
				continue;
			}

			else{
				cout << "Unrecognized command sequence: " << line_str << endl;
			}

		} //read data-info line	

	} // while
	inputfile.close();
	inputfile.clear();

	// read synapse definition input file
	if (syn_filename.empty()){ // default file name
		syn_filename = in_filename;
		syn_filename.append("_syn"); // .ygin_syn
	}
	inputfile.open(syn_filename, ifstream::in); // close, clear and then open: reusing is a bad practice?
	if (inputfile){}
	else {
		cout << "Error: cannot open syn file " << syn_filename << endl;
		return 0;
	}
	
	// temporary data storage
	vector< vector<int> > I_temp, J_temp; //for chemical synapses
	vector< vector<double> > K_temp, D_temp; //for chemical synapses
	vector< vector<int> > IJKD_chem_info; // {Type(0:AMPA/1:GABAa/2:NMDA),Pre_pop_ind(0/1/...),Post_pop_ind(0/1/...)}
	
	// read data
	while (getline(inputfile, line_str)){
		istringstream line_ss(line_str);
		if (line_str.empty()){continue;}
		else if (line_str.front() == commentor){continue;}
		else if (line_str.front() == indicator){// read data-info line
			
			// read connection data: [I;J;K;D] (chemical)
			found = line_str.find("INIT006");
			if (found != string::npos){// if found match
				cout << "\t Reading chemical connection..." << endl;
				// [type, pre_pop_ind, post_pop_ind]
				IJKD_chem_info.resize(IJKD_chem_info.size()+1);
				read_next_line_as_vector(IJKD_chem_info.back());
				// I_temp
				I_temp.resize(I_temp.size()+1);
				read_next_line_as_vector(I_temp.back());
				// J_temp
				J_temp.resize(J_temp.size()+1);
				read_next_line_as_vector(J_temp.back());
				// K_temp
				K_temp.resize(K_temp.size()+1);
				read_next_line_as_vector(K_temp.back());
				// D_temp
				D_temp.resize(D_temp.size()+1);		
				read_next_line_as_vector(D_temp.back());

				continue; // move to next line
			}
		}//read data-info line	

		
	} // while

	inputfile.close();
	inputfile.clear();



	// build NeuronNetwork based on data imported
	network = NeuronNetwork(N_array, dt, step_tot);
	cout << "\t Network created." << endl;
	cout << "\t Initialising neuron populations...";
	for (unsigned int ind = 0; ind < N_array.size(); ++ind){
		network.NeuronPopArray.push_back(Neurons(ind, N_array[ind], network.dt, network.step_tot));
		network.NeuronPopArray.back().set_para(pop_para[ind], delim);
		cout << ind+1 << "...";
	}
	cout << "done." << endl;

	// chemical connections
	if (I_temp.size() != 0){
		cout << "\t Initialising chemical synapses... ";
		for (unsigned int ind = 0; ind < I_temp.size(); ++ind){
			int type = IJKD_chem_info[ind][0];
			int i_pre = IJKD_chem_info[ind][1];
			int j_post = IJKD_chem_info[ind][2];
			network.ChemicalSynapsesArray.push_back(ChemicalSynapses(network.dt, network.step_tot));
			network.ChemicalSynapsesArray.back().init(type, i_pre, j_post, network.N_array[i_pre], network.N_array[j_post], I_temp[ind], J_temp[ind], K_temp[ind], D_temp[ind]);
			network.ChemicalSynapsesArray.back().set_para(syn_para, delim);
			cout << ind+1 << "...";
		}
		cout << "done." << endl;
	}

	// external spike setting (int pop_ind, int type_ext, double K_ext, int Num_ext, 
	// vector<double> rate_ext_t, int ia, int ib)
	if (ext_spike_settings.size() != 0){
		cout << "\t External spike settings...";
		for (unsigned int ind = 0; ind < ext_spike_settings.size(); ++ind){
			int j_post = int(ext_spike_settings[ind][0]);
			int type_ext = int(ext_spike_settings[ind][1]);
			double K_ext = ext_spike_settings[ind][2];
			int Num_ext = int(ext_spike_settings[ind][3]);
			int ia = int(ext_spike_settings[ind][4]);
			int ib = int(ext_spike_settings[ind][5]);
			network.ChemicalSynapsesArray.push_back(ChemicalSynapses(network.dt, network.step_tot));
			network.ChemicalSynapsesArray.back().init(type_ext, j_post, network.N_array[j_post], K_ext, Num_ext, rate_ext_t[ind], ia, ib);
			network.ChemicalSynapsesArray.back().set_para(syn_para, delim);
			cout << ind+1 << "...";
		}
		cout << "done." << endl;
	}

	// external current setting
	if (ext_current_settings_pop_ind.size() != 0){
		cout << "\t External current settings...";
		for (unsigned int ind = 0; ind < ext_current_settings_pop_ind.size(); ++ind){
			int pop_ind = ext_current_settings_pop_ind[ind];
			vector<double> mean = ext_current_settings[ ind*2 ];
			vector<double> std = ext_current_settings[ ind*2+1 ];
			network.NeuronPopArray[pop_ind].set_gaussian_I_ext(mean, std);
			cout << ind+1 << "...";
		}
		cout << "done." << endl;
	}

	// random initial condition settings (double p_fire)
	if (init_condition_settings.size() != 0){
		cout << "\t Random initial condition settings...";
		for (unsigned int pop_ind = 0; pop_ind < init_condition_settings.size(); ++pop_ind){
			double p_fire = init_condition_settings[pop_ind];
			network.NeuronPopArray[pop_ind].random_V(p_fire);
			cout << pop_ind+1 << "...";
		}
		cout << "done." << endl;
	}
	
	// perturbation setting
	if (step_perturb_setting.size() != 0){
		cout << "\t Perturbation settings...";
		for (unsigned int i = 0; i < step_perturb_setting.size(); ++i){
			int pop_ind = step_perturb_setting[i][0];
			int step_perturb = step_perturb_setting[i][1];
			network.NeuronPopArray[pop_ind].add_perturbation(step_perturb);
			cout << pop_ind+1 << "...";
		}
		cout << "done." << endl;
	}
	
	
	
	if (STD_setting.size() != 0){
		cout << "\t Short-term depression settings...";
		for (unsigned int i = 0; i < STD_setting.size(); ++i){
			int pop_ind_pre = STD_setting[i][0];
			int pop_ind_post = STD_setting[i][1];
			int STD_on_step = STD_setting[i][2];

			int syn_match = false;
			for (unsigned int s = 0; s < network.ChemicalSynapsesArray.size(); ++s){
				if (network.ChemicalSynapsesArray[s].pop_ind_pre == pop_ind_pre &&
				network.ChemicalSynapsesArray[s].pop_ind_post == pop_ind_post){ // find the right synapse object
					network.ChemicalSynapsesArray[s].add_short_term_depression(STD_on_step);
					cout << i+1 << "...";
					syn_match = true;
				}
			}
			if (!syn_match){
				cout << "(no match found!)...";
			}
		}
		cout << "done." << endl;
	}
	
	if (inh_STDP_setting.size() != 0){
		cout << "\t Inhibitory STDP settings...";
		
		for (unsigned int i = 0; i < inh_STDP_setting.size(); ++i){
			int pop_ind_pre = inh_STDP_setting[i][0];
			int pop_ind_post = inh_STDP_setting[i][1];
			int inh_STDP_on_step = inh_STDP_setting[i][2];

			int syn_match = false;
			for (unsigned int s = 0; s < network.ChemicalSynapsesArray.size(); ++s){
				if (network.ChemicalSynapsesArray[s].pop_ind_pre == pop_ind_pre &&
				network.ChemicalSynapsesArray[s].pop_ind_post == pop_ind_post){ // find the right synapse object
					network.ChemicalSynapsesArray[s].add_inh_STDP(inh_STDP_on_step);
					cout << i+1 << "...";
					syn_match = true;
				}
			}
			if (!syn_match){
				cout << "(no match found!)...";
			}
		}
		cout << "done." << endl;
	}
	
	// neuron data sampling settings
	if (neuron_sample_pop_ind.size() != 0){
		cout << "\t Neuron data sampling settings...";
		for (unsigned int ind = 0; ind < neuron_sample_pop_ind.size(); ++ind){
			int pop_ind = neuron_sample_pop_ind[ind];
			network.NeuronPopArray[pop_ind].add_sampling(neuron_sample_neurons[ind], neuron_sample_type[ind], neuron_sample_time_points[ind]);
			cout << ind+1 << "...";
		}
		cout << "done." << endl;
	}	
	
	// neuron V mean std record settings
	if (neuron_stats_setting.size() != 0){
		cout << "\t Neuron stats record settings...";
		for (unsigned int ind = 0; ind < neuron_stats_setting.size(); ++ind){
			int pop_ind = neuron_stats_setting[ind];
			network.NeuronPopArray[pop_ind].start_stats_record();
			cout << ind+1 << "...";
		}
		cout << "done." << endl;
	}	
	
	
	
	// syn data sampling settings
	if (syn_sample_pop_ind.size() != 0){
		cout << "\t Synapse data sampling settings...";
		for (unsigned int ind = 0; ind < syn_sample_pop_ind.size(); ++ind){
			
			int pop_ind_pre = syn_sample_pop_ind[ind][0];
			int pop_ind_post = syn_sample_pop_ind[ind][1];
			int syn_type = syn_sample_pop_ind[ind][2];
			
			int syn_sample_match = false;
			for (unsigned int s = 0; s < network.ChemicalSynapsesArray.size(); ++s){
				if (network.ChemicalSynapsesArray[s].pop_ind_pre == pop_ind_pre &&
					network.ChemicalSynapsesArray[s].pop_ind_post == pop_ind_post &&
				network.ChemicalSynapsesArray[s].synapses_type == syn_type){ // find the right synapse object
					network.ChemicalSynapsesArray[s].add_sampling(syn_sample_neurons[ind], syn_sample_time_points[ind]);
					cout << ind+1 << "...";
					syn_sample_match = true;
				}
			}
			if (!syn_sample_match){
				cout << "(no match found!)...";
			}
		}
		cout << "done." << endl;
	}	



	// syn I mean std record settings
	if (syn_stats_setting.size() != 0){
		cout << "\t Synapse stats record settings...";
		for (unsigned int ind = 0; ind < syn_stats_setting.size(); ++ind){
			
			int pop_ind_pre = syn_stats_setting[ind][0];
			int pop_ind_post = syn_stats_setting[ind][1];
			int syn_type = syn_stats_setting[ind][2];
			
			int syn_I_match = false;
			for (unsigned int s = 0; s < network.ChemicalSynapsesArray.size(); ++s){
				if (network.ChemicalSynapsesArray[s].pop_ind_pre == pop_ind_pre &&
					network.ChemicalSynapsesArray[s].pop_ind_post == pop_ind_post &&
				network.ChemicalSynapsesArray[s].synapses_type == syn_type){ // find the right synapse object
					network.ChemicalSynapsesArray[s].start_stats_record();
					cout << ind+1 << "...";
					syn_I_match = true;
				}
			}
			if (!syn_I_match){
				cout << "(no match found!)...";
			}
		}
		cout << "done." << endl;
	}	




	// initialise runaway-killer
	if (runaway_killer_setting.size() != 0){
		cout << "\t Runaway killer licensing for pop...";
		for (unsigned int ind = 0; ind < runaway_killer_setting.size(); ++ind){
			int pop_ind = int(runaway_killer_setting[ind][0]);
			network.NeuronPopArray[pop_ind].init_runaway_killer(runaway_killer_setting[ind][1], runaway_killer_setting[ind][2], runaway_killer_setting[ind][3]);
			cout << pop_ind+1 << "...";
		}
		cout << "done." << endl;
		cout << "\t \t No women, no kids." << endl;
	}


	cout << "Importing done." << endl;
	return 1;
}

void SimulatorInterface::simulate(){
	// simulate
	for (int i = 0; i < network.step_tot; ++i){
		network.update(i);
	}
	cout << "Simulation done." << endl;
}

void SimulatorInterface::output_results(){


	// creat output file
	out_filename = gen_out_filename();
	output_file.open(out_filename);
	
	
	// write data
	cout << "Outputting results into file..." << endl;
	
	// KILL002 # step at which runaway activity is killed
	output_file << indicator << " KILL002" << endl;
	output_file << network.step_killed << delim << endl;

	// dump population data
	for (int i = 0; i < network.Num_pop; i++){
		network.NeuronPopArray[i].output_results(output_file, delim, indicator);
	}

	// dump synapse data
	for (unsigned int i = 0; i < network.ChemicalSynapsesArray.size(); i++){
		network.ChemicalSynapsesArray[i].output_results(output_file, delim, indicator);
	}

	// attach input file (.ygin) to the output file for data completeness
 	ifstream in_file_attach( in_filename ) ;
        output_file << in_file_attach.rdbuf() ;


	cout << "Outputting done." << endl << "------------------------------------------------------------" << endl;
	// Write data file name to stdout and use "grep ygout" to extract it!
	cout << "Data file name is: " << endl;
	cout << "	" << out_filename << endl;
}



string SimulatorInterface::gen_out_filename(){
	// creat output file name using some pre-defined format
	ostringstream convert_temp;   // stream used for the conversion
	// Make use of the input file path and name
	string in_filename_trim;
	istringstream in_filename_ss(in_filename);
	getline(in_filename_ss, in_filename_trim, '.');
	// Using time since epoch to stamp the output file
	chrono::high_resolution_clock::duration tse = chrono::high_resolution_clock::now().time_since_epoch();
	chrono::milliseconds tse_ms = chrono::duration_cast<chrono::milliseconds>(tse);
	unsigned long long time_stamp = tse_ms.count(); // ms from epoch, better than "time_t time_stamp = time(0);"
	// Combine all the parts of the output file path and name
	convert_temp << in_filename_trim << "_" << time_stamp << output_suffix; // insert the textual representation of 'Number' in the characters in the stream
	return convert_temp.str(); // set 'Result' to the contents of the stream
}




template < typename Type > Type SimulatorInterface::read_next_entry(istringstream &line_ss){
	string entry_str;
	Type entry;
	getline(line_ss, entry_str, delim);
	if (entry_str.empty()){
		cout << "ERROR: SimulatorInterface::read_next_entry: empty content!" << endl;
	}
	stringstream(entry_str) >> entry;
	return entry;
}

template < typename Type, typename A > void SimulatorInterface::read_next_line_as_vector(vector<Type, A> &vec){
	
	string line_str, entry_str;
	Type entry;
	getline(inputfile, line_str); istringstream line_ss(line_str); // read next line
	if (line_str.empty()){
		cout << "ERROR: SimulatorInterface::read_next_line_as_vector: empty string!" << endl;
	}
	while (getline(line_ss, entry_str, delim)){
		stringstream(entry_str) >> entry;
		vec.push_back(entry);
	}

}

void SimulatorInterface::write2file(vector<int>& v){
	if (!v.empty()){
		//for (int f : v){ output_file << f << delim; } // range-based "for" in C++11
		for (unsigned int i = 0; i < v.size(); ++i){
			output_file << v[i] << delim;
		}
		output_file << endl;
	}
	else {output_file << " " << endl;}
}




