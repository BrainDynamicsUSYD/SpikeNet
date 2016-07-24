//#include "NeuroNet.h"
#include "SimuInterface.h"
#include <chrono> // #include <boost/chrono.hpp>
#include <ctime>

using namespace std;

SimuInterface::SimuInterface(){
	// define default output data file format
	output_suffix = ".ygout"; // filename extension
	delim = ',';
	indicator = '>';
	commentor = '#';
	
}

bool SimuInterface::import(string in_filename_input){
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
	vector< vector<double> > ext_spike_settings; // (int pop_ind, int type_ext, double K_ext, int Num_ext, double rate_ext)
	vector< vector<double> > rate_ext_t; // vector<duoble> rate_ext(t)
	vector< vector<bool> > ext_spike_neurons; //
	
	vector<int> ext_current_settings_pop_ind; // (int pop_ind)
	vector< vector<double> > ext_current_settings; // (vector<double> mean, vector<double> std)
	vector<int> ext_conductance_settings_pop_ind; // (int pop_ind)
	vector< vector<double> > ext_conductance_settings; // (vector<double> mean, vector<double> std)
	
	vector<double> init_condition_settings_depre; // (double p_fire)
	vector< vector<double> >  init_condition_settings; // (double r_V0, double p_fire)
	
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
	vector<int> spike_freq_adpt_setting; //
	
	vector<int> LFP_record_pop_ind; //
	vector< vector < vector<bool> > > LFP_record_setting; //
	
	string syn_filename; // name of file that defines synaptic connection
	vector< vector<double> > runaway_killer_setting; // [pop_ind, min_ms, runaway_Hz, Hz_ms]
	
	vector< vector<int> > step_perturb_setting; // [pop_ind step_perturb]

	int synapse_model_choice = 0;
	
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



			// Random initial distributions for membrane potentials (deprecated)
			found = line_str.find("INIT003");
			if (found != string::npos){// if found match
				cout << "\t Reading random initial distributions for membrane potentials..." << endl;
				read_next_line_as_vector(init_condition_settings_depre); // p_fire
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
				// [pop_ind,type_ext,K_ext,Num_ext]
				ext_spike_settings.resize(ext_spike_settings.size()+1);
				read_next_line_as_vector(ext_spike_settings.back());
				//	ext_spike_neurons
				ext_spike_neurons.resize(ext_spike_neurons.size()+1);
				read_next_line_as_vector(ext_spike_neurons.back());
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
			
			found = line_str.find("INIT010");
			if (found != string::npos){// if found match
				cout << "\t Reading spike-frequency adaptation settings..." << endl;
				// [pop_ind]
				getline(inputfile, line_str);istringstream line_ss(line_str);// Read next line 
				spike_freq_adpt_setting.push_back(read_next_entry<int>(line_ss));			
				continue; // move to next line
			}
			
			// Random initial distributions for membrane potentials and initial firing probabilities
			found = line_str.find("INIT011");
			if (found != string::npos){// if found match
				cout << "\t Reading random initial distributions for membrane potentials and initial firing probabilities..." << endl;
				init_condition_settings.resize(2);
				read_next_line_as_vector(init_condition_settings[0]); // r_V0
				read_next_line_as_vector(init_condition_settings[1]); // p_fire
				continue; // move to next line
			}
			
			// read external conductance setting
			found = line_str.find("INIT012");
			if (found != string::npos){// if found match
				cout << "\t Reading external conductance settings..." << endl;
				getline(inputfile, line_str);istringstream line_ss(line_str);// Read next line 
				ext_conductance_settings_pop_ind.push_back(read_next_entry<int>(line_ss));
				ext_conductance_settings.resize(ext_conductance_settings.size()+1);
				read_next_line_as_vector(ext_conductance_settings.back()); // vector<double> mean
				ext_conductance_settings.resize(ext_conductance_settings.size()+1);
				read_next_line_as_vector(ext_conductance_settings.back()); // vector<double> std
				continue; // move to next line
			}
			
			// read synapse model choice
			found = line_str.find("INIT013");
			if (found != string::npos){// if found match
				cout << "\t Reading synapse model choice..." << endl;
				getline(inputfile, line_str);istringstream line_ss(line_str);// Read next line
				synapse_model_choice = read_next_entry<int>(line_ss); 
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
			
			// LFP record
			found = line_str.find("SAMP005");
			if (found != string::npos){// if found match
				cout << "\t Reading LFP record settings..." << endl;
				getline(inputfile, line_str);istringstream line_ss(line_str);// Read next line
				// int pop_ind, sample_number;
				LFP_record_pop_ind.push_back(read_next_entry<int>(line_ss));
				int n = read_next_entry<int>(line_ss);
				LFP_record_setting.resize(LFP_record_setting.size()+1);
				for (int ind = 0; ind < n; ++ind){
					LFP_record_setting.back().resize(LFP_record_setting.back().size()+1);
					read_next_line_as_vector(LFP_record_setting.back().back());
				}
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

	out_filename = gen_out_filename();


	// build NeuroNet based on data imported
	network = NeuroNet(N_array, dt, step_tot, delim, indicator);
	cout << "\t Network created." << endl;
	cout << "\t Initialising neuron populations...";
	for (unsigned int ind = 0; ind < N_array.size(); ++ind){
		network.NeuroPopArray.push_back(new NeuroPop(ind, N_array[ind], network.dt, network.step_tot, delim, indicator));
		network.NeuroPopArray.back()->set_para(pop_para[ind]);
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
			network.ChemSynArray.push_back(new ChemSyn(network.dt, network.step_tot, delim, indicator));
			network.ChemSynArray.back()->init(type, i_pre, j_post, network.N_array[i_pre], network.N_array[j_post], I_temp[ind], J_temp[ind], K_temp[ind], D_temp[ind]);
			network.ChemSynArray.back()->set_para(syn_para);
			if (synapse_model_choice != 0){
				network.ChemSynArray.back()->set_synapse_model(synapse_model_choice);
			}
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
			network.ChemSynArray.push_back(new ChemSyn(network.dt, network.step_tot, delim, indicator));
			network.ChemSynArray.back()->init(type_ext, j_post, network.N_array[j_post], K_ext, Num_ext, rate_ext_t[ind], ext_spike_neurons[ind]);
			network.ChemSynArray.back()->set_para(syn_para);
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
			network.NeuroPopArray[pop_ind]->set_gaussian_I_ext(mean, std);
			cout << ind+1 << "...";
		}
		cout << "done." << endl;
	}

	// external conductance setting
	if (ext_conductance_settings_pop_ind.size() != 0){
		cout << "\t External conductance settings...";
		for (unsigned int ind = 0; ind < ext_conductance_settings_pop_ind.size(); ++ind){
			int pop_ind = ext_conductance_settings_pop_ind[ind];
			vector<double> mean = ext_conductance_settings[ ind*2 ];
			vector<double> std = ext_conductance_settings[ ind*2+1 ];
			network.NeuroPopArray[pop_ind]->set_gaussian_g_ext(mean, std);
			cout << ind+1 << "...";
		}
		cout << "done." << endl;
	}
	
	// random initial condition settings (double p_fire)
	if (init_condition_settings_depre.size() != 0){
		cout << "\t Random initial condition settings...";
		for (unsigned int pop_ind = 0; pop_ind < init_condition_settings_depre.size(); ++pop_ind){
			double p_fire = init_condition_settings_depre[pop_ind];
			network.NeuroPopArray[pop_ind]->random_V(p_fire);
			cout << pop_ind+1 << "...";
		}
		cout << "done." << endl;
	}
	
	// random initial condition settings (double r_V0, double p_fire)
	if (init_condition_settings.size() != 0){
		cout << "\t Random initial condition settings...";
		for (unsigned int pop_ind = 0; pop_ind < init_condition_settings.size(); ++pop_ind){
			double r_V0 = init_condition_settings[0][pop_ind];
			double p_fire = init_condition_settings[1][pop_ind];
			network.NeuroPopArray[pop_ind]->set_init_condition(r_V0, p_fire);
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
			network.NeuroPopArray[pop_ind]->add_perturbation(step_perturb);
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
			for (unsigned int s = 0; s < network.ChemSynArray.size(); ++s){
				if (network.ChemSynArray[s]->get_pop_ind_pre() == pop_ind_pre &&
				network.ChemSynArray[s]->get_pop_ind_post() == pop_ind_post){ // find the right synapse object
					network.ChemSynArray[s]->add_short_term_depression(STD_on_step);
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
			for (unsigned int s = 0; s < network.ChemSynArray.size(); ++s){
				if (network.ChemSynArray[s]->get_pop_ind_pre() == pop_ind_pre &&
				network.ChemSynArray[s]->get_pop_ind_post() == pop_ind_post){ // find the right synapse object
					network.ChemSynArray[s]->add_inh_STDP(inh_STDP_on_step);
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
	
	if (spike_freq_adpt_setting.size() != 0){
		cout << "\t Spike-frequency adaptation settings...";
		for (unsigned int i = 0; i < spike_freq_adpt_setting.size(); ++i){
			int pop_ind = spike_freq_adpt_setting[i];
			network.NeuroPopArray[pop_ind]->add_spike_freq_adpt();
			cout << pop_ind+1 << "...";
		}
		cout << "done." << endl;
	}
	
	// neuron data sampling settings
	if (neuron_sample_pop_ind.size() != 0){
		cout << "\t Neuron data sampling settings...";
		for (unsigned int ind = 0; ind < neuron_sample_pop_ind.size(); ++ind){
			int pop_ind = neuron_sample_pop_ind[ind];
			network.NeuroPopArray[pop_ind]->add_sampling_real_time(neuron_sample_neurons[ind], neuron_sample_type[ind], neuron_sample_time_points[ind], out_filename);
			cout << ind+1 << "...";
		}
		cout << "done." << endl;
	}	
	
	// neuron V mean std record settings
	if (neuron_stats_setting.size() != 0){
		cout << "\t Neuron stats record settings...";
		for (unsigned int ind = 0; ind < neuron_stats_setting.size(); ++ind){
			int pop_ind = neuron_stats_setting[ind];
			network.NeuroPopArray[pop_ind]->start_stats_record();
			cout << ind+1 << "...";
		}
		cout << "done." << endl;
	}	
	
	// LFP record settings
	if (LFP_record_pop_ind.size() != 0){
		cout << "\t LFP record settings...";
		for (unsigned int ind = 0; ind < LFP_record_pop_ind.size(); ++ind){
			int pop_ind = LFP_record_pop_ind[ind];
			network.NeuroPopArray[pop_ind]->start_LFP_record(LFP_record_setting[ind]);
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
			for (unsigned int s = 0; s < network.ChemSynArray.size(); ++s){
				if (network.ChemSynArray[s]->get_pop_ind_pre() == pop_ind_pre &&
					network.ChemSynArray[s]->get_pop_ind_post() == pop_ind_post &&
				network.ChemSynArray[s]->get_syn_type() == syn_type){ // find the right synapse object
					network.ChemSynArray[s]->add_sampling(syn_sample_neurons[ind], syn_sample_time_points[ind]);
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
			for (unsigned int s = 0; s < network.ChemSynArray.size(); ++s){
				if (network.ChemSynArray[s]->get_pop_ind_pre() == pop_ind_pre &&
					network.ChemSynArray[s]->get_pop_ind_post() == pop_ind_post &&
				network.ChemSynArray[s]->get_syn_type() == syn_type){ // find the right synapse object
					network.ChemSynArray[s]->start_stats_record();
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
			network.NeuroPopArray[pop_ind]->init_runaway_killer(runaway_killer_setting[ind][1], runaway_killer_setting[ind][2], runaway_killer_setting[ind][3]);
			cout << pop_ind+1 << "...";
		}
		cout << "done." << endl;
		cout << "\t \t No women, no kids." << endl;
	}


	cout << "Importing done." << endl;
	return 1;
}

void SimuInterface::simulate(){
	clock_t begin = clock();
	  
	// simulate
	for (int step_current = 0; step_current < network.step_tot; ++step_current){
		network.update(step_current);
		/*---------------------------------------------------------------------*/
		// Countdown
		if (step_current == 0){
			cout << "Commencing countdown, engines on..." << endl << flush;	
			// if not "flush", output will be delayed in buffer
		}
		int steps_left = network.step_tot - step_current - 1;
		if ( (steps_left % (network.step_tot / 10)) == 0 ){
			clock_t end = clock();
			char str_min[80];
			double elapsed_mins = double(end - begin) / CLOCKS_PER_SEC / 60.0;
			sprintf(str_min, " (%0.1f min)...\n", elapsed_mins);
			cout <<  "\t " << steps_left / (network.step_tot / 10) << str_min << flush;
		}
		/*---------------------------------------------------------------------*/
	}
	cout << "Simulation done." << endl;

	// output results into HDF5 file
#ifdef HDF5
	if(in_filename.substr(in_filename.find_last_of(".") + 1) == "h5") 
	{
		output_results_HDF5();
	} 
	else if (in_filename.substr(in_filename.find_last_of(".") + 1) == "ygin")
	{
		output_results();
	}
#else
	// output results into text file
	output_results();
#endif
	
}

void SimuInterface::output_results(){
	// output results into text file
	ofstream output_file;
	cout << "Creating text output file...";
	output_file.open(out_filename.append(output_suffix));
	network.output_results(output_file);
	// attach input file (.ygin) to the output file for data completeness
 	ifstream in_file_attach( in_filename );
    output_file << in_file_attach.rdbuf();

	cout << "done." << endl;
	cout << "Data file name is: " << endl;
	cout << "	" << out_filename << endl;
	cout << "------------------------------------------------------------" << endl;
	// Write data file name to stdout and use "grep ygout" to extract it!
	
		
}

string SimuInterface::gen_out_filename(){
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
	convert_temp << in_filename_trim << "_" << time_stamp; // insert the textual representation of 'Number' in the characters in the stream
	return convert_temp.str(); // set 'Result' to the contents of the stream
}


template < typename Type > Type SimuInterface::read_next_entry(istringstream &line_ss){
	string entry_str;
	Type entry;
	getline(line_ss, entry_str, delim);
	if (entry_str.empty()){
		cout << "ERROR: SimuInterface::read_next_entry: empty content!" << endl;
	}
	stringstream(entry_str) >> entry;
	return entry;
}

template < typename Type, typename A > void SimuInterface::read_next_line_as_vector(vector<Type, A> &vec){
	
	string line_str, entry_str;
	Type entry;
	getline(inputfile, line_str); istringstream line_ss(line_str); // read next line
	if (line_str.empty()){
		cout << "ERROR: SimuInterface::read_next_line_as_vector: empty string!" << endl;
	}
	while (getline(line_ss, entry_str, delim)){
		stringstream(entry_str) >> entry;
		vec.push_back(entry);
	}

}




#ifdef HDF5
void SimuInterface::output_results_HDF5(){
	// output results into HDF5 file
	cout << "Outputting results into HDF5 file...";
	
	
	H5File file_HDF5;
	string file_name_HDF5 = out_filename.append("_out.h5");
	file_HDF5 = H5File( file_name_HDF5.c_str(), H5F_ACC_TRUNC );
	
	Group group_tmp = file_HDF5.createGroup(string("/config_filename"));
	write_string_HDF5(group_tmp, in_filename, string("config_filename"));
	
	network.output_results(file_HDF5);
	cout << "done." << endl;
	// Write data file name to stdout and use "grep ygout" to extract it!
	cout << "Data file name is: " << endl;
	cout << "	" << out_filename << endl;
	cout << "------------------------------------------------------------" << endl;	
}


bool SimuInterface::import_HDF5(string in_filename_input){
	in_filename = in_filename_input;
	
	const H5std_string file_name( in_filename );
	H5File file( file_name, H5F_ACC_RDONLY );
	
	out_filename = gen_out_filename();
	
	// build NeuroNet based on data imported
	if (true){
		vector<int> N_array, step_tot_v;
		read_vector_HDF5(file, string("/config/Net/INIT001/N"), N_array);
		int step_tot = read_scalar_HDF5<int>(file, string("/config/Net/INIT002/step_tot"));
		double dt = read_scalar_HDF5<double>(file, string("/config/Net/INIT002/dt"));
		network = NeuroNet( N_array, dt, step_tot, delim, indicator);
		cout << "\t Network created." << endl;
		
		
		for (unsigned int ind = 0; ind < N_array.size(); ++ind){
			network.NeuroPopArray.push_back(new NeuroPop(ind, N_array[ind], network.dt, network.step_tot, delim, indicator));
			cout << "\t Initialising neuron pop " << ind+1 << "..." << endl;

			string pop_n = "/config/pops/pop" + to_string(ind);
			
			// parameter
			if (group_exist_HDF5(in_filename, pop_n + string("/PARA001"))){
				vector<int> v_tmp;
				read_vector_HDF5(file, pop_n + string("/PARA001/para_str_ascii"), v_tmp);
				// convert ascii to string
				string para_str;
				for (unsigned i = 0; i < v_tmp.size(); i++){
					para_str += char(v_tmp[i]);
				}
				network.NeuroPopArray[ind]->set_para(para_str);
			}
			
			// external current setting
			if (group_exist_HDF5(in_filename, pop_n + string("/INIT004"))){
				cout << "\t\t External current settings...";
				vector<double> mean, std;
				read_vector_HDF5(file, pop_n + string("/INIT004/mean"), mean);
				read_vector_HDF5(file, pop_n + string("/INIT004/std"), std);
				network.NeuroPopArray[ind]->set_gaussian_I_ext(mean, std);
				cout << "done." << endl;
			}

			// external conductance setting
			if (group_exist_HDF5(in_filename, pop_n + string("/INIT012"))){
				cout << "\t\t External conductance settings...";
				vector<double> mean, std;
				read_vector_HDF5(file, pop_n + string("/INIT012/mean"), mean);
				read_vector_HDF5(file, pop_n + string("/INIT012/std"), std);
				network.NeuroPopArray[ind]->set_gaussian_g_ext(mean, std);
				cout << "done." << endl;
			}

			// random initial condition settings (double p_fire)
			if (group_exist_HDF5(in_filename, pop_n + string("/INIT003"))){
				cout << "\t\t Random initial condition settings...";
				double p_fire = read_scalar_HDF5<double>(file, pop_n + string("/INIT003/p_fire"));
				network.NeuroPopArray[ind]->random_V(p_fire);
				cout << "done." << endl;
			}
			

	
			// random initial condition settings (double r_V0, double p_fire)
			if (group_exist_HDF5(in_filename, pop_n + string("/INIT011"))){
				cout << "\t\t Random initial condition settings...";
				double r_V0 = read_scalar_HDF5<double>(file, pop_n + string("/INIT011/r_V0"));
				double p_fire = read_scalar_HDF5<double>(file, pop_n + string("/INIT011/p_fire"));
				network.NeuroPopArray[ind]->set_init_condition(r_V0, p_fire);
				cout << "done." << endl;
			}
	
			// perturbation setting
			if (group_exist_HDF5(in_filename, pop_n + string("/INIT007"))){
				cout << "\t\t Perturbation settings...";
				int step_perturb = read_scalar_HDF5<int>(file, pop_n + string("/INIT007/step_perturb"));
				network.NeuroPopArray[ind]->add_perturbation(step_perturb);
				cout << "done." << endl;
			}
			
			
			// spike-frequency adaptation setting
			if (group_exist_HDF5(in_filename, pop_n + string("/INIT010"))){
				cout << "\t\t Spike-frequency adaptation settings...";
				//int spike_freq_adpt = read_scalar_HDF5<int>(file, pop_n + string("/INIT010/spike_freq_adpt"));
				network.NeuroPopArray[ind]->add_spike_freq_adpt();
				cout << "done." << endl;
			}
			
	
			// initialise runaway-killer
			if (group_exist_HDF5(in_filename, pop_n + string("/KILL001"))){
				cout << "\t\t Runaway killer licensing for pop...";
				double min_ms = read_scalar_HDF5<double>(file, pop_n + string("/KILL001/min_ms"));
				double runaway_Hz = read_scalar_HDF5<double>(file, pop_n + string("/KILL001/runaway_Hz"));
				double Hz_ms = read_scalar_HDF5<double>(file, pop_n + string("/KILL001/Hz_ms"));
				network.NeuroPopArray[ind]->init_runaway_killer(min_ms, runaway_Hz, Hz_ms);
				cout << "done." << endl;
			}
	
			// neuron stats record settings
			if (group_exist_HDF5(in_filename, pop_n + string("/SAMP003"))){
				cout << "\t\t Neuron stats record settings...";
				network.NeuroPopArray[ind]->start_stats_record();
				cout << "done." << endl;
			}
			
			// LFP record settings
			if (group_exist_HDF5(in_filename, pop_n + string("/SAMP005"))){
				cout << "\t\t LFP record settings...";
				vector<bool> v_tmp;
				read_vector_HDF5(file, pop_n + string("/SAMP005/LFP_neurons"), v_tmp);
				// reshape
				int n_LFP = int(v_tmp.size()) / network.N_array[ind];
				vector< vector<bool> > LFP_neurons;
				LFP_neurons.resize(n_LFP);
				int ctr = 0;
				for (int n = 0; n < n_LFP; n++){
					for (int dummy = 0; dummy < network.N_array[ind]; dummy++){
						LFP_neurons[n].push_back(v_tmp[ctr]);
						ctr++;
					}
				}
				network.NeuroPopArray[ind]->start_LFP_record(LFP_neurons);
				cout << "done." << endl;
			}
			

			// neuron data sampling settings
			if (group_exist_HDF5(in_filename, pop_n + string("/SAMP001"))){
				cout << "\t\t Neuron data sampling settings...";
				vector<bool> time_points, type;
				vector<int> neurons;
				read_vector_HDF5(file, pop_n + string("/SAMP001/neurons"),neurons);
				read_vector_HDF5(file, pop_n + string("/SAMP001/time_points"),time_points);
				
				type.resize(8);
				type[0] = read_scalar_HDF5<bool>(file, pop_n + string("/SAMP001/data_type/V"));
				type[1] = read_scalar_HDF5<bool>(file, pop_n + string("/SAMP001/data_type/I_leak"));
				type[2] = read_scalar_HDF5<bool>(file, pop_n + string("/SAMP001/data_type/I_AMPA"));
				type[3] = read_scalar_HDF5<bool>(file, pop_n + string("/SAMP001/data_type/I_GABA"));
				type[4] = read_scalar_HDF5<bool>(file, pop_n + string("/SAMP001/data_type/I_NMDA"));
				type[5] = read_scalar_HDF5<bool>(file, pop_n + string("/SAMP001/data_type/I_GJ"));
				type[6] = read_scalar_HDF5<bool>(file, pop_n + string("/SAMP001/data_type/I_ext"));
				type[7] = read_scalar_HDF5<bool>(file, pop_n + string("/SAMP001/data_type/I_K"));
				
				network.NeuroPopArray[ind]->add_sampling_real_time_HDF5(neurons, type, time_points, out_filename);
				cout << "done." << endl;
			}	

			
			// cout << "\t done." << endl;
		}
		
	}

	
	if (true){
		
		int n_syns = read_scalar_HDF5<int>(file, string("/config/syns/n_syns"));
		for (int ind = 0; ind < n_syns; ++ind){
			cout << "\t Initialising chemical synapses ";
			
			string syn_n = "/config/syns/syn" + to_string(ind);
			// chemical connections
			if (group_exist_HDF5(in_filename, syn_n + string("/INIT006"))){
				
				cout << ind+1 << "..." << endl;
				
				int type = read_scalar_HDF5<int>(file, syn_n + string("/INIT006/type"));
				int i_pre = read_scalar_HDF5<int>(file, syn_n + string("/INIT006/i_pre"));
				int j_post = read_scalar_HDF5<int>(file, syn_n + string("/INIT006/j_post"));
				vector<int> I, J;
				vector<double> K, D;
				read_vector_HDF5(file, syn_n + string("/INIT006/I"), I);
				read_vector_HDF5(file, syn_n + string("/INIT006/J"), J);
				read_vector_HDF5(file, syn_n + string("/INIT006/K"), K);
				read_vector_HDF5(file, syn_n + string("/INIT006/D"), D);
				network.ChemSynArray.push_back(new ChemSyn(network.dt, network.step_tot, delim, indicator));
				network.ChemSynArray.back()->init(type, i_pre, j_post, network.N_array[i_pre], network.N_array[j_post], I, J, K, D);
				
				// parameters
				if (group_exist_HDF5(in_filename, string("/config/syns/PARA002"))){
					cout << "\t\t Synapse parameter settings...";
					vector<int> v_tmp;
					read_vector_HDF5(file, string("/config/syns/PARA002/para_str_ascii"), v_tmp);
					// convert ascii to string
					string para_str;
					for (unsigned i = 0; i < v_tmp.size(); i++){
						para_str += char(v_tmp[i]);
					}
					network.ChemSynArray.back()->set_para(para_str);
					cout << "done." << endl;
				}
					
				if (group_exist_HDF5(in_filename, string("/config/syns/INIT013"))){
					cout << "\t\t Synaptic dynamics model choice...";
					int model_choice = read_scalar_HDF5<int>(file, string("/config/syns/INIT013/model_choice"));
					network.ChemSynArray.back()->set_synapse_model(model_choice);
					cout << "done." << endl;
				}
				
				// STD
				if (group_exist_HDF5(in_filename, syn_n + string("/INIT008"))){
					cout << "\t\t Short-term depression settings...";
					int STD_on_step = read_scalar_HDF5<int>(file, syn_n + string("/INIT008/STD_on_step"));
					network.ChemSynArray.back()->add_short_term_depression(STD_on_step);
					cout << "done." << endl;
				}
				
				
				// Inhibitory STDP
				if (group_exist_HDF5(in_filename, syn_n + string("/INIT009"))){
					cout << "\t\t Inhibitory STDP settings...";
					int inh_STDP_on_step = read_scalar_HDF5<int>(file, syn_n + string("/INIT009/inh_STDP_on_step"));
					network.ChemSynArray.back()->add_inh_STDP(inh_STDP_on_step);
					cout << "done." << endl;
				}
				
				// syn I mean std record settings
				if (group_exist_HDF5(in_filename, syn_n + string("/SAMP004"))){
					cout << "\t\t Synapse stats record settings...";
					network.ChemSynArray.back()->start_stats_record();
					cout << "done." << endl;
				}	
				
				// syn data sampling settings
				if (group_exist_HDF5(in_filename, syn_n + string("/SAMP002"))){
					cout << "\t\t Synapse data sampling settings...";
					vector<int> neurons;
					vector<bool> time_points;
					read_vector_HDF5(file, syn_n + string("/SAMP002/neurons"), neurons);
					read_vector_HDF5(file, syn_n + string("/SAMP002/time_points"), time_points);
					network.ChemSynArray.back()->add_sampling(neurons, time_points);
					cout << "done." << endl;
				}	
				
				
			}
			// external spike setting
			else if (group_exist_HDF5(in_filename, syn_n + string("/INIT005"))){
				cout << ind+1 << " (ext spike)..." << endl;
				
				int j_post = read_scalar_HDF5<int>(file, syn_n + string("/INIT005/pop_ind"));
				int type_ext = read_scalar_HDF5<int>(file, syn_n + string("/INIT005/type_ext"));
				double K_ext = read_scalar_HDF5<double>(file, syn_n + string("/INIT005/K_ext"));
				int Num_ext = read_scalar_HDF5<int>(file, syn_n + string("/INIT005/Num_ext"));
				vector<bool> neurons;
				read_vector_HDF5(file, syn_n + string("/INIT005/neurons"), neurons);
				vector<double> rate_ext_t;
				read_vector_HDF5(file, syn_n + string("/INIT005/rate_ext_t"), rate_ext_t);
				network.ChemSynArray.push_back(new ChemSyn(network.dt, network.step_tot, delim, indicator));
				network.ChemSynArray.back()->init(type_ext, j_post, network.N_array[j_post], K_ext, Num_ext, rate_ext_t, neurons);
				// network.ChemSynArray.back()->set_para(syn_para);
				
			}

			// cout << "\t done." << endl;
			
		}
		
	}

	
	
	cout << "Importing done." << endl;
	return 1;
}


/*
void SimuInterface::read_matrix_HDF5(const H5File & file, const string & name, vector< vector <double> > & m_tmp){
	const H5std_string dataset_name( name );
	DataSet dataset = file.openDataSet( dataset_name );
	DataSpace dataspace = dataset.getSpace();
	
	
	hsize_t dims_out[2];
	dataspace.getSimpleExtentDims( dims_out, NULL);
	cout << dims_out[0] << "," <<  dims_out[1] << endl;

    hsize_t fdims[2];            // new data dimensions 
	fdims[0] = 1;
	fdims[1] = dims_out[1];
	
	m_tmp.resize(dims_out[0]);
	for (int i = 0; i < int(dims_out[0]); i++ ){
		hsize_t offset[2];
		offset[0] = i;
		offset[1] = 0;
		
		// selectHyperslab not working properly? 
		cout << fdims[0]<< "," <<  fdims[1] << endl;
		cout << offset[0] << "," <<  offset[1] << endl;
		
		DataSpace memspace = dataset.getSpace();
		memspace.selectHyperslab( H5S_SELECT_SET, fdims, offset );
		
		hsize_t m_out[2];
		memspace.getSimpleExtentDims(  m_out, NULL);
		cout <<  m_out[0] << "," <<   m_out[1] << endl;
		// selectHyperslab not working properly? 
		
		m_tmp[i].resize(dims_out[1]);
		cout << "here3" << endl;
		
		dataset.read(  m_tmp[i].data(), PredType::NATIVE_DOUBLE, memspace, dataspace );	
		cout << "here4" << endl;
	}
	
	cout << "here" << endl;
	for (int i = 0; i < int(m_tmp.size()); i++ ){
		for (int j = 0; j < int(m_tmp[i].size()); j++ ){
			cout << m_tmp[i][j] << ","; cout.flush();
		}
		cout << endl;
	} 
}
*/



#endif


