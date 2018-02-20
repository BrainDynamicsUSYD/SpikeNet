//#include "NeuroNet.h"
#include "SimuInterface.h"
#include <chrono> // #include <boost/chrono.hpp>
#include <ctime>

using namespace std;

SimuInterface::SimuInterface() {
	// define default output data file format
	output_suffix = ".ygout"; // filename extension

}

void SimuInterface::simulate() {
	clock_t begin = clock();

	// simulate
	for (int step_current = 0; step_current < network.step_tot; ++step_current) {

		/*---------------------------------------------------------------------*/
		// Countdown
		if (step_current == 0) {
			cout << "Commencing countdown, engines on..." << endl << flush;
			// if not "flush", output will be delayed in buffer
		}
		int steps_left = network.step_tot - step_current - 1;
		if ( (steps_left % (network.step_tot / 100)) == 0 ) {
			clock_t end = clock();
			char str_min[80];
			double elapsed_mins = double(end - begin) / CLOCKS_PER_SEC / 60.0;
			sprintf(str_min, " (CPU time elasped: %0.1f min; ", elapsed_mins);
			cout <<  "\t " << steps_left / (network.step_tot / 100) << str_min << flush;
			time_t rawtime;
			time (&rawtime);
			char * timeinfo;
			timeinfo = asctime(localtime(&rawtime));
			timeinfo[strlen(timeinfo) - 1] = '\0'; // get rid of endline
			cout << "Current local time: " << timeinfo << ")..." << endl;
		}
		/*---------------------------------------------------------------------*/

		network.update(step_current);
	}
	cout << "Simulation done." << endl;

	// output results into HDF5 file
	output_results_HDF5();
	export_restart_HDF5();



}

string SimuInterface::gen_out_filename() {
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

bool SimuInterface::import_restart_HDF5(string in_filename_input) {

	out_filename = gen_out_filename();
	in_filename = in_filename_input;

	const H5std_string file_name( in_filename );
	H5File file( file_name, H5F_ACC_RDONLY );

	out_filename = gen_out_filename();
	cout << "Importing Restart Config File\n";

	vector<int> N_array;
	read_vector_HDF5(file, string("/Net/N_array"), N_array);
	int step_tot = read_scalar_HDF5<int>(file, string("/Net/step_tot"));
	double dt = read_scalar_HDF5<double>(file, string("/Net/dt"));

	network = NeuroNet( N_array, dt, step_tot);
	cout << "\t Network created." << endl;
	network.import_restart(file, out_filename);
	return 1;

}

void SimuInterface::export_restart_HDF5() {
	H5File file_HDF5;
	int restart_no = 1;
	ostringstream convert_temp;
	string restart_filename_trim;

	cout << "Creating Restart Config File\n";

	istringstream restart_filename_ss(out_filename);
	getline(restart_filename_ss, restart_filename_trim, '.');

	const H5std_string file_name( in_filename );
	H5File file( file_name, H5F_ACC_RDWR );

	Group group_restart;
	string str = "/Restart/";
	if (group_exist_HDF5(file, str)) {
		restart_no = read_scalar_HDF5<int>(file, str + string("no_children"));
		restart_no++;
		//+++TODO POSSIBEL PROBLEM HERE FOR SCRIPTS DOING LOTS OF JOBS THIS NEEDS TO BE LOCKED
		hsize_t dims[1];
		dims[0] = 1;
		DataSpace fspace(1, dims);
		DataSet dataset = file.openDataSet("/Restart/no_children");
		dataset.write( &restart_no, PredType::NATIVE_INT, fspace, fspace );
		// convert_temp<< restart_filename_trim<<"_restart_1";
	}
	else {
		group_restart = file.createGroup(string("/Restart/"));
		write_scalar_HDF5(group_restart, restart_no, string("no_children"));
	}

	if (out_filename.find("restart") != string::npos) {
		convert_temp << restart_filename_trim.substr(0, restart_filename_trim.find_last_of('_')) << "_" << restart_no;
	}
	else {
		convert_temp << restart_filename_trim.substr(0, restart_filename_trim.find_last_of('_')) << "_restart_" << restart_no;
	}
	convert_temp << ".h5";

	file_HDF5 = H5File( convert_temp.str(), H5F_ACC_TRUNC );

	Group group_SimuInterface = file_HDF5.createGroup(string("/SimuInterface/"));
	write_string_HDF5(group_SimuInterface, in_filename, string("in_filename"));
	write_string_HDF5(group_SimuInterface, out_filename, string("out_filename"));

	network.export_restart(file_HDF5, restart_no);


}

void SimuInterface::output_results_HDF5() {
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

bool SimuInterface::import_HDF5(string in_filename_input) {
	in_filename = in_filename_input;

	const H5std_string file_name( in_filename );
	H5File file( file_name, H5F_ACC_RDONLY );

	out_filename = gen_out_filename();

	// build NeuroNet based on data imported
	if (true) {
		vector<int> N_array, step_tot_v;
		read_vector_HDF5(file, string("/config/Net/INIT001/N"), N_array);
		int step_tot = read_scalar_HDF5<int>(file, string("/config/Net/INIT002/step_tot"));
		double dt = read_scalar_HDF5<double>(file, string("/config/Net/INIT002/dt"));
		network = NeuroNet( N_array, dt, step_tot);
		cout << "\t Network created." << endl;


		for (unsigned int ind = 0; ind < N_array.size(); ++ind) {
			network.NeuroPopArray.push_back(new NeuroPop(ind, N_array[ind], network.dt, network.step_tot));
			cout << "\t Initialising neuron pop " << ind + 1 << "..." << endl;

			string pop_n = "/config/pops/pop" + to_string(ind);

			// parameter
			if (group_exist_HDF5(in_filename, pop_n + string("/PARA001"))) {
				vector<int> v_tmp;
				read_vector_HDF5(file, pop_n + string("/PARA001/para_str_ascii"), v_tmp);
				// convert ascii to string
				string para_str;
				for (unsigned i = 0; i < v_tmp.size(); i++) {
					para_str += char(v_tmp[i]);
				}
				network.NeuroPopArray[ind]->set_para(para_str);
			}

			// seed
			if (group_exist_HDF5(in_filename, pop_n + string("/SEED001"))) {
				cout << "\t\t RNG seed setting...";
				int seed = read_scalar_HDF5<int>(file, pop_n + string("/SEED001/seed"));
				network.NeuroPopArray[ind]->set_seed(seed);
				cout << "done." << endl;
			}

			// external current setting
			if (group_exist_HDF5(in_filename, pop_n + string("/INIT004"))) {
				cout << "\t\t External current settings...";
				vector<double> mean, std;
				read_vector_HDF5(file, pop_n + string("/INIT004/mean"), mean);
				read_vector_HDF5(file, pop_n + string("/INIT004/std"), std);
				network.NeuroPopArray[ind]->set_gaussian_I_ext(mean, std);
				cout << "done." << endl;
			}

			// external conductance setting
			if (group_exist_HDF5(in_filename, pop_n + string("/INIT012"))) {
				cout << "\t\t External conductance settings...";
				vector<double> mean, std;
				read_vector_HDF5(file, pop_n + string("/INIT012/mean"), mean);
				read_vector_HDF5(file, pop_n + string("/INIT012/std"), std);
				network.NeuroPopArray[ind]->set_gaussian_g_ext(mean, std);
				cout << "done." << endl;
			}

			// random initial condition settings (double p_fire)
			if (group_exist_HDF5(in_filename, pop_n + string("/INIT003"))) {
				cout << "\t\t Random initial condition settings...";
				double p_fire = read_scalar_HDF5<double>(file, pop_n + string("/INIT003/p_fire"));
				network.NeuroPopArray[ind]->random_V(p_fire);
				cout << "done." << endl;
			}



			// random initial condition settings (double r_V0, double p_fire)
			if (group_exist_HDF5(in_filename, pop_n + string("/INIT011"))) {
				cout << "\t\t Random initial condition settings...";
				double r_V0 = read_scalar_HDF5<double>(file, pop_n + string("/INIT011/r_V0"));
				double p_fire = read_scalar_HDF5<double>(file, pop_n + string("/INIT011/p_fire"));
				network.NeuroPopArray[ind]->set_init_condition(r_V0, p_fire);
				cout << "done." << endl;
			}
			// external initial membrane potential setting, whose prior is higher than random settings
			if (group_exist_HDF5(in_filename, pop_n + string("/SETINITV"))) {				
				cout << "\t\t External initial V settings...";
				vector<double> external_init_V;
				read_vector_HDF5(file, pop_n + string("/SETINITV/external_init_V"), external_init_V);
				network.NeuroPopArray[ind]->set_init_V_external(external_init_V);
				cout << "done." << endl;
			}



			// perturbation setting
			if (group_exist_HDF5(in_filename, pop_n + string("/INIT007"))) {
				cout << "\t\t Perturbation settings...";
				int step_perturb = read_scalar_HDF5<int>(file, pop_n + string("/INIT007/step_perturb"));
				network.NeuroPopArray[ind]->add_perturbation(step_perturb);
				cout << "done." << endl;
			}


			// spike-frequency adaptation setting
			if (group_exist_HDF5(in_filename, pop_n + string("/INIT010"))) {
				cout << "\t\t Spike-frequency adaptation settings...";
				//int spike_freq_adpt = read_scalar_HDF5<int>(file, pop_n + string("/INIT010/spike_freq_adpt"));
				network.NeuroPopArray[ind]->add_spike_freq_adpt();
				if (dataset_exist_HDF5(file, pop_n + string("/INIT010/dg_K"))) {
					cout << "reading dg_K...";
					double dg_K = read_scalar_HDF5<double>(file, pop_n + string("/INIT010/dg_K"));
					network.NeuroPopArray[ind]->set_spike_freq_adpt_para(dg_K); 
				}
				// read in heterogenous spike_freq_adpt para
				if (dataset_exist_HDF5(file, pop_n + string("/INIT010/dg_K_heter"))) {					
					cout << "reading dg_K_heter...";
					vector<double> dg_K_heter;
					int heter_SFA_start_step = read_scalar_HDF5<int>(file, pop_n + string("/INIT010/start_step"));
					int heter_SFA_end_step = read_scalar_HDF5<int>(file, pop_n + string("/INIT010/end_step"));
					read_vector_HDF5(file, pop_n + string("/INIT010/dg_K_heter"), dg_K_heter);
					network.NeuroPopArray[ind]->set_spike_freq_adpt_para_heter(dg_K_heter,heter_SFA_start_step,heter_SFA_end_step); 
				}
				cout << "done." << endl;
			}


			// initialise runaway-killer
			if (group_exist_HDF5(in_filename, pop_n + string("/KILL001"))) {
				cout << "\t\t Runaway killer licensing for pop...";
				double min_ms = read_scalar_HDF5<double>(file, pop_n + string("/KILL001/min_ms"));
				double runaway_Hz = read_scalar_HDF5<double>(file, pop_n + string("/KILL001/runaway_Hz"));
				double Hz_ms = read_scalar_HDF5<double>(file, pop_n + string("/KILL001/Hz_ms"));
				network.NeuroPopArray[ind]->init_runaway_killer(min_ms, runaway_Hz, Hz_ms);
				cout << "done." << endl;
			}

			// neuron stats record settings
			if (group_exist_HDF5(in_filename, pop_n + string("/SAMP003"))) {
				cout << "\t\t Neuron stats record settings...";
				network.NeuroPopArray[ind]->start_stats_record();
				cout << "done." << endl;
			}

			// neuron cov record settings
			if (group_exist_HDF5(in_filename, pop_n + string("/SAMP103"))) {
				cout << "\t\t Neuron cov record settings...";
				int time_start = read_scalar_HDF5<int>(file, pop_n + string("/SAMP103/time_start"));
				int time_end = read_scalar_HDF5<int>(file, pop_n + string("/SAMP103/time_end"));
				network.NeuroPopArray[ind]->start_cov_record(time_start, time_end);
				cout << "done." << endl;
			}


			// LFP record settings
			if (group_exist_HDF5(in_filename, pop_n + string("/SAMP005"))) {
				cout << "\t\t LFP record settings...";
				vector< vector<double> > LFP_neurons;
				read_matrix_HDF5(file, pop_n + string("/SAMP005/LFP_neurons"), LFP_neurons);
				network.NeuroPopArray[ind]->start_LFP_record(LFP_neurons);
				cout << "done." << endl;
			}


			// neuron data sampling settings
			if (group_exist_HDF5(in_filename, pop_n + string("/SAMP001"))) {
				cout << "\t\t Neuron data sampling settings...";
				vector<bool> time_points, type;
				vector<int> neurons;
				read_vector_HDF5(file, pop_n + string("/SAMP001/neurons"), neurons);
				read_vector_HDF5(file, pop_n + string("/SAMP001/time_points"), time_points);

				type.resize(9);
				type[0] = read_scalar_HDF5<bool>(file, pop_n + string("/SAMP001/data_type/V"));
				type[1] = read_scalar_HDF5<bool>(file, pop_n + string("/SAMP001/data_type/I_leak"));
				type[2] = read_scalar_HDF5<bool>(file, pop_n + string("/SAMP001/data_type/I_AMPA"));
				type[3] = read_scalar_HDF5<bool>(file, pop_n + string("/SAMP001/data_type/I_GABA"));
				type[4] = read_scalar_HDF5<bool>(file, pop_n + string("/SAMP001/data_type/I_NMDA"));
				type[5] = read_scalar_HDF5<bool>(file, pop_n + string("/SAMP001/data_type/I_GJ"));
				type[6] = read_scalar_HDF5<bool>(file, pop_n + string("/SAMP001/data_type/I_ext"));
				type[7] = read_scalar_HDF5<bool>(file, pop_n + string("/SAMP001/data_type/I_K"));
				type[8] = 0;
				if (dataset_exist_HDF5(file, pop_n + string("/SAMP001/data_type/rhat"))) {
					type[8] = read_scalar_HDF5<bool>(file, pop_n + string("/SAMP001/data_type/rhat"));
				}
				network.NeuroPopArray[ind]->add_sampling_real_time_HDF5(neurons, type, time_points, out_filename);
				cout << "done." << endl;
			}

			// choose neuron model
			if (dataset_exist_HDF5(file, pop_n + string("/neuron_model"))) {
				int n_mod;
				n_mod = read_scalar_HDF5<int>(file, pop_n + string("/neuron_model"));
				network.NeuroPopArray[ind]->set_neuron_model(n_mod);
			}
			if (group_exist_HDF5(in_filename, pop_n + string("/ELIF"))) {
				cout << "\t\t Exponential LIF settings...";
				double V_T, delT;
				V_T = read_scalar_HDF5<double>(file, pop_n + string("/ELIF/ELIF_VT"));
				delT = read_scalar_HDF5<double>(file, pop_n + string("/ELIF/ELIF_delT"));
				network.NeuroPopArray[ind]->set_ELIF_Params(delT, V_T);

				cout << "done." << endl;
			}

			// Poisson spiking population
			if (group_exist_HDF5(in_filename, pop_n + string("/poisson_pop"))) {
				cout << "\t\t Poisson population input settings...";
				double rate;
				rate = read_scalar_HDF5<double>(file,  pop_n + string("/poisson_pop/rate"));
				network.NeuroPopArray[ind]->init_poisson_pop(rate);
				cout << "done." << endl;
			}

			//ext_spike_input from HDF5 file
			if (group_exist_HDF5(in_filename, pop_n + string("/file_spike_input"))) {
				cout << "\t\t File spike input settings...";
				string fname;
				read_string_HDF5(file,  pop_n + string("/file_spike_input/fname"), fname);
				network.NeuroPopArray[ind]->load_file_spike_input(fname);
				cout << "done." << endl;
			}

			//ext_current_input from HDF5 file
			if (group_exist_HDF5(in_filename, pop_n + string("/file_current_input"))) {
				cout << "\t\t File current input settings...";
				string fname;
				read_string_HDF5(file,  pop_n + string("/file_current_input/fname"), fname);
				network.NeuroPopArray[ind]->load_file_current_input(fname);
				cout << "done." << endl;
			}

			// cout << "\t done." << endl;
		}

	}


	if (group_exist_HDF5(in_filename, string("/config/syns"))) {
		int n_syns = read_scalar_HDF5<int>(file, string("/config/syns/n_syns"));
		for (int ind = 0; ind < n_syns; ++ind) {
			cout << "\t Initialising chemical synapses ";

			string syn_n = "/config/syns/syn" + to_string(ind);
			// chemical connections
			if (group_exist_HDF5(in_filename, syn_n + string("/INIT006"))) {

				cout << ind + 1 << "..." << endl;

				int type = read_scalar_HDF5<int>(file, syn_n + string("/INIT006/type"));
				int i_pre = read_scalar_HDF5<int>(file, syn_n + string("/INIT006/i_pre"));
				int j_post = read_scalar_HDF5<int>(file, syn_n + string("/INIT006/j_post"));
				vector<int> I, J;
				vector<double> K, D;
				read_vector_HDF5(file, syn_n + string("/INIT006/I"), I);
				read_vector_HDF5(file, syn_n + string("/INIT006/J"), J);
				read_vector_HDF5(file, syn_n + string("/INIT006/K"), K);
				read_vector_HDF5(file, syn_n + string("/INIT006/D"), D);
				network.ChemSynArray.push_back(new ChemSyn(network.dt, network.step_tot));
				network.ChemSynArray.back()->init(type, i_pre, j_post, network.N_array[i_pre], network.N_array[j_post], I, J, K, D);

				// parameters
				if (group_exist_HDF5(in_filename, string("/config/syns/PARA002"))) {
					cout << "\t\t Synapse parameter settings...";
					vector<int> v_tmp;
					read_vector_HDF5(file, string("/config/syns/PARA002/para_str_ascii"), v_tmp);
					// convert ascii to string
					string para_str;
					for (unsigned i = 0; i < v_tmp.size(); i++) {
						para_str += char(v_tmp[i]);
					}
					network.ChemSynArray.back()->set_para(para_str);
					cout << "done." << endl;
				}

				if (group_exist_HDF5(in_filename, string("/config/syns/INIT013"))) {
					cout << "\t\t Synaptic dynamics model choice...";
					int model_choice = read_scalar_HDF5<int>(file, string("/config/syns/INIT013/model_choice"));
					network.ChemSynArray.back()->set_synapse_model(model_choice);
					cout << "done." << endl;
				}

				// seed
				if (group_exist_HDF5(in_filename, syn_n + string("/SEED002"))) {
					cout << "\t\t RNG seed setting...";
					int seed = read_scalar_HDF5<int>(file, syn_n + string("/SEED002/seed"));
					network.ChemSynArray.back()->set_seed(seed);
					cout << "done." << endl;
				}

				// STD
				if (group_exist_HDF5(in_filename, syn_n + string("/INIT008"))) {
					cout << "\t\t Short-term depression settings...";
					int STD_on_step = read_scalar_HDF5<int>(file, syn_n + string("/INIT008/STD_on_step"));
					network.ChemSynArray.back()->add_short_term_depression(STD_on_step);
					cout << "done." << endl;
				}


				// Inhibitory STDP
				if (group_exist_HDF5(in_filename, syn_n + string("/INIT009"))) {
					cout << "\t\t Inhibitory STDP settings...";
					int inh_STDP_on_step = read_scalar_HDF5<int>(file, syn_n + string("/INIT009/inh_STDP_on_step"));
					network.ChemSynArray.back()->add_inh_STDP(inh_STDP_on_step);
					cout << "done." << endl;
				}

				// JH Learning
				if (group_exist_HDF5(in_filename, syn_n + string("/INIT016"))) {
					cout << "\t\t JH Learning settings...";
					int infsteps = read_scalar_HDF5<int>(file, syn_n + string("/INIT016/infsteps"));
					double infscale = read_scalar_HDF5<double>(file, syn_n + string("/INIT016/infscale"));
					double learnrate = read_scalar_HDF5<double>(file, syn_n + string("/INIT016/learnrate"));
					double learnrateall = read_scalar_HDF5<double>(file, syn_n + string("/INIT016/learnrateall"));
					double noise = read_scalar_HDF5<double>(file, syn_n + string("/INIT016/noise"));
					int tau = read_scalar_HDF5<int>(file, syn_n + string("/INIT016/tau"));
					int ntype_pre = read_scalar_HDF5<int>(file, syn_n + string("/INIT016/ntype_pre"));
					int ntype_post = read_scalar_HDF5<int>(file, syn_n + string("/INIT016/ntype_post"));
					vector<double> dropout;
					read_vector_HDF5(file, syn_n + string("/INIT016/dropout"),dropout);
					network.NeuroPopArray[i_pre]->add_JH_Learn(dropout[0]);
					network.NeuroPopArray[j_post]->add_JH_Learn(dropout[1]);
					int direction = read_scalar_HDF5<int>(file, syn_n + string("/INIT016/direction"));
					network.ChemSynArray.back()->add_JH_Learning(network.NeuroPopArray,infsteps,infscale,learnrate,learnrateall,tau,noise,ntype_pre,ntype_post,direction);
					cout << "done." << endl;
				}


				// syn stat record settings
				if (group_exist_HDF5(in_filename, syn_n + string("/SAMP004"))) {
					cout << "\t\t Synapse stats record settings...";
					network.ChemSynArray.back()->start_stats_record();
					cout << "done." << endl;
				}

				// syn cov record settings
				if (group_exist_HDF5(in_filename, syn_n + string("/SAMP104"))) {
					cout << "\t\t Synapse stats record settings...";
					int time_start = read_scalar_HDF5<int>(file, syn_n + string("/SAMP104/time_start"));
					int time_end = read_scalar_HDF5<int>(file, syn_n + string("/SAMP104/time_end"));
					network.ChemSynArray.back()->start_cov_record(time_start, time_end);
					cout << "done." << endl;
				}


				// syn data sampling settings
				if (group_exist_HDF5(in_filename, syn_n + string("/SAMP002"))) {
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
			else if (group_exist_HDF5(in_filename, syn_n + string("/INIT005"))) {
				cout << ind + 1 << " (ext spike)..." << endl;

				int j_post = read_scalar_HDF5<int>(file, syn_n + string("/INIT005/pop_ind"));
				int type_ext = read_scalar_HDF5<int>(file, syn_n + string("/INIT005/type_ext"));
				double K_ext = read_scalar_HDF5<double>(file, syn_n + string("/INIT005/K_ext"));
				int Num_ext = read_scalar_HDF5<int>(file, syn_n + string("/INIT005/Num_ext"));
				vector<bool> neurons;
				read_vector_HDF5(file, syn_n + string("/INIT005/neurons"), neurons);
				vector<double> rate_ext_t;
				read_vector_HDF5(file, syn_n + string("/INIT005/rate_ext_t"), rate_ext_t);
				network.ChemSynArray.push_back(new ChemSyn(network.dt, network.step_tot));
				network.ChemSynArray.back()->init(type_ext, j_post, network.N_array[j_post], K_ext, Num_ext, rate_ext_t, neurons);
				// network.ChemSynArray.back()->set_para(syn_para);

			}
			
			// external spike setting (with different time invariant rates for differnt neurons) 
			else if (group_exist_HDF5(in_filename, syn_n + string("/INIT017"))) {
				cout << ind + 1 << " (ext spike t inv)..." << endl;

				int j_post = read_scalar_HDF5<int>(file, syn_n + string("/INIT017/pop_ind"));
				int type_ext = read_scalar_HDF5<int>(file, syn_n + string("/INIT017/type_ext"));
				double K_ext = read_scalar_HDF5<double>(file, syn_n + string("/INIT017/K_ext"));
				int Num_ext = read_scalar_HDF5<int>(file, syn_n + string("/INIT017/Num_ext"));
				cout << "done 1" << flush;
				vector<double> rate_ext_neuron;
				read_vector_HDF5(file, syn_n + string("/INIT017/rate_ext_neuron"), rate_ext_neuron);
				vector<bool> rate_ext_on;
				read_vector_HDF5(file, syn_n + string("/INIT017/rate_ext_on"), rate_ext_on);
				network.ChemSynArray.push_back(new ChemSyn(network.dt, network.step_tot));
				network.ChemSynArray.back()->init(type_ext, j_post, network.N_array[j_post], K_ext, Num_ext, rate_ext_on, rate_ext_neuron);
				// network.ChemSynArray.back()->set_para(syn_para);

			}
			

			// cout << "\t done." << endl;

		}

	}



	cout << "Importing done." << endl;
	return 1;
}





// string SimuInterface::gen_restart_filename(int child_no){
// 	// creat output file name using some pre-defined format
// 	ostringstream convert_temp;   // stream used for the conversion
// 	// Make use of the input file path and name
// 	string in_filename_trim;
// 	istringstream in_filename_ss(in_filename);
// 	getline(in_filename_ss, in_filename_trim, '.');

// 	// Combine all the parts of the output file path and name
// 	if(filename.find("restart")=std::string::npos){
// 		convert_temp << in_filename_trim << "_restart_0"; // insert the textual representation of 'Number' in the characters in the stream
// 	}
// 	else
// 	{
// 		in_filename_trim.pop_back();
// 		convert_temp << in_filename_trim << "_"<<child_no;
// 	}
// 	return convert_temp.str(); // set 'Result' to the contents of the stream
// }



