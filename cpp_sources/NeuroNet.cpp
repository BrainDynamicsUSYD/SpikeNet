#include "NeuroNet.h"
//#include <iostream>
//#include <random>
//#include <vector>
#include <iterator>     // std::back_inserter
#include <algorithm>    // std::for_each, copy
#include <numeric>      // std::accumulate
#include <string.h>     // for memcpy
#include <ctime>

using namespace std;

NeuroNet::NeuroNet(vector<int> N_array_input, double dt_input, int step_tot_input){

	N_array = N_array_input;
	dt = dt_input;
	step_tot = step_tot_input;
	Num_pop = N_array.size();
	runaway_killed = false;
	step_killed = -1;
	
	
};


void NeuroNet::update(int step_current){

	if (runaway_killed == false){ // if not runaway_killed
		/*------------------------------------------------------------------------------------------------------*/
	
		// Update neuron states
		// Pre-coupling update
	
		for (int pop_ind = 0; pop_ind < Num_pop; ++pop_ind){
			// Update spikes and nonref
			NeuroPopArray[pop_ind]->update_spikes(step_current); // this must always be the first operation!!!
			// Update action potential for electrical coupling
			// ElectricalSynapsesArray[i_pre][i_pre].update_action_potential(NeuroPopArray[i_pre], step_current);

			// check runaway activity
			runaway_killed |= NeuroPopArray[pop_ind]->get_runaway_killed(); // Bitwise OR assignment, a |= b meaning  a = a | b
			if (runaway_killed == true){
				step_killed = step_current;
			}
		}

		/*------------------------------------------------------------------------------------------------------*/
		// Chemical Coupling
		for (unsigned int syn_ind = 0; syn_ind < ChemSynArray.size(); ++syn_ind){
				ChemSynArray[syn_ind]->recv_pop_data(NeuroPopArray);
				// recv_data should be more optimized if using MPI!
				ChemSynArray[syn_ind]->update(step_current);
				ChemSynArray[syn_ind]->send_pop_data(NeuroPopArray);
		}
		
		// Electrical coupling
		
		/*------------------------------------------------------------------------------------------------------*/
// JH Learning Scheme
		for (int pop_ind = 0; pop_ind < Num_pop; ++pop_ind){
			NeuroPopArray[pop_ind]->reset_rhat();
		}

		for (unsigned int syn_ind = 0; syn_ind < ChemSynArray.size(); ++syn_ind){	
			ChemSynArray[syn_ind]->record_V_post_JH_Learn(NeuroPopArray);
			ChemSynArray[syn_ind]->update_post_spike_hist_JH_Learn(); 
			ChemSynArray[syn_ind]->new_post_spikes_JH_Learn(); 
			ChemSynArray[syn_ind]->old_pre_spikes();
			if(ChemSynArray[syn_ind]->get_direction()==0){	
				ChemSynArray[syn_ind]->get_rhat_spiking(); //include all synpases? CHANGE?
				ChemSynArray[syn_ind]->wchange_non_Hebbian_outgoing(NeuroPopArray);	
				ChemSynArray[syn_ind]->update_Vint_JH_Learn();//goes at end	
			}
			else if(ChemSynArray[syn_ind]->get_direction()==1){
				//TODO
			}	
		}
 // cout <<"DEBUG3\n";
		for (unsigned int syn_ind = 0; syn_ind < ChemSynArray.size(); ++syn_ind){
			if(ChemSynArray[syn_ind]->get_direction()==0){
				//Apply U rule
				ChemSynArray[syn_ind]->wchange_Hebbian_outgoing();  
			}
			ChemSynArray[syn_ind]->new_pre_spikes_JH_Learn(); // goes after old_pre_spikes which increments t_hist	
		}
		// for JH learn scheme rhat sampling (like all sampling) happens during update_V
		// so need to update rhat's before the call to update_V below
		for (int pop_ind = 0; pop_ind < Num_pop; ++pop_ind){
			NeuroPopArray[pop_ind]->get_all_rhat_JHLearn(ChemSynArray, step_current);
		}
		/*------------------------------------------------------------------------------------------------------*/
		// Update membrane potential
		// Post-coupling update
		for (int pop_ind = 0; pop_ind < Num_pop; ++pop_ind){
			// Update membrane potential
			NeuroPopArray[pop_ind]->update_V(step_current); 
		}

	} // if not runaway_killed

}

void NeuroNet::import_restart(H5File & file,string out_filename){

	Num_pop=read_scalar_HDF5<int>(file,string("/Net/Num_pop"));
	for (unsigned int ind = 0; ind < N_array.size(); ++ind){
		NeuroPopArray.push_back(new NeuroPop(ind, N_array[ind], dt, step_tot));
		cout << "\t Initialising neuron pop " << ind+1 << "..." << endl;

		NeuroPopArray[ind]->import_restart(file,ind,out_filename);
	}

	string syn_str = "/syns/";
	int n_syns=read_scalar_HDF5<int>(file,syn_str+"n_syns");
	for (int ind = 0; ind < n_syns; ++ind){

		ChemSynArray.push_back(new ChemSyn(dt, step_tot));
		// network.ChemSynArray.back()->init(type, i_pre, j_post, N_array[i_pre], N_array[j_post], I, J, K, D);

		ChemSynArray[ind]->import_restart(file,ind);
	}

}

void NeuroNet::export_restart(H5File & file, int restart_no){
	Group group_restart = file.createGroup(string("/Restart/"));
	write_scalar_HDF5(group_restart,restart_no,string("child_no_of_parent")); 
	write_scalar_HDF5(group_restart,0,string("no_children")); 

	Group group_Net = file.createGroup(string("/Net"));
	write_vector_HDF5(group_Net,N_array,string("N_array"));
	write_scalar_HDF5(group_Net,step_tot,string("step_tot")); 
	write_scalar_HDF5(group_Net,dt,string("dt")); 
	write_scalar_HDF5(group_Net,Num_pop,string("Num_pop")); 
		
	string pops_str = "/pops/";
	Group group_pops = file.createGroup(pops_str);
	for (int ind = 0; ind < Num_pop; ++ind){
		NeuroPopArray[ind]->export_restart(group_pops);
	}
	string syn_str = "/syns/";
	Group group_syns = file.createGroup(syn_str);
	int n_syns=static_cast<int>(ChemSynArray.size());
	write_scalar_HDF5(group_syns,n_syns,"n_syns");
	for (int ind = 0; ind < n_syns; ++ind){
		ChemSynArray[ind]->export_restart(group_syns,ind);
	}

}

void NeuroNet::output_results(H5File & file_HDF5){

	// KILL002 # step at which runaway activity is killed
	Group group_tmp = file_HDF5.createGroup(string("/run_away_killed"));
	vector<int> v_tmp; v_tmp.push_back(step_killed);
	write_vector_HDF5(group_tmp, v_tmp, string("step"));
		
	// dump population data
	for (int i = 0; i < Num_pop; i++){
		NeuroPopArray[i]->output_results(file_HDF5);
	}

	// dump synapse data
	for (unsigned int i = 0; i < ChemSynArray.size(); i++){
		ChemSynArray[i]->output_results(file_HDF5, i);
	}
	
}
