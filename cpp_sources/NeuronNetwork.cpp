#include "NeuronNetwork.h"
//#include <iostream>
//#include <random>
//#include <vector>
#include <iterator>     // std::back_inserter
#include <algorithm>    // std::for_each, copy
#include <numeric>      // std::accumulate
#include <string.h>     // for memcpy
#include <ctime>

using namespace std;

NeuronNetwork::NeuronNetwork(vector<int> N_array_input, double dt_input, int step_tot_input, char delim_input, char indicator_input){

	N_array = N_array_input;
	dt = dt_input;
	step_tot = step_tot_input;
	delim = delim_input;
	indicator = indicator_input;
	Num_pop = N_array.size();
	runaway_killed = false;
	step_killed = -1;
	
	
};


void NeuronNetwork::update(int step_current){

	if (runaway_killed == false){ // if not runaway_killed
		/*------------------------------------------------------------------------------------------------------*/
	
		// Update neuron states
		// Pre-coupling update
	
		for (int pop_ind = 0; pop_ind < Num_pop; ++pop_ind){
			// Update spikes and nonref
			NeuronPopArray[pop_ind]->update_spikes(step_current); // this must always be the first operation!!!
			// Update action potential for electrical coupling
			// ElectricalSynapsesArray[i_pre][i_pre].update_action_potential(NeuronPopArray[i_pre], step_current);

			// check runaway activity
			runaway_killed |= NeuronPopArray[pop_ind]->runaway_killed; // Bitwise OR assignment, a |= b meaning  a = a | b
			if (runaway_killed == true){
				step_killed = step_current;
			}
		}


	
		/*------------------------------------------------------------------------------------------------------*/
		// Chemical Coupling
		for (unsigned int syn_ind = 0; syn_ind < ChemicalSynapsesArray.size(); ++syn_ind){
				ChemicalSynapsesArray[syn_ind]->recv_pop_data(NeuronPopArray);
				// recv_data should be more optimized if using MPI!
				ChemicalSynapsesArray[syn_ind]->update(step_current);
				ChemicalSynapsesArray[syn_ind]->send_pop_data(NeuronPopArray);
		}
		
		// Electrical coupling
		

		/*------------------------------------------------------------------------------------------------------*/
		// Update membrane potential
		// Post-coupling update
		for (int pop_ind = 0; pop_ind < Num_pop; ++pop_ind){
			// Update membrane potential
			NeuronPopArray[pop_ind]->update_V(step_current); 
		}

	} // if not runaway_killed

}

void NeuronNetwork::output_results(ofstream& output_file){

	// write data
	cout << "Outputting results into file..." << endl;
	
	// KILL002 # step at which runaway activity is killed
	output_file << indicator << " KILL002" << endl;
	output_file << step_killed << delim << endl;

	// dump population data
	for (int i = 0; i < Num_pop; i++){
		NeuronPopArray[i]->output_results(output_file);
	}

	// dump synapse data
	for (unsigned int i = 0; i < ChemicalSynapsesArray.size(); i++){
		ChemicalSynapsesArray[i]->output_results(output_file);
	}
	
}



