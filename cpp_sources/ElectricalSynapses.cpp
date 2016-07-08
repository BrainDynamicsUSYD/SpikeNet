#include "ElectricalSynapses.h"
//#include <iostream>
#include <math.h> // for round, sqrt, ceil


ElectricalSynapses::ElectricalSynapses(double dt_input) :dt(dt_input){
	// Default parameters
	spike_peak = 0, // (mV), height of spike, artificial value only to model EPSC/IPSC
	spike_duration = 0.5; // (ms), duration of spike
	boundary = 1; // 1 = "Periodic", 2 = "Free"
};

void ElectricalSynapses::update_action_potential(Neurons &pop, int step_current){
	// Model the spike by action potential (another approach is modelling current instead of potential)

	if (!C.empty() && g_GJ > 0){ // If synapses exist and strength is not zero
		// Check that spike duration is shorter than refractory period
		if (spike_steps > pop.ref_steps){
			cout << "Spike duration cannot be longer than refractory period!" << endl;
			cout << "spike_steps is: " << spike_steps << endl;
			cout << "pop.ref_steps is :" << pop.ref_steps << endl;
		}

		// The shape of action potential is assumed to be square (impulse-like) for simplicity
		int i_firing;
		int t_ring_firing = int(step_current % pop.history_steps);
		for (int ind = 0; ind < pop.num_spikes_pop[step_current]; ++ind){
			i_firing = pop.spikes_pop[t_ring_firing][ind];
			pop.V[i_firing] = spike_peak;
			// recall that only the potential of a non-refractory neuron will be updated
			// So, within refractory period, once set to peak, it remains at peak unless reset by the following code, rendering it square-shaped
		}

		// Reset the neuron to rest potential if action potential has finished
		if (step_current - spike_steps >= 0){
			int i_reset;
			int t_ring_reset = int((step_current - spike_steps) % pop.history_steps);
			for (int ind = 0; ind < pop.num_spikes_pop[step_current - spike_steps]; ++ind){
				i_reset = pop.spikes_pop[t_ring_reset][ind];
				pop.V[i_reset] = pop.V_rt;
			}
		}
	}
}

void ElectricalSynapses::update_electrical_coupling(Neurons &pre_pop, Neurons &post_pop){
	// I_GJ_post = g_GJ*(V_pre-V_post)�� summation over all the synapses

	// For heterogeneous coupling
	if (!C.empty() && g_GJ > 0){ // If synapses exist and strength is not zero
		int j_post;
		for (unsigned int i_pre = 0; i_pre < C.size(); ++i_pre){
			// Note that something funny is happening here:
			// interaction via gap junction is essentially "diffusion" of membrane potential where ion number conserves.
			// But the interaction will be no longer symmetric if happening between one refractory neuron and one non-ref neuron,
			// where the latter "feels" it but the former not,
			// As a result, ions are not conservative during the interaction!
			for (unsigned int ind = 0; ind < C[i_pre].size(); ++ind){
				j_post = C[i_pre][ind];
				post_pop.I_GJ[j_post] += g_GJ*(pre_pop.V[i_pre] - post_pop.V[j_post]);
			}
		}
	}


	// // For homogeneous coupling
	//if (!C.empty() && g_GJ > 0 ){ // If synapses exist and strength is not zero
	//	int num_synapses = C[0].size(); // recall that the connection is homogeneous in our model, i.e., same number of synapses for every neuron.
	//	for (int i_pre = 0; i_pre < C.size(); ++i_pre){
	//		// Note that something funny is happening here:
	//		// interaction via gap junction is essentially "diffusion" of membrane potential where ion number conserves.
	//		// But the interaction will be no longer symmetric if happening between one refractory neuron and one non-ref neuron,
	//		// where the latter "feels" it but the former not,
	//		// As a result, ions are not conservative during the interaction!
	//		for (int j_post = 0; j_post < num_synapses; ++j_post){
	//			post_pop.I_GJ[C[i_pre][j_post]] += pre_pop.V[i_pre];
	//		}
	//	}
	//	for (int j_post = 0; j_post < post_pop.N; ++j_post){
	//		post_pop.I_GJ[j_post] = g_GJ*(post_pop.I_GJ[j_post] - num_synapses*post_pop.V[j_post]);
	//		// Generally speaking, the more gap junctions one neuron has, the faster the diffusion is.
	//		// On the other hand, coupling by spike and leaky current "dissipate" or "creat" membrane potential.
	//	}
	//	//cout << "sample I_GJ is " << post_pop.I_GJ[60 * 50] << endl;
	//}
}

/*void ElectricalSynapses::init_homo_electical_synapses(double extent_input, double g_GJ_input, int a){
	spike_steps = int(round(spike_duration / dt));

	// read parameter
	extent = extent_input;
	g_GJ = g_GJ_input;

	// build connection seed
	int R = (int)ceil(extent); // radius
	if (boundary == 1 && R >= 0.5*a){ cout << "init_homo_electical_synapses: Coupling range too big!!" << "R is " << R << " and a is " << a << endl; } // radius upper limit check for periodic boundary
	double dist_temp;
	// search points within the circle from the circumscribing square
	vector<int> x_seed, y_seed;
	for (int i = -R; i <= R; ++i){
		for (int j = -R; j <= R; ++j){
			dist_temp = sqrt(i*i + j*j);
			if (dist_temp <= extent && dist_temp != 0){
				x_seed.push_back(i);
				y_seed.push_back(j);
			}
		}
	}

	// build connection matrix
	C.resize(a*a);
	int xs, ys, i0;
	for (int x0 = 0; x0 < a; ++x0){
		for (int y0 = 0; y0 < a; ++y0){
			i0 = x0*a + y0;// from subscript to index
			C[i0].resize(x_seed.size());
			for (unsigned int s = 0; s < x_seed.size(); ++s){
				xs = x0 + x_seed[s];
				ys = y0 + y_seed[s];
				if (boundary == 1){ // Periodic boundary condition
					if (xs < 0){ xs = xs + a; };
					if (xs >= a){ xs = xs - a; };
					if (ys < 0){ ys = ys + a; };
					if (ys >= a){ ys = ys - a; };
				}
				if (boundary == 2){
					cout << "Free boundary condition is to be defined!" << endl;
				}
				C[i0][s] = xs*a + ys; // from subscript to index
			}
		}
	}

}*/
