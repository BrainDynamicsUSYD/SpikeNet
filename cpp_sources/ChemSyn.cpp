/// \file
//#include <iostream>
#include <cmath> // always cmath instead of math.h ?
#include <algorithm>
#include <sstream>
//#include <fstream>

#include "ChemSyn.h"

void ChemSyn::add_JH_Learning(vector<NeuroPop*> &NeuronPopArray,int isteps, double iscaleE, double iscaleI,double lrate_E, double lrateall_E,double lrate_I, double lrateall_I,int intau, double innoise,int type_pre,int type_post){
	
	jh_learn_syn.on=true; //indicates this learning is to be used
	jh_learn_syn.inf_steps=isteps;
	jh_learn_syn.post_V_hist.resize(isteps+1, vector<double> (N_post,-70));  //+++TODO UPDATE init value to something else???
	jh_learn_syn.post_R_hist.resize(isteps+1, vector<int> (N_post,0)); 
	jh_learn_syn.t_ind=isteps+1; // will start by incrementing to 0
	jh_learn_syn.post_hist_len=isteps +2; //+1 for ring buffer empty space 
	jh_learn_syn.post_t_hist.resize(isteps+2,vector<int>(N_post,0)); 
	jh_learn_syn.ind_post_new.resize(N_post,isteps+1);
	jh_learn_syn.ind_post_old.resize(N_post,0);
	jh_learn_syn.pre_hist_len=isteps+2; //+1 for storing 0 to -inf_Steps, +1 for ring buffer empty space
	jh_learn_syn.pre_t_hist.resize(isteps+2,vector<int>(N_pre,0)); 
	jh_learn_syn.ind_pre_new.resize(N_pre,isteps+1);
	jh_learn_syn.ind_pre_old.resize(N_pre,0);		
	jh_learn_syn.Vint.resize(N_post,0);
	jh_learn_syn.Vint_ctr.resize(N_post,intau);
	jh_learn_syn.Q_pre.resize(N_pre,0);
	jh_learn_syn.inf_scaleI=iscaleI;
	jh_learn_syn.inf_scaleE=iscaleE;
	jh_learn_syn.learn_rate_I=lrate_I;
	jh_learn_syn.learn_rate_all_I=lrateall_I;
	jh_learn_syn.learn_rate_E=lrate_E;
	jh_learn_syn.learn_rate_all_E=lrateall_E;
	jh_learn_syn.tau=intau;
	jh_learn_syn.C=NeuronPopArray[pop_ind_post]->Cm;
	jh_learn_syn.noise=innoise;
	jh_learn_syn.ntype_pre=type_pre;
	jh_learn_syn.ntype_post=type_post;
	//+++TODO REUSE j_2_i FROM Inh_STDP
	jh_learn_syn.j_2_i.resize(N_post);
	jh_learn_syn.j_2_syn_ind.resize(N_post);
	int j_post;
	for (int i_pre = 0; i_pre < N_pre; ++i_pre){ 
		for (unsigned int syn_ind = 0; syn_ind < C[i_pre].size(); ++syn_ind){
			j_post = C[i_pre][syn_ind];
			jh_learn_syn.j_2_i[j_post].push_back( i_pre );
			jh_learn_syn.j_2_syn_ind[j_post].push_back( int(syn_ind) );
		}
	}
}


void ChemSyn::record_V_post_JH_Learn(vector<NeuroPop*> &NeuronPopArray){
	if(jh_learn_syn.on){
		jh_learn_syn.t_ind=(jh_learn_syn.t_ind+1)%jh_learn_syn.post_V_hist.size();
		jh_learn_syn.post_V_hist[jh_learn_syn.t_ind]=V_post;
		jh_learn_syn.post_R_hist[jh_learn_syn.t_ind]=NeuronPopArray[pop_ind_post]->ref_step_left;
	}
}

void ChemSyn::update_post_spike_hist_JH_Learn(){
	if(jh_learn_syn.on){
		int j_post,i;
		for(j_post=0;j_post<N_post;j_post++){
			for(i=jh_learn_syn.ind_post_old[j_post];i!=(jh_learn_syn.ind_post_new[j_post]+1)%jh_learn_syn.post_hist_len;i=(i+1)%jh_learn_syn.post_hist_len){
				jh_learn_syn.post_t_hist[i][j_post]++;
				if(jh_learn_syn.post_t_hist[i][j_post]>=jh_learn_syn.inf_steps){
					//if spike is too old increment index past it
					jh_learn_syn.ind_post_old[j_post]=(jh_learn_syn.ind_post_old[j_post]+1)%jh_learn_syn.post_hist_len;
				}
			}
		}
	}
}


void ChemSyn::update_Vint_JH_Learn(){
	if(jh_learn_syn.on){
		int j_post;
		for(j_post=0;j_post<N_post;j_post++){
			jh_learn_syn.Vint[j_post]*=exp(-1/jh_learn_syn.tau);
			if(jh_learn_syn.post_R_hist[jh_learn_syn.t_ind][j_post]==0){

				jh_learn_syn.Vint_ctr[j_post]++;
				if(syn_type==0){
					jh_learn_syn.Vint[j_post]-=(V_post[j_post]-V_ex);
				}
				else{
					jh_learn_syn.Vint[j_post]+=(V_post[j_post]-V_in);
				}
			}
		}
	}
}


void ChemSyn::new_post_spikes_JH_Learn(){
	if(jh_learn_syn.on){
		int  j_post, i_pre;
		unsigned int syn_ind,ind, i;
		//update new post spikes
		for(i=0; i<spikes_post.size(); i++){
			j_post=spikes_post[i];
			jh_learn_syn.ind_post_new[j_post]=(jh_learn_syn.ind_post_new[j_post]+1)%jh_learn_syn.post_hist_len;
			jh_learn_syn.post_t_hist[jh_learn_syn.ind_post_new[j_post]][j_post]=0;
			
			//update weights with Vint
			for( ind=0; ind< jh_learn_syn.j_2_i[j_post].size(); ind++){
				i_pre=jh_learn_syn.j_2_i[j_post][ind];
				syn_ind=jh_learn_syn.j_2_syn_ind[j_post][ind];
				if(syn_type==0){
					K[i_pre][syn_ind]-=jh_learn_syn.learn_rate_E*(1.0/jh_learn_syn.Vint_ctr[j_post])*(1.0/jh_learn_syn.C)*jh_learn_syn.Vint[j_post];
				}
				else{
					K[i_pre][syn_ind]-=jh_learn_syn.learn_rate_I*(1.0/jh_learn_syn.Vint_ctr[j_post])*(1.0/jh_learn_syn.C)*jh_learn_syn.Vint[j_post];
				}

				if(K[i_pre][syn_ind]<0){
					K[i_pre][syn_ind]=0;
				}
			}
			jh_learn_syn.Vint[j_post]=0; //reset Vint
			jh_learn_syn.Vint_ctr[j_post]=0; //reset Vint
		}	
		if(((jh_learn_syn.learn_rate_all_E>0)&&(syn_type==0))||((jh_learn_syn.learn_rate_all_I>0)&&(syn_type==1))){
			//also modify weights on nonspiking neurons
			for(j_post=0; j_post<N_post; j_post++){
				jh_learn_syn.ind_post_new[j_post]=(jh_learn_syn.ind_post_new[j_post]+1)%jh_learn_syn.post_hist_len;
				jh_learn_syn.post_t_hist[jh_learn_syn.ind_post_new[j_post]][j_post]=0;
				
				//update weights with Vint
				for( ind=0; ind< jh_learn_syn.j_2_i[j_post].size(); ind++){
					i_pre=jh_learn_syn.j_2_i[j_post][ind];
					syn_ind=jh_learn_syn.j_2_syn_ind[j_post][ind];
					
					if(syn_type==0){
						K[i_pre][syn_ind]-=jh_learn_syn.learn_rate_all_E;
					}
					else{
						K[i_pre][syn_ind]-=jh_learn_syn.learn_rate_all_I;
					}
					if(K[i_pre][syn_ind]<0){
						K[i_pre][syn_ind]=0;
					}
				}
		
			}
		}
	}
}


void ChemSyn::new_pre_spikes_JH_Learn(){
	if(jh_learn_syn.on){
		int  i_pre;
		unsigned int i;
		//update new post spikes
		for(i=0; i<spikes_pre.size(); i++){
			i_pre=spikes_pre[i];
			jh_learn_syn.ind_pre_new[i_pre]=(jh_learn_syn.ind_pre_new[i_pre]+1)%jh_learn_syn.pre_hist_len;
			jh_learn_syn.pre_t_hist[jh_learn_syn.ind_pre_new[i_pre]][i_pre]=0;
		}
	}
}

void ChemSyn::old_pre_spikes_Q_JH_Learn(vector<NeuroPop*> &NeuronPopArray){
	if(jh_learn_syn.on){
		int i_pre,j_post,i,inf_t_ind,age;
		unsigned int syn_ind;
		double expdecay,Q;
		vector <double> post_h;
		jh_learn_syn.old_pre.clear();
		for(i_pre=0;i_pre<N_pre;i_pre++){
			for(i=jh_learn_syn.ind_pre_old[i_pre]; i!=(jh_learn_syn.ind_pre_new[i_pre]+1)%jh_learn_syn.pre_hist_len;i=(i+1)%jh_learn_syn.pre_hist_len){
				//increment time of history pre spikes
				jh_learn_syn.pre_t_hist[i][i_pre]++;
			}

			if(	(((jh_learn_syn.ind_pre_new[i_pre]+1)%jh_learn_syn.pre_hist_len!=jh_learn_syn.ind_pre_old[i_pre]) &&
				(jh_learn_syn.pre_t_hist[jh_learn_syn.ind_pre_old[i_pre]][i_pre]>=jh_learn_syn.inf_steps))	){
				// IF random noise spike generated for learning OR
				// the oldest pre spike is too old 
				//THEN  update weights 
				jh_learn_syn.old_pre.push_back(i_pre);
				// Calculate Q value
				Q=0;
				post_h.resize(C[i_pre].size(),0);
				for (syn_ind = 0; syn_ind < C[i_pre].size(); ++syn_ind){
					j_post = C[i_pre][syn_ind];
					if(jh_learn_syn.ind_post_old[j_post]!=(jh_learn_syn.ind_post_new[j_post]+1)%jh_learn_syn.post_hist_len){
						//there is a j_post spike in history
						age=jh_learn_syn.post_t_hist[jh_learn_syn.ind_post_old[j_post]][j_post];
						inf_t_ind=(jh_learn_syn.t_ind+jh_learn_syn.post_V_hist.size()-jh_learn_syn.inf_steps)%(jh_learn_syn.post_V_hist.size()); 
						if(jh_learn_syn.post_R_hist[inf_t_ind][j_post]==0){
							expdecay=exp(-(jh_learn_syn.inf_steps-1-age)/jh_learn_syn.tau);

							if(syn_type==0){
								post_h[syn_ind]=-(jh_learn_syn.post_V_hist[inf_t_ind][j_post]-V_ex)*expdecay;
							}
							else
							{
								post_h[syn_ind]=-(jh_learn_syn.post_V_hist[inf_t_ind][j_post]-V_in)*expdecay;
							}
							Q+=post_h[syn_ind]*K[i_pre][syn_ind];
						}
					}
				}
				if(jh_learn_syn.ntype_post==0){
					NeuronPopArray[pop_ind_pre]->jh_learn_pop.QE[i_pre]+=(1/jh_learn_syn.C)*Q;
				}
				else{
					NeuronPopArray[pop_ind_pre]->jh_learn_pop.QI[i_pre]+=(1/jh_learn_syn.C)*Q;
				}
			}

			if(((jh_learn_syn.ind_pre_new[i_pre]+1)%jh_learn_syn.pre_hist_len!=jh_learn_syn.ind_pre_old[i_pre]) &&
			(jh_learn_syn.pre_t_hist[jh_learn_syn.ind_pre_old[i_pre]][i_pre]>=jh_learn_syn.inf_steps))
			{
				// there is a spike in history AND the oldest pre spike is too old and needs to be removed and weights updated
				// dont do this if it was a noise spike -- so not in the IF statement above
				jh_learn_syn.ind_pre_old[i_pre]=(jh_learn_syn.ind_pre_old[i_pre]+1)%jh_learn_syn.pre_hist_len;
			}
		}
	}
}


void ChemSyn::old_pre_spikes_K_JH_Learn(vector<NeuroPop*> &NeuronPopArray){
	if(jh_learn_syn.on){
		int i_pre,j_post,inf_t_ind,age;
		unsigned int syn_ind,i;
		double expdecay,Q;
		vector <double> post_h;
		
		for(i=0;i<jh_learn_syn.old_pre.size();i++){
			i_pre=jh_learn_syn.old_pre[i];	
			if((syn_type==0)&(jh_learn_syn.ntype_post==0)){
				Q=jh_learn_syn.inf_scaleE*NeuronPopArray[pop_ind_pre]->jh_learn_pop.QE[i_pre];
			}
			else if ((syn_type==0)&(jh_learn_syn.ntype_post==1)){
				Q=jh_learn_syn.inf_scaleE*NeuronPopArray[pop_ind_pre]->jh_learn_pop.QI[i_pre];
			}
			else if ((syn_type==1)&(jh_learn_syn.ntype_post==0)){
				Q=-jh_learn_syn.inf_scaleI*NeuronPopArray[pop_ind_pre]->jh_learn_pop.QE[i_pre];
			}
			else if ((syn_type==1)&(jh_learn_syn.ntype_post==1)){
				Q=-jh_learn_syn.inf_scaleI*NeuronPopArray[pop_ind_pre]->jh_learn_pop.QI[i_pre];
			}	

			// cout<<"\r"<<"syn type: "<<syn_type<<" pop pre:"<<pop_ind_pre<<" pop post:"<<pop_ind_post<<" Q: "<<Q;

			// If Q is less than noise then make it noise level
			if(Q<jh_learn_syn.noise){
				Q=jh_learn_syn.noise;
			}

			post_h.resize(C[i_pre].size(),0);
			for (syn_ind = 0; syn_ind < C[i_pre].size(); ++syn_ind){
				j_post = C[i_pre][syn_ind];
				if(jh_learn_syn.ind_post_old[j_post]!=(jh_learn_syn.ind_post_new[j_post]+1)%jh_learn_syn.post_hist_len){
					//there is a j_post spike in history
					age=jh_learn_syn.post_t_hist[jh_learn_syn.ind_post_old[j_post]][j_post];
					inf_t_ind=(jh_learn_syn.t_ind+jh_learn_syn.post_V_hist.size()-jh_learn_syn.inf_steps)%(jh_learn_syn.post_V_hist.size()); 
					if(jh_learn_syn.post_R_hist[inf_t_ind][j_post]==0){

						expdecay=exp(-(jh_learn_syn.inf_steps-1-age)/jh_learn_syn.tau); 
						if(syn_type==0){
							post_h[syn_ind]=-(jh_learn_syn.post_V_hist[inf_t_ind][j_post]-V_ex)*expdecay;
							K[i_pre][syn_ind]+=jh_learn_syn.learn_rate_E*(1.0/jh_learn_syn.C)*post_h[syn_ind]/Q;
						}
						else
						{
							post_h[syn_ind]=(jh_learn_syn.post_V_hist[inf_t_ind][j_post]-V_in)*expdecay;
							K[i_pre][syn_ind]+=jh_learn_syn.learn_rate_I*(1.0/jh_learn_syn.C)*post_h[syn_ind]/Q;	
						}
						if(K[i_pre][syn_ind]<0){
							K[i_pre][syn_ind]=0;
						}	
					}
				}
				if(((jh_learn_syn.learn_rate_all_E>0)&&(syn_type==0))||((jh_learn_syn.learn_rate_all_I>0)&&(syn_type==1))){
					//update weights to non spiking neurons as well
					if(syn_type==0){
						K[i_pre][syn_ind]+=jh_learn_syn.learn_rate_all_E/Q;
					}
					else{
						K[i_pre][syn_ind]+=jh_learn_syn.learn_rate_all_I/Q;
					}
					if(K[i_pre][syn_ind]<0){
						K[i_pre][syn_ind]=0;
					}				
				}	
			}
		}
	}
}

ChemSyn::ChemSyn(const double dt_input, const int step_tot_input, const char delim_input, const char indicator_input){
	
	dt = dt_input;
	step_tot = step_tot_input;
	delim = delim_input;
	indicator = indicator_input;
	
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

	// default settting
	stats.record = false;
	STD.on = false; 
	STD.on_step = -1;
	inh_STDP.on_step = -1;
	inh_STDP.on = false;
	synapse_model = 0; // default model
	
}



void ChemSyn::init(const int syn_type_input, const int i_pre, const int j_post, const int N_pre_input, const int N_post_input, const vector<int> &C_i, const vector<int> &C_j, const vector<double> &K_ij, const vector<double> &D_ij){

	// read parameter
	syn_type = syn_type_input;
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


	
	// parameter-dependent initialisation
	init();

}


void ChemSyn::init(const int syn_type_input, const int j_post, const int N_post_input, const double K_ext, const int Num_ext, const vector<double> &rate_ext_t, const vector<bool> &neurons){

	// Initialise chemical synapses for simulating external neuron population
	syn_type = syn_type_input;
	pop_ind_pre = -1; // -1 for external noisy population
	pop_ind_post = j_post;
	N_pre = 1; // just for initialization
	N_post = N_post_input;
	max_delay_steps = 0; // no delay;


	// Parameters for noise generation
	ext_noise.K_ext = K_ext;
	ext_noise.Num_ext = Num_ext;
	ext_noise.rate_ext_t = rate_ext_t;	
	ext_noise.neurons = neurons;
	
	// Random seed (random engine should be feed with DIFFERENT seed at every implementation)
	random_device rd; // random number from operating system for seed
	my_seed = rd(); // record seed
	//my_seed = 321;
	//cout << "My_seed is: " << my_seed << endl;
	
	// parameter-dependent initialisation
	init();
}



void ChemSyn::init(){
	// parameter-dependent initialisation

	
	// Initialise chemical synapse parameters
	if (syn_type == 0){
		tau_decay = tau_decay_AMPA;
		tau_rise = Dt_trans_AMPA;
	}
	else if (syn_type == 1){	
		tau_decay = tau_decay_GABA;
		tau_rise = Dt_trans_GABA;
	}
	else if (syn_type == 2){
		tau_decay = tau_decay_NMDA;
		tau_rise = Dt_trans_NMDA;
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
	steps_trans = int(round(tau_rise / dt));

	// Initialise exp_step
	exp_step_decay = exp(-dt / tau_decay); // single step
	exp_step_rise = exp(-dt / tau_rise);
	I.assign(N_post, 0);
	gs_sum.assign(N_post, 0);
	
	// transmitter_strength
	K_trans.assign(N_pre, 1.0 / steps_trans); // be careful! 1 / transmitter steps gives zero (int)!!
	
	// Synapse model choice
	// model 0, the default model
	if (synapse_model == 0){ 
		// Initialize pre- and post-synaptic dynamic variables
		gsm_0.buffer_steps = max_delay_steps + steps_trans + 1; // the +1 is vital
		gsm_0.s.assign(N_pre, 0);
		gsm_0.d_gs_sum_buffer.resize(gsm_0.buffer_steps);
		for (int i = 0; i < gsm_0.buffer_steps; ++i){
			gsm_0.d_gs_sum_buffer[i].assign(N_post, 0);
		}
		gsm_0.trans_left.assign(N_pre, 0);
	}
	else if (synapse_model == 1){
		// model 1
		gsm_1.buffer_steps = max_delay_steps + 1; // the +1 is vital
		gsm_1.gs_rise_sum.assign(N_post, 0);
		gsm_1.gs_decay_sum.assign(N_post, 0);
		gsm_1.d_gs_rd_sum_buffer.resize(gsm_1.buffer_steps);
		for (int i = 0; i < gsm_1.buffer_steps; ++i){
			gsm_1.d_gs_rd_sum_buffer[i].assign(N_post, 0);
		}
		// clear model 0
		gsm_0.s.clear();
		gsm_0.d_gs_sum_buffer.clear(); // clear() clears all of its components recursively
		gsm_0.trans_left.clear();
	}
}

const int & ChemSyn::get_syn_type()
{
	return syn_type;
}
const int & ChemSyn::get_pop_ind_pre()
{
	return pop_ind_pre;
}
const int & ChemSyn::get_pop_ind_post()
{
	return pop_ind_post;
}

void ChemSyn::set_synapse_model(const int synapse_model_input){
	if (synapse_model_input != 0){
		synapse_model = synapse_model_input;
		init(); // initialize again
	}
}

	

void ChemSyn::update(const int step_current){

	if (synapse_model == 0){ 
		// short-term depression
		update_STD(step_current);
		
		// inhibitory STDP
		update_inh_STDP(step_current);
		
		// Update transmitter dynamics
		update_gs_sum_model_0(step_current);
	}
	else if (synapse_model == 1){
		update_gs_sum_model_1(step_current);
	}
	
	// Calculate chemical currents
	calc_I();

	// sample data
	sample_data(step_current);
	
	//
	record_stats();
	
}// update



void ChemSyn::update_STD(const int step_current){
	// STD modifies K_trans
	
	if (STD.on_step == step_current){STD.on = true;}
	// short-term depression
	if (STD.on == true){
		for (unsigned int i = 0; i < spikes_pre.size(); ++i){ 
			K_trans[spikes_pre.at(i)] = 1.0 / steps_trans * STD.f_ves[i];
			STD.f_ves[spikes_pre.at(i)] *= 1.0 - STD.p_ves; // decrease at spikes
		}
		for (int i = 0; i < N_pre; ++i){
			STD.f_ves[i] = 1.0 - STD.exp_ves * (1.0 - STD.f_ves[i]); // decay to 1.0
		}
	}

}


void ChemSyn::add_inh_STDP(const int inh_STDP_on_step_input){
	if (syn_type != 1){
		cout << "Warning: initializing inhibitory STDP on non-GABA synapses!" << endl;
	}
	if (synapse_model != 0){
		cout << "Warning: inhibitory STDP is not supported for synapse models other than 0 yet!" << endl;
	}
	
	inh_STDP.on_step = inh_STDP_on_step_input;
	if (inh_STDP.on_step == 0){
		inh_STDP.on = true;
	}
	
	inh_STDP.x_trace_pre.assign(N_pre, 0.0);
	inh_STDP.x_trace_post.assign(N_post, 0.0);
	
	inh_STDP.tau = 20; // ms
	inh_STDP.exp_step = exp(-dt / inh_STDP.tau);
	inh_STDP.eta = 0.0001; // learning rate, 0.0001 is the published value but requires 60min of simulation
	inh_STDP.rho_0 = 0.010; //kHz
	inh_STDP.alpha = 2.0 * inh_STDP.rho_0 * inh_STDP.tau; // depression factor

	// j_2_i and j_2_syn_ind
	inh_STDP.j_2_i.resize(N_post);
	inh_STDP.j_2_syn_ind.resize(N_post);
	int j_post;
	for (int i_pre = 0; i_pre < N_pre; ++i_pre){ 
		for (unsigned int syn_ind = 0; syn_ind < C[i_pre].size(); ++syn_ind){
			j_post = C[i_pre][syn_ind];
			inh_STDP.j_2_i[j_post].push_back( i_pre );
			inh_STDP.j_2_syn_ind[j_post].push_back( int(syn_ind) );
		}
	}
}	

void ChemSyn::calc_I(){
	// need V from post population
	// need gs_sum from previous calculations
	
	if (syn_type == 0){ //AMPA
		for (int j = 0; j < N_post; ++j){
			I[j] = -gs_sum[j] * (V_post.at(j) - V_ex);
		}
	}	
	else if (syn_type == 1){ //GABA
		for (int j = 0; j < N_post; ++j){
			I[j] = -gs_sum[j] * (V_post.at(j) - V_in);
			// For inhibition, every equation is in the same form as excitation. 
			// Only "V_in" encodes its inhibitory nature.
		}
	}
	else if (syn_type == 2){ //NMDA
		for (int j = 0; j < N_post; ++j){
			I[j] = -gs_sum[j] * B[(int)round((V_post.at(j) - B_V_min) / B_dV)] * (V_post.at(j) - V_ex);
		}
	}

	// decay gs_sum
	// numerical error of this integration scheme should be less 1.5%
	for (int j = 0; j < N_post; ++j){ gs_sum[j] *= exp_step_decay; };
	
}

void ChemSyn::update_gs_sum_model_0(const int step_current){
	// See Gu, Yifan, Gong, Pulin, 2016, The dynamics of memory retrieval in hierarchical networks: a modeling study
	// this function updates gs_sum
	if (pop_ind_pre >= 0){
		for (unsigned int i = 0; i < spikes_pre.size(); ++i){ // add spikes (transmitter release)
			gsm_0.trans_left[spikes_pre.at(i)] += steps_trans;
		}
		for (int i_pre = 0; i_pre < N_pre; ++i_pre){
			if (gsm_0.trans_left[i_pre] > 0){
				// conduction delay: put into d_gs_sum_buffer
				int j_post, delay_step, t_ring;
				for (unsigned int syn_ind = 0; syn_ind < C[i_pre].size(); ++syn_ind){ //loop through all the post-synapses
					j_post = C[i_pre][syn_ind]; // index of the post-synaptic neuron
					delay_step = D[i_pre][syn_ind]; // delay in steps for this post-synaptic neuron
					t_ring = int( (step_current + delay_step) % gsm_0.buffer_steps ); // index in the gs_buffer
					gsm_0.d_gs_sum_buffer[t_ring][j_post] += K_trans[i_pre] * (1.0 - gsm_0.s[i_pre]) * K[i_pre][syn_ind];
				}
				gsm_0.trans_left[i_pre] -= 1;
				gsm_0.s[i_pre] += K_trans[i_pre] * (1.0 - gsm_0.s[i_pre]);
			}
		}
		// decay pre-synaptic dynamics
		for (int i = 0; i < N_pre; ++i){ gsm_0.s[i] *= exp_step_decay; };
	}
	else if (pop_ind_pre == -1){ // if external noisy population
		// Contribution of external spikes, assuming square pulse transmitter release
		// Generate current random number generator, note that rate_ext_t is in Hz
		gen.seed(my_seed + step_current);// reseed random engine!!!
		poisson_distribution<int> dist(ext_noise.Num_ext * ext_noise.rate_ext_t[step_current] * (dt / 1000.0));		
		auto ext_spikes = bind(dist, gen);

		// Post-synaptic dynamics
		int t_ring;
		for (int t_trans = 0; t_trans < steps_trans; ++t_trans){
			t_ring = int( (step_current + t_trans) % gsm_0.buffer_steps );
			for (int j_post = 0; j_post < N_post; ++j_post){
				if (ext_noise.neurons[j_post]){
					gsm_0.d_gs_sum_buffer[t_ring][j_post] += K_trans[0] * ext_noise.K_ext * ext_spikes(); 
				}
			}
		}
	}
	// update post-synaptic dynamics
	int t_ring = int( step_current % gsm_0.buffer_steps );
	for (int j_post = 0; j_post < N_post; ++j_post){
		gs_sum[j_post] += gsm_0.d_gs_sum_buffer[t_ring][j_post];
		// should I decay gs_sum here??
	}
	// immediately reset the current buffer to zeros after being used!!
	fill(gsm_0.d_gs_sum_buffer[t_ring].begin(), gsm_0.d_gs_sum_buffer[t_ring].end(), 0.0);
}

void ChemSyn::update_gs_sum_model_1(const int step_current){
	// See Keane, A., Gong, P., 2015, Propagating Waves Can Explain Irregular Neural Dynamics
	// this function updates gs_sum
	for (unsigned int ind = 0; ind < spikes_pre.size(); ++ind){ // loop through all the spikes
		int i_pre = spikes_pre.at(ind);
		int j_post, delay_step, t_ring;
		for (unsigned int syn_ind = 0; syn_ind < C[i_pre].size(); ++syn_ind){ //loop through all the post-synapses
			j_post = C[i_pre][syn_ind]; // index of the post-synaptic neuron
			delay_step = D[i_pre][syn_ind]; // delay in steps for this post-synaptic neuron
			t_ring = int( (step_current + delay_step) % gsm_1.buffer_steps ); // index in the gs_buffer
			gsm_1.d_gs_rd_sum_buffer[t_ring][j_post] += K[i_pre][syn_ind];  // the peak value is linear to the initial impulse. 
		}
	}
	int t_ring = int( step_current % gsm_1.buffer_steps );
	for (int j_post = 0; j_post < N_post; ++j_post){ // Check the error of the following numerical scheme!
		gsm_1.gs_rise_sum[j_post] *= exp_step_rise;
		gsm_1.gs_decay_sum[j_post] *= exp_step_decay;
		gs_sum[j_post] = (gsm_1.gs_decay_sum[j_post] - gsm_1.gs_rise_sum[j_post]) / (tau_decay - tau_rise); 
		gsm_1.gs_rise_sum[j_post] += gsm_1.d_gs_rd_sum_buffer[t_ring][j_post];
		gsm_1.gs_decay_sum[j_post] += gsm_1.d_gs_rd_sum_buffer[t_ring][j_post];
	}
	// immediately reset the current buffer to zeros after being used!!
	fill(gsm_1.d_gs_rd_sum_buffer[t_ring].begin(), gsm_1.d_gs_rd_sum_buffer[t_ring].end(), 0.0);
}


void ChemSyn::add_sampling(const vector<int> & sample_neurons_input, const vector<bool> & sample_time_points_input){
	sample.neurons = sample_neurons_input;
	sample.time_points = sample_time_points_input;
	
	// initialise
	int sample_time_points_tot = 0;// count non zero elements in sample_time_points
	for (unsigned int i = 0; i < sample.time_points.size(); ++i){
		if (sample.time_points[i]){
			sample_time_points_tot += 1;
		}
	}
	int sample_neurons_tot = sample.neurons.size();// count non zero elements in sample_time_points

	sample.data.resize(sample_neurons_tot);
	for (int i = 0; i < sample_neurons_tot; ++i){
		sample.data[i].reserve(sample_time_points_tot); // reserve and push_back so that it won't be affected by adapting step_tot
	}

}

void ChemSyn::add_short_term_depression(const int STD_on_step_input){
	if (syn_type != 0){
		cout << "Warning: initializing STD on non-AMPA synapses!" << endl;
	}
	if (synapse_model != 0){
		cout << "Warning: STD is not supported for synapse models other than 0 yet!" << endl;
	}
	
	STD.on_step = STD_on_step_input;
	if (STD.on_step == 0){
		STD.on = true;
	}
	// short term depression
	STD.p_ves =  0.4; // see X. Wang, 1999, The Journal of Neuroscience
	STD.tau_ves =  700; // ms
	STD.f_ves.assign(N_pre, 1.0);
	STD.exp_ves = exp(-dt / STD.tau_ves);
}




void ChemSyn::update_inh_STDP(const int step_current){
	// inh_STDP modifies K
	
	if (inh_STDP.on_step == step_current){inh_STDP.on = true;}
	if (inh_STDP.on == true){
		// update x_trace
		for (unsigned int i = 0; i < spikes_pre.size(); ++i){
			inh_STDP.x_trace_pre[spikes_pre.at(i)] += 1.0;
		}
		for (unsigned int j = 0; j < spikes_post.size(); ++j){
			inh_STDP.x_trace_post[spikes_post.at(j)] += 1.0;
		}
		// update K
		int i_pre, j_post;
		for (unsigned int ind_spike = 0; ind_spike < spikes_pre.size(); ++ind_spike){
			i_pre = spikes_pre.at(ind_spike);
			for (unsigned int syn_ind = 0; syn_ind < C[i_pre].size(); ++syn_ind){
				j_post = C[i_pre][syn_ind];
				K[i_pre][syn_ind] += inh_STDP.eta * ( inh_STDP.x_trace_post[j_post] - inh_STDP.alpha );
			}
		}
		int syn_ind;
		for (unsigned int ind_spike = 0; ind_spike < spikes_post.size(); ++ind_spike){
			j_post = spikes_post.at(ind_spike);
			for (unsigned int ind = 0; ind < inh_STDP.j_2_i[j_post].size(); ++ind){
				// j_2_i and j_2_syn_ind together serve as the "inverse function" of j_post = C[i_pre][syn_ind]
				i_pre = inh_STDP.j_2_i[j_post][ind];
				syn_ind = inh_STDP.j_2_syn_ind[j_post][ind];
				K[i_pre][syn_ind] += inh_STDP.eta * inh_STDP.x_trace_pre[i_pre];
			}
		}
		// for testing
		//if (tmp_data.size() == 0){
		//	tmp_data.resize(2);
		//}
		//tmp_data[0].push_back(K[0][0]);
		//tmp_data[1].push_back(K[100][0]);
		for (int i = 0; i < N_pre; ++i){ inh_STDP.x_trace_pre[i] *= inh_STDP.exp_step; }
		for (int j = 0; j < N_post; ++j){ inh_STDP.x_trace_post[j] *= inh_STDP.exp_step; }
	}
}


void ChemSyn::sample_data(const int step_current){
	if (!sample.neurons.empty()){
		if (sample.time_points[step_current]){ // push_back is amazing
			for (unsigned int i = 0; i < sample.neurons.size(); ++i){ // performance issue when sampling many neurons?
				int ind_temp = sample.neurons[i];
				sample.data[i].push_back( I[ind_temp] );
			}
		}
	}
}



void ChemSyn::set_para(string para_str){
	if (!para_str.empty()){
		istringstream para(para_str);
		string para_name, para_value_str; 
		double para_value;
		while (getline(para, para_name, delim)){
			getline(para, para_value_str, delim); // get parameter value (assume double)
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


string ChemSyn::dump_para(){
	stringstream dump;


	dump << "pop_ind_pre" << delim << pop_ind_pre << delim << endl;
	dump << "pop_ind_post" << delim << pop_ind_post << delim << endl;
	dump << "syn_type" << delim << syn_type << delim << endl;

	dump << "V_ex" << delim << V_ex << delim << endl;
	dump << "V_in" << delim << V_in << delim << endl;

	dump << "seed" << delim << my_seed << delim << endl;

	dump << "synapse_model" << delim << synapse_model << endl;
	
	if (syn_type == 0){
		dump << "Dt_trans_AMPA" << delim << Dt_trans_AMPA << delim << endl;
		dump << "tau_decay_AMPA" << delim << tau_decay_AMPA << delim << endl;
	}
	else if (syn_type == 1){
		dump << "Dt_trans_GABA" << delim << Dt_trans_GABA << delim << endl;
		dump << "tau_decay_GABA" << delim << tau_decay_GABA << delim << endl;
	}
	else if (syn_type == 2){
		dump << "Dt_trans_NMDA" << delim << Dt_trans_NMDA << delim << endl;
		dump << "tau_decay_NMDA" << delim << tau_decay_NMDA << delim << endl;
	}

		
	return dump.str();
}

void ChemSyn::start_stats_record(){
	stats.record = true;
	stats.I_mean.reserve(step_tot);
	stats.I_std.reserve(step_tot);
}

void ChemSyn::output_results(ofstream& output_file){
	// SYND001 # synapse parameters
	// count number of variables
	stringstream dump_count;
	string para_str = dump_para();
	dump_count << para_str;
	string str_temp;
	int var_number = 0;
	while(getline(dump_count,str_temp)){++var_number;} // count number of variables
	output_file << indicator << " SYND001" << endl;
	output_file << var_number << delim << endl;
	output_file << para_str;
	
	
	// SYND002 # sampled synapse data
	if (!sample.neurons.empty()){
		output_file << indicator << " SYND002" << endl;
		output_file << pop_ind_pre << delim << pop_ind_post << delim << syn_type << delim << sample.neurons.size() << delim << endl;
		write2file(output_file, sample.data); // 2D matrix
	}

	
	// SYND003 # currents mean and std
	if (stats.record){
		output_file << indicator << " SYND003" << endl;
		output_file << pop_ind_pre << delim << pop_ind_post << delim << syn_type << delim << endl;
		write2file(output_file, stats.I_mean);
		write2file(output_file, stats.I_std);
	}
	
	// tmp data
	if (tmp_data.size() != 0){
		output_file << indicator << " SYND004" << endl;
		output_file << pop_ind_pre << delim << pop_ind_post << delim << syn_type << delim << tmp_data.size() << delim << endl;
		write2file(output_file, tmp_data);	
	}
}



void ChemSyn::recv_pop_data(vector<NeuroPop*> &NeuronPopArray){
	// get current spikes from pre-pop
	if (pop_ind_pre >= 0){
		spikes_pre = NeuronPopArray[pop_ind_pre]->get_spikes_current(); // This might be problematic!!!
		spikes_post = NeuronPopArray[pop_ind_post]->get_spikes_current();
	}
	// get current V from post-pop
	V_post = NeuronPopArray[pop_ind_post]->get_V(); // This  might be  problematic!!!
}


void ChemSyn::send_pop_data(vector<NeuroPop*> &NeuronPopArray){
	
	NeuronPopArray[pop_ind_post]->recv_I(I, pop_ind_pre, syn_type);

}

void ChemSyn::record_stats(){
	if (stats.record){
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
		stats.I_mean.push_back(mean_tmp);
		stats.I_std.push_back(std_tmp);
	}
}


#ifdef HDF5

void ChemSyn::import_restart(H5File& file, int syn_ind){

	string str;

	string syn_str = "/syns/syn" + to_string(syn_ind)+"/";

	dt=read_scalar_HDF5<double>(file, syn_str+"dt");
	step_tot=read_scalar_HDF5<int>(file,syn_str+"step_tot");
	pop_ind_pre=read_scalar_HDF5<int>(file,syn_str+"pop_ind_pre");
	pop_ind_post=read_scalar_HDF5<int>(file,syn_str+"pop_ind_post");
	N_pre=read_scalar_HDF5<int>(file, syn_str+"N_pre");
	N_post=read_scalar_HDF5<int>(file, syn_str+"N_post");
	syn_type=read_scalar_HDF5<int>(file, syn_str+"syn_type");
	V_ex=read_scalar_HDF5<double>(file, syn_str+"V_ex");
	V_in=read_scalar_HDF5<double>(file, syn_str+"V_in");
	max_delay_steps=read_scalar_HDF5<int>(file, syn_str+"max_delay_steps");
	read_vector_HDF5(file, syn_str+"V_post",V_post);
	read_vector_HDF5(file, syn_str+"spikes_pre",spikes_pre);
	read_vector_HDF5(file,syn_str+"spikes_post",spikes_post);
	read_vector_HDF5(file, syn_str+"I",I);

		

	str = syn_str+"/Stats/";
	if(group_exist_HDF5(file,str)){
		stats.record=read_scalar_HDF5<bool>(file, str+"record");
		read_vector_HDF5(file, str+"I_mean",stats.I_mean);
		read_vector_HDF5(file, str+"I_std",stats.I_std);
		start_stats_record();
	}

	str =syn_str+"/Sample/";
	if(group_exist_HDF5(file,str)){
		read_vector_HDF5(file, str+"neurons", sample.neurons);
		read_vector_HDF5(file,str+ "time_points", sample.time_points);
		add_sampling(sample.neurons, sample.time_points);
		// read_matrix_HDF5(file, str+ "data",sample.data);
	}

	Dt_trans_AMPA=read_scalar_HDF5<double>(file, syn_str+"Dt_trans_AMPA");
	Dt_trans_GABA=read_scalar_HDF5<double>(file, syn_str+"Dt_trans_GABA");
	Dt_trans_NMDA=read_scalar_HDF5<double>(file, syn_str+"Dt_trans_NMDA");
	tau_decay_AMPA=read_scalar_HDF5<double>(file, syn_str+"tau_decay_AMPA");
	tau_decay_GABA=read_scalar_HDF5<double>(file, syn_str+"tau_decay_GABA");
	tau_decay_NMDA=read_scalar_HDF5<double>(file,syn_str+ "tau_decay_NMDA");
	tau_rise=read_scalar_HDF5<double>(file, syn_str+"tau_rise");
	tau_decay=read_scalar_HDF5<double>(file,syn_str+ "tau_decay");
	steps_trans=read_scalar_HDF5<double>(file,syn_str+ "steps_trans");
	read_vector_HDF5(file, syn_str+"K_trans",K_trans);
	exp_step_decay=read_scalar_HDF5<double>(file, syn_str+"exp_step_decay");
	exp_step_rise=read_scalar_HDF5<double>(file,syn_str+ "exp_step_rise");
	miuMg_NMDA=read_scalar_HDF5<double>(file, syn_str+"miuMg_NMDA");
	gamma_NMDA=read_scalar_HDF5<double>(file, syn_str+"gamma_NMDA");
	B_V_min=read_scalar_HDF5<double>(file, syn_str+"B_V_min");
	B_V_max=read_scalar_HDF5<double>(file, syn_str+"B_V_max");
	B_dV=read_scalar_HDF5<double>(file,syn_str+ "B_dV");
	read_vector_HDF5(file, syn_str+"B",B);

	str =syn_str+"/Inh_STDP/";
	if(group_exist_HDF5(file,str)){
		inh_STDP.on=read_scalar_HDF5<bool>(file, str+ "on");
		read_vector_HDF5(file, str+ "x_trace_pre",inh_STDP.x_trace_pre);
		read_vector_HDF5(file,  str+"x_trace_post",inh_STDP.x_trace_post);
		inh_STDP.tau=read_scalar_HDF5<double>(file,  str+"tau");
		inh_STDP.exp_step=read_scalar_HDF5<double>(file, str+ "exp_step");
		inh_STDP.eta=read_scalar_HDF5<double>(file, str+ "eta");
		inh_STDP.rho_0=read_scalar_HDF5<double>(file, str+ "rho_0");
		inh_STDP.alpha=read_scalar_HDF5<double>(file,  str+"alpha");
		inh_STDP.on_step=read_scalar_HDF5<int>(file,  str+"on_step");
		read_matrix_HDF5(file,  str+"j_2_i",inh_STDP.j_2_i);
		read_matrix_HDF5(file, str+ "j_2_syn_ind",inh_STDP.j_2_syn_ind);
	}

	str =syn_str+"/Std/";
	if(group_exist_HDF5(file,str)){
		STD.on=read_scalar_HDF5<bool>(file, str+ "on");
		STD.p_ves=read_scalar_HDF5<double>(file,  str+"p_ves");
		STD.tau_ves=read_scalar_HDF5<double>(file,  str+"tau_ves");
		STD.exp_ves=read_scalar_HDF5<double>(file,  str+"exp_ves");
		STD.on_step=read_scalar_HDF5<bool>(file,  str+"on_step");
		read_vector_HDF5(file, str+ "f_ves",STD.f_ves);
	}

	synapse_model=read_scalar_HDF5<int>(file,syn_str+ "synapse_model");
	read_vector_HDF5(file, syn_str+"gs_sum",gs_sum);

	str =syn_str+"/Gsm_0/";		
	if(group_exist_HDF5(file,str)){
		gsm_0.buffer_steps=read_scalar_HDF5<int>(file,  str+"buffer_steps");
		read_vector_HDF5(file,  str+"s",gsm_0.s);
		read_vector_HDF5(file,  str+"trans_left",gsm_0.trans_left);
		read_matrix_HDF5(file,  str+"d_gs_sum_buffer",gsm_0.d_gs_sum_buffer);
	}

	str =syn_str+"/Gsm_1/";
	if(group_exist_HDF5(file,str)){
		gsm_1.buffer_steps=read_scalar_HDF5<int>(file,  str+"buffer_steps");
		read_vector_HDF5(file, str+ "gs_rise_sum",gsm_1.gs_rise_sum);
		read_vector_HDF5(file,  str+"gs_decay_sum",gsm_1.gs_decay_sum);
		read_matrix_HDF5(file,  str+"d_gs_rd_sum_buffer",gsm_1.d_gs_rd_sum_buffer);
	}

	read_matrix_HDF5(file,syn_str+ "C",C);
	read_matrix_HDF5(file, syn_str+"D",D);
	read_matrix_HDF5(file,syn_str+"K",K);

	str =syn_str+"/Ext_noise/";
	if(group_exist_HDF5(file,str)){
		ext_noise.K_ext=read_scalar_HDF5<double>(file, str+ "K_ext");
		ext_noise.Num_ext=read_scalar_HDF5<int>(file,  str+"Num_ext");
		read_vector_HDF5(file,  str+"neurons", ext_noise.neurons);
		read_vector_HDF5(file,  str+"rate_ext_t", ext_noise.rate_ext_t);
	}
	my_seed=read_scalar_HDF5<int>(file, syn_str+"my_seed");
	// +++ TODO BASE_GENERATOR_TYP
	
	// JH Learning
	str = syn_str+"/JH_Learn/";
	if(group_exist_HDF5(file,str)){		
		jh_learn_syn.on=read_scalar_HDF5<bool>(file,str+ "on");
		read_matrix_HDF5(file,str+ "post_t_hist",jh_learn_syn.post_t_hist);
		read_vector_HDF5(file,str+ "ind_post_new",jh_learn_syn.ind_post_new);
		read_vector_HDF5(file,str+ "ind_post_old",jh_learn_syn.ind_post_old);
		read_matrix_HDF5(file,str+ "pre_t_hist",jh_learn_syn.pre_t_hist);
		read_vector_HDF5(file,str+ "ind_pre_new",jh_learn_syn.ind_pre_new);
		read_vector_HDF5(file,str+ "ind_pre_old",jh_learn_syn.ind_pre_old);
		read_vector_HDF5(file,str+ "old_pre",jh_learn_syn.old_pre);
		jh_learn_syn.post_hist_len=read_scalar_HDF5<int>(file,str+ "post_hist_len");
		jh_learn_syn.pre_hist_len=read_scalar_HDF5<int>(file,str+ "pre_hist_len");
		jh_learn_syn.ntype_pre=read_scalar_HDF5<int>(file,str+ "ntype_pre");
		jh_learn_syn.ntype_post=read_scalar_HDF5<int>(file,str+ "ntype_post");
		read_vector_HDF5(file,str+ "Vint",jh_learn_syn.Vint);
		read_vector_HDF5(file,str+ "Vint_ctr",jh_learn_syn.Vint_ctr);
		read_vector_HDF5(file,str+ "Q_pre",jh_learn_syn.Q_pre);
		read_matrix_HDF5(file,str+ "post_V_hist",jh_learn_syn.post_V_hist);
		read_matrix_HDF5(file,str+ "post_R_hist",jh_learn_syn.post_R_hist);
		jh_learn_syn.t_ind=read_scalar_HDF5<int>(file,str+ "t_ind");
		jh_learn_syn.inf_steps=read_scalar_HDF5<int>(file,str+ "inf_steps");
		jh_learn_syn.inf_scaleI=read_scalar_HDF5<double>(file,str+ "inf_scaleI");
		jh_learn_syn.inf_scaleE=read_scalar_HDF5<double>(file,str+ "inf_scaleE");
		jh_learn_syn.learn_rate_E=read_scalar_HDF5<double>(file,str+ "learn_rate_E");
		jh_learn_syn.learn_rate_all_E=read_scalar_HDF5<double>(file,str+ "learn_rate_all_E");
		jh_learn_syn.learn_rate_I=read_scalar_HDF5<double>(file,str+ "learn_rate_I");
		jh_learn_syn.learn_rate_all_I=read_scalar_HDF5<double>(file,str+ "learn_rate_all_I");
		jh_learn_syn.tau=read_scalar_HDF5<double>(file,str+ "tau");
		jh_learn_syn.C=read_scalar_HDF5<double>(file, str+"C");
		jh_learn_syn.noise=read_scalar_HDF5<double>(file,str+ "noise");
		read_matrix_HDF5(file,str+"j_2_i",jh_learn_syn.j_2_i);
		read_matrix_HDF5(file,str+"j_2_syn_ind",jh_learn_syn.j_2_syn_ind);
	}
}

void ChemSyn::export_restart(Group& group, int syn_ind){
	string syn_str = "/syns/syn" + to_string(syn_ind);
	Group group_syn = group.createGroup(syn_str);

	write_scalar_HDF5(group_syn,dt, "dt");
	write_scalar_HDF5(group_syn,step_tot, "step_tot");
	write_scalar_HDF5(group_syn,pop_ind_pre,"pop_ind_pre");
	write_scalar_HDF5(group_syn,pop_ind_post,"pop_ind_post");
	write_scalar_HDF5(group_syn,N_pre, "N_pre");
	write_scalar_HDF5(group_syn,N_post, "N_post");
	write_scalar_HDF5(group_syn,syn_type, "syn_type");
	write_scalar_HDF5(group_syn,V_ex, "V_ex");
	write_scalar_HDF5(group_syn,V_in, "V_in");
	write_scalar_HDF5(group_syn,max_delay_steps, "max_delay_steps");
	write_vector_HDF5(group_syn,V_post, "V_post");
	write_vector_HDF5(group_syn,spikes_pre, "spikes_pre");
	write_vector_HDF5(group_syn,spikes_post, "spikes_post");
	write_vector_HDF5(group_syn,I, "I");

		
	if(stats.record){
		string str = syn_str+"/Stats";
		Group group_stats = group_syn.createGroup(str);
		write_scalar_HDF5(group_stats,stats.record, "record");
		// write_vector_HDF5(group_stats,stats.I_mean, "I_mean");
		// write_vector_HDF5(group_stats,stats.I_std, "I_std");
	}

	if(!sample.time_points.empty()){
		string str =syn_str+"/Sample";
		Group group_sample = group_syn.createGroup(str);
		write_vector_HDF5(group_sample, sample.neurons, "neurons");
		write_vector_HDF5(group_sample, sample.time_points, "time_points");
		// write_matrix_HDF5(group_sample, sample.data, "data");
	}

	write_scalar_HDF5(group_syn,Dt_trans_AMPA, "Dt_trans_AMPA");
	write_scalar_HDF5(group_syn,Dt_trans_GABA, "Dt_trans_GABA");
	write_scalar_HDF5(group_syn,Dt_trans_NMDA, "Dt_trans_NMDA");
	write_scalar_HDF5(group_syn,tau_decay_AMPA, "tau_decay_AMPA");
	write_scalar_HDF5(group_syn,tau_decay_GABA, "tau_decay_GABA");
	write_scalar_HDF5(group_syn,tau_decay_NMDA, "tau_decay_NMDA");
	write_scalar_HDF5(group_syn,tau_rise, "tau_rise");
	write_scalar_HDF5(group_syn,tau_decay, "tau_decay");
	write_scalar_HDF5(group_syn,steps_trans, "steps_trans");
	write_vector_HDF5(group_syn,K_trans, "K_trans");
	write_scalar_HDF5(group_syn,exp_step_decay, "exp_step_decay");
	write_scalar_HDF5(group_syn,exp_step_rise, "exp_step_rise");
	write_scalar_HDF5(group_syn,miuMg_NMDA, "miuMg_NMDA");
	write_scalar_HDF5(group_syn,gamma_NMDA, "gamma_NMDA");
	write_scalar_HDF5(group_syn,B_V_min, "B_V_min");
	write_scalar_HDF5(group_syn,B_V_max, "B_V_max");
	write_scalar_HDF5(group_syn,B_dV, "B_dV");
	write_vector_HDF5(group_syn,B, "B");

	if(inh_STDP.on){
		string str =syn_str+"/Inh_STDP";
		Group group_Inh_STDP = group_syn.createGroup(str);
		write_scalar_HDF5(group_Inh_STDP,inh_STDP.on, "on");
		write_vector_HDF5(group_Inh_STDP,inh_STDP.x_trace_pre, "x_trace_pre");
		write_vector_HDF5(group_Inh_STDP,inh_STDP.x_trace_post, "x_trace_post");
		write_scalar_HDF5(group_Inh_STDP,inh_STDP.tau, "tau");
		write_scalar_HDF5(group_Inh_STDP,inh_STDP.exp_step, "exp_step");
		write_scalar_HDF5(group_Inh_STDP,inh_STDP.eta, "eta");
		write_scalar_HDF5(group_Inh_STDP,inh_STDP.rho_0, "rho_0");
		write_scalar_HDF5(group_Inh_STDP,inh_STDP.alpha, "alpha");
		write_scalar_HDF5(group_Inh_STDP,inh_STDP.on_step, "on_step");
		write_matrix_HDF5(group_Inh_STDP,inh_STDP.j_2_i, "j_2_i");
		write_matrix_HDF5(group_Inh_STDP,inh_STDP.j_2_syn_ind, "j_2_syn_ind");
	}

	if(STD.on){
		string str =syn_str+"/Std";
		Group group_STD = group_syn.createGroup(str);
		write_scalar_HDF5(group_STD,STD.on, "on");
		write_scalar_HDF5(group_STD,STD.p_ves, "p_ves");
		write_scalar_HDF5(group_STD,STD.tau_ves, "tau_ves");
		write_scalar_HDF5(group_STD,STD.exp_ves, "exp_ves");
		write_scalar_HDF5(group_STD,STD.on_step, "on_step");
		write_vector_HDF5(group_STD,STD.f_ves, "f_ves");
	}

	write_scalar_HDF5(group_syn,synapse_model, "synapse_model");
	write_vector_HDF5(group_syn,gs_sum, "gs_sum");
		
	if(synapse_model==0){
		string str =syn_str+"/Gsm_0";
		Group group_Gsm_0 = group_syn.createGroup(str);
		write_scalar_HDF5(group_Gsm_0,gsm_0.buffer_steps, "buffer_steps");
		write_vector_HDF5(group_Gsm_0,gsm_0.s, "s");
		write_vector_HDF5(group_Gsm_0,gsm_0.trans_left, "trans_left");
		write_matrix_HDF5(group_Gsm_0,gsm_0.d_gs_sum_buffer, "d_gs_sum_buffer");
	}
	
	if(synapse_model==1){
		string str =syn_str+"/Gsm_1";
		Group group_Gsm_1 = group_syn.createGroup(str);
		write_scalar_HDF5(group_Gsm_1,gsm_1.buffer_steps, "buffer_steps");
		write_vector_HDF5(group_Gsm_1,gsm_1.gs_rise_sum, "gs_rise_sum");
		write_vector_HDF5(group_Gsm_1,gsm_1.gs_decay_sum, "gs_decay_sum");
		write_matrix_HDF5(group_Gsm_1,gsm_1.d_gs_rd_sum_buffer, "d_gs_rd_sum_buffer");
	}

	write_matrix_HDF5(group_syn,C, "C");
	write_matrix_HDF5(group_syn,D, "D");
	write_matrix_HDF5(group_syn,K, "K");

	if(!ext_noise.neurons.empty()){
		string str =syn_str+"/Ext_noise";
		Group group_Ext_noise = group_syn.createGroup(str);
		write_scalar_HDF5(group_Ext_noise, ext_noise.K_ext, "K_ext");
		write_scalar_HDF5(group_Ext_noise, ext_noise.Num_ext, "Num_ext");
		write_vector_HDF5(group_Ext_noise, ext_noise.neurons, "neurons");
		
		write_vector_HDF5(group_Ext_noise, ext_noise.rate_ext_t, "rate_ext_t");
	}
	write_scalar_HDF5(group_syn, my_seed, "my_seed");
	// +++ TODO BASE_GENERATOR_TYP

	// JH Learning
	if (jh_learn_syn.on){
		string str = syn_str+"/JH_Learn/";
		Group group_JH_Learn = group_syn.createGroup(str);
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.on, "on");
		write_matrix_HDF5(group_JH_Learn,jh_learn_syn.post_t_hist, "post_t_hist");
		write_vector_HDF5(group_JH_Learn,jh_learn_syn.ind_post_new, "ind_post_new");
		write_vector_HDF5(group_JH_Learn,jh_learn_syn.ind_post_old, "ind_post_old");
		write_matrix_HDF5(group_JH_Learn,jh_learn_syn.pre_t_hist, "pre_t_hist");
		write_vector_HDF5(group_JH_Learn,jh_learn_syn.ind_pre_new, "ind_pre_new");
		write_vector_HDF5(group_JH_Learn,jh_learn_syn.ind_pre_old, "ind_pre_old");
		write_vector_HDF5(group_JH_Learn,jh_learn_syn.old_pre, "old_pre");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.post_hist_len, "post_hist_len");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.pre_hist_len, "pre_hist_len");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.ntype_pre, "ntype_pre");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.ntype_post, "ntype_post");
		write_vector_HDF5(group_JH_Learn,jh_learn_syn.Vint, "Vint");
		write_vector_HDF5(group_JH_Learn,jh_learn_syn.Vint_ctr, "Vint_ctr");
		write_vector_HDF5(group_JH_Learn,jh_learn_syn.Q_pre, "Q_pre");
		write_matrix_HDF5(group_JH_Learn,jh_learn_syn.post_V_hist, "post_V_hist");
		write_matrix_HDF5(group_JH_Learn,jh_learn_syn.post_R_hist, "post_R_hist");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.t_ind, "t_ind");

		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.inf_steps, "inf_steps");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.inf_scaleI, "inf_scaleI");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.inf_scaleE, "inf_scaleE");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.learn_rate_E, "learn_rate_E");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.learn_rate_all_E, "learn_rate_all_E");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.learn_rate_I, "learn_rate_I");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.learn_rate_all_I, "learn_rate_all_I");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.tau, "tau");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.C, "C");
		write_scalar_HDF5(group_JH_Learn,jh_learn_syn.noise, "noise");
		write_matrix_HDF5(group_JH_Learn,jh_learn_syn.j_2_i,"j_2_i");
		write_matrix_HDF5(group_JH_Learn,jh_learn_syn.j_2_syn_ind ,"j_2_syn_ind");
	}
}
void ChemSyn::output_results(H5File& file, int syn_ind){
	// new group
	stringstream group_name;
	group_name << "/syn_result_"  << syn_ind;
	Group group_syn = file.createGroup(group_name.str());
	
	write_string_HDF5(group_syn, dump_para(), string("syn_para"));
		
	if (!sample.neurons.empty()){
		write_matrix_HDF5(group_syn, sample.data, string("sample_data"));
	}
	
	if (stats.record){
		write_vector_HDF5(group_syn, stats.I_mean, string("stats_I_mean"));
		write_vector_HDF5(group_syn, stats.I_std, string("stats_I_std"));
	}
}

#endif
