function main_test_HDF5(varargin)

% clc;clear;close all;
% addpath(genpath(cd));
% cd ~/tmp_data
    
loop_num = 2;

% seed the matlab rand function! The seed is global.
[FID ] = new_ygin_files_and_randseedHDF5(loop_num);


k = 1*2.4e-3; % miuSiemens


% Basic parameters
dt = 0.1;
step_tot = 10000;
N = [10; 12; ];
writeBasicParaHDF5(FID, dt, step_tot, N)
Num_pop = length(N);
discard_transient = 0; % ms

% write pop para
writePopParaHDF5(FID, 1,  'tau_ref', 3.1);
writePopParaHDF5(FID, 2,  'tau_ref', 3.2);
% write synapse para
writeSynParaHDF5(FID, 'tau_decay_AMPA', 2.0, 'Dt_trans_AMPA', 0.5);

% external current settings (int pop_ind, double mean, double std)
% I_ext_strength = 10; %1.4; % nA
% writeExtCurrentSettingsHDF5(FID, 1, I_ext_strength*ones(1,N(1)), 0*ones(1,N(1)));
% 
% g_ext_strength = 0.5;
% writeExtConductanceSettingsHDF5(FID, 2, g_ext_strength*ones(1,N(2)), 0*ones(1,N(2)));

% external spike settings (neuron-invariant)
writeExtSpikeSettingsHDF5(FID, 2, 1, k,  20, 10*ones(1,step_tot), ones(1, N(1)) );

% external spike settings (time-invariant)
writeExtSpikeTinvSettingsHDF5(FID, 1, 1, k,  20, 10*randn(1,N(1)))
 
% neuronal data sampling
sample_steps = zeros(1,step_tot);
sample_steps(1:1:step_tot) =  true;
writeNeuronSamplingHDF5(FID, 1, ones(1,8), 1:2:10 , sample_steps);
writeNeuronSamplingHDF5(FID, 2, ones(1,8), 1 , sample_steps);


%%%%%%%%%%%%%%%%%%% Chemical Connections %%%%%%%%%%%%%%%%%%%%%%%
% type(1:AMAP, 2:GABAa, 3:NMDA)

pop_type = [1 2];
for i_pre = 1:Num_pop
    for j_post = 1:Num_pop
        [I, J, ~] = find(MyRandomGraphGenerator('E_R_pre_post','N_pre',N(i_pre),'N_post',N(j_post),'p',rand));
        K = ones(size(I))*k;
        D = ones(size(I))*1;
        writeChemicalConnectionHDF5(FID,  pop_type(i_pre), i_pre,  j_post,  I,J,K,D); % (FID, type, i_pre, j_post, I, J, K, D)
    end
end



% synapse data sampling
writeSynSamplingHDF5(FID, 1,  2, 1,  1, sample_steps)

writeInhSTDPHDF5(FID, 2, 1, 10)

writeSTDHDF5(FID, 1, 1, 10)

% Add LFP sampling
writeLFPRecordHDF5(FID, 1, ones(3, 10));


% V mean and std record
writePopStatsRecordHDF5(FID, 1);

% I mean and std record
writeSynStatsRecordHDF5(FID, 1,  2, 1);

%%%% Define the initial condition
p_fire = [0 0]; % initial firing probabilities for both populations
% set initial V distribution to be [V_rt, V_rt + (V_th-V_rt)*r_V0] 
r_V0 = [0 0];
writeInitCondHDF5(FID, r_V0, p_fire)


writeSpikeFreqAdptHDF5(FID, 1);

% %%%%%%% write runaway killer
min_ms = 10*1000; % 10 sec
runaway_Hz = 20; % ??
Hz_ms = 1000; % ms
run_away_pop = 1;
writeRunawayKillerHDF5(FID, run_away_pop, min_ms, runaway_Hz, Hz_ms);
% %%%%%%%%%%%%%%%%%%%%%%%
                
                

% % Explanatory (ExplVar) and response variables (RespVar) for cross-simulation data gathering and post-processing
% % Record explanatory variables, also called "controlled variables"
% writeExplVar(FID, 'discard_transient', discard_transient, 'k', k);
% % writeExplVar(FID, 'I_ext_strength', I_ext_strength);
% writeExplVar(FID, 'comments', 'calibration');
    
% append this file self into .ygin for future reference
% appendThisMatlabFile(FID)


end


