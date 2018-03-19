function main_test_single_EPSP(varargin)

% clc;clear;close all;
% addpath(genpath(cd));
%cd ~/tmp_data
    
loop_num = 1;


% seed the matlab rand function! The seed is global.
[FID, FID_syn] = new_ygin_files_and_randseed(loop_num);


k = 2.4*10^-3*0.6; % miuSiemens


% Basic parameters
dt = 0.1;
step_tot = 1000;
N = [1; 1];
writeBasicPara(FID, dt, step_tot, N)
Num_pop = length(N);
discard_transient = 50; % ms

% write pop para
writePopPara(FID, 1,  'tau_ref', 2);
writePopPara(FID, 2,  'tau_ref', 2);
% write synapse para
writeSynPara(FID, 'tau_decay_GABA', 3);

% external current settings (int pop_ind, double mean, double std)
I_ext_strength = 5; %1.4; % nA
writeExtCurrentSettings(FID, 1, I_ext_strength, 0);


% neuronal data sampling
sample_steps = zeros(1,step_tot);
sample_steps(1:end) =  true;
writeNeuronSampling(FID, 2, ones(1,7), 1, sample_steps);
writeNeuronSampling(FID, 1, ones(1,7), 1, sample_steps);

% % initial firing rate
% writeInitV(FID,[0.2,0.2]);

       
%%%%%%%%%%%%%%%%%%% Chemical Connections %%%%%%%%%%%%%%%%%%%%%%%
% type(1:AMAP, 2:GABAa, 3:NMDA)

[I, J, ~] = find(MyRandomGraphGenerator('E_R_pre_post','N_pre',N(1),'N_post',N(2),'p',1));
K = ones(size(I))*k;
D = ones(size(I))*1;
writeChemicalConnection(FID_syn, 1,  1, 2,   I,J,K,D); % (FID, type, i_pre, j_post, I, J, K, D)

% Explanatory (ExplVar) and response variables (RespVar) for cross-simulation data gathering and post-processing
% Record explanatory variables, also called "controlled variables"
writeExplVar(FID, 'discard_transient', discard_transient, 'k', k);
% writeExplVar(FID, 'I_ext_strength', I_ext_strength);
writeExplVar(FID, 'comments', 'calibration');
    
% append this file self into .ygin for future reference
appendThisMatlabFile(FID)


end



% This function must be here!
function appendThisMatlabFile(FID)
breaker = repmat('#',1,80);
fprintf(FID, '%s\n', breaker);
fprintf(FID, '%s\n', '# MATLAB script generating this file: ');
fprintf(FID, '%s\n', breaker);
Fself = fopen([mfilename('fullpath'),'.m'],'r');
while ~feof(Fself)
    tline = fgetl(Fself);
    fprintf(FID, '%s\n', tline);
end
fprintf(FID, '%s\n', breaker);
fprintf(FID, '%s\n', breaker);
fprintf(FID, '%s\n', breaker);
end

