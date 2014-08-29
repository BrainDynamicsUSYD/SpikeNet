function main_test_code(varargin)

clc;clear all;close all;
cd /import/yossarian1/yifan/Project1/
addpath(genpath(cd));
cd tmp_data
    

loop_num = 1;

% Creat ygin file
% using loop_num in filename to ensure unique naming!
% Otherwise overwriting may occur when using PBS.
name = [sprintf('%03g-', loop_num), datestr(now,'yyyymmdd-HHMM-SSFFF'), '.ygin'];
FID = fopen(name, 'w'); % creat file
FID_syn = fopen([name,'_syn'], 'w'); % creat file
%%%%%%%%%%%%%%%%%%%%%%%%

k = 1*2.4e-3; % miuSiemens


% Basic parameters
dt = 0.1;
step_tot = 500;
N = [2; 6];
writeBasicPara(FID, dt, step_tot, N)
Num_pop = length(N);

% write pop para
writePopPara(FID, 1,  'tau_ref', 3.1);
writePopPara(FID, 2,  'tau_ref', 3.2);
% write synapse para
writeSynPara(FID, 'tau_decay_AMPA', 3.3);

% external current settings (int pop_ind, double mean, double std)
I_ext_strength = 1.4; % nA
writeExtCurrentSettings(FID, 1, I_ext_strength, 0);

% external spike settings
writeExtSpikeSettings(FID, 1, 1, k,  20, 10*ones(1,step_tot));

% neuronal data sampling
writeNeuronSampling(FID, 1, ones(1,7), 1:2);
writeNeuronSampling(FID, 2, ones(1,7), 1:6);

% populational data sampling
writePopSampling(FID, 1, ones(1,step_tot));
writePopSampling(FID, 2, ones(1,step_tot));

% initial firing rate
writeInitV(FID,[0.2,0.2]);

%%%%%%%%%%%%%%%%%%% Chemical Connections %%%%%%%%%%%%%%%%%%%%%%%
% type(1:AMAP, 2:GABAa, 3:NMDA)

[I, J, ~] = find(MyRandomGraphGenerator('E_R_pre_post','N_pre',N(1),'N_post',N(2),'p',1));
K = ones(size(I))*k;
D = ones(size(I))*1;
writeChemicalConnection(FID_syn, 1,  1,2,   I,J,K,D); % (FID, type, i_pre, j_post, I, J, K, D)

% Explanatory (ExplVar) and response variables (RespVar) for cross-simulation data gathering and post-processing
% Record explanatory variables, also called "controlled variables"
writeExplVar(FID, 'k', k);
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

