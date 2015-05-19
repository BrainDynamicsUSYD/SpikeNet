function main_cc_separated_modules(varargin)
% Do it!!!
% Find it!!!
% Hunt it down!!!


% varargin is for PBS arrary job
if nargin == 0
    clc;clear all;close all;
    cd /import/yossarian1/yifan/Project1/source
    addpath(genpath(cd));
    cd ../tmp_data
%     cd ~/tmp_data;
end % Basic parameters

N_I = 1000;
N_E_module = 500;
Mnum = 8; %!!!
N = [ N_E_module*ones(1, Mnum)  N_I];
Num_pop = length(N);


dt = 0.1;
sec = round(10^3/dt); % 1*(10^3/dt) = 1 sec
step_tot = 100*sec; % use 10 second!

% Loop number for PBS array job
loop_num = 0;


discard_transient = 500; % ms
EE_factor = 0.6;
II_factor = 0.8;
rr = 0.7; % this is different from 0.6!!
kk = 1; %2:5; % use 2 to roughly compensate synaptic saturation

lesion_1 = 1.1;  %+(-0.05:0.05:0.05) %1.1:0.1:1.4 % range [0-1]
lesion_2 = 1;    %+(-0.05:0.05:0.05)
lesion_3 = 1;    %+(-0.05:0.05:0.05)
lesion_4 = 0.6;  %+(-0.05:0.05:0.05)

for I_ext_strength = 1.5 * ones(1,50) % 1.3-1.35
    
    
    loop_num = loop_num + 1;
    
    % For PBS array job
    if nargin ~= 0
        PBS_ARRAYID = varargin{1};
        if loop_num ~=  PBS_ARRAYID
            continue;
        end
    end
    
    % seed the matlab rand function! The seed is global.
    [FID, FID_syn] = new_ygin_files_and_randseed(loop_num);
    
    % write basic parameters
    writeBasicPara(FID, dt, step_tot, N);
    
    
    % write pop para
    for pop_ind = 1:Num_pop
        writePopPara(FID, pop_ind,  'tau_ref', 2);
        writeExtCurrentSettings(FID, pop_ind, I_ext_strength, 0);
    end

    % write synapse para
    writeSynPara(FID, 'tau_decay_GABA', 3);
    
    %%%%%%% write runaway killer
    min_ms = 5*1000; % 5 sec
    runaway_Hz = 15; % ??
    Hz_ms = 1000; % ms
    writeRunawayKiller(FID, 1, min_ms, runaway_Hz, Hz_ms);
    %%%%%%%%%%%%%%%%%%%%%%%

    
    
    %%%%%%% data sampling
    sample_pop = 1;
%     sample_neurons = 1:2; %40:40:400;
%     sample_steps = zeros(1,step_tot);
%     sample_steps(1:100:end) = 1;
    %writeNeuronSampling(FID, sample_pop, [1,1,1,1,0,0,1], sample_neurons, sample_steps);
    writePopStatsRecord(FID, sample_pop);
    for pop_ind_pre = 1:Num_pop
        pop_ind_post = sample_pop;
        if pop_ind_pre == Num_pop
            syn_type = 2;
        else
            syn_type = 1;
        end
        %writeSynSampling(FID, pop_ind_pre, pop_ind_post, syn_type, sample_neurons, sample_steps)
        writeSynStatsRecord(FID, pop_ind_pre, pop_ind_post, syn_type)
    end
    
    
    %%%%%%% random initial condition settings (int pop_ind, double p_fire)
    p_fire = 0.00*ones(size(N)); % between [0,1], 0.05
    writeInitV(FID, p_fire);
    
    %%%%%%%%%%%%%%%%%%% Chemical Connections %%%%%%%%%%%%%%%%%%%%%%%
    % type(1:AMAP, 2:GABAa, 3:NMDA)
    
    %       Pmat = [0.2 0.5;
    %               0.5 0.5];
    % Generate inter modular connection probability matrix
    P0 = 0.2;
    Msize = N_E_module*ones(Mnum, 1);
    [P, CL] = inter_module_Pmatrix(Msize, P0, rr);
    P(CL==1) = P(CL==1)*lesion_1;
    P(CL==2) = P(CL==2)*lesion_2;
    P(CL==3) = P(CL==3)*lesion_3;
    P(CL==4) = P(CL==4)*lesion_4;
    P_mat = 0.5*ones(Num_pop);    
    P_mat(1:Mnum, 1:Mnum) = P;
    clear P;
    
    
    Type_mat = ones(Num_pop);
    Type_mat(end, :) = 2;
    
    K_mat = ones(Num_pop);
    K_mat(1:Mnum, 1:Mnum) = 2.4*EE_factor*kk*10^-3;
    K_mat(1:Mnum, end) = 1.4*kk*10^-3;
    K_mat(end, 1:Mnum) = 4.5*kk*10^-3;
    K_mat(end,end) = 5.7*II_factor*kk*10^-3;
    
    %         K_mat = [2.4*EE_factor  1.4;
    %             4.5  5.7*II_factor]*kk*10^-3; % miuSiemens

    for i_pre = 1:Num_pop
        for j_post = 1:Num_pop
            if i_pre == j_post % no self-connection!!!!!!!!
                [I, J, ~] = find(MyRandomGraphGenerator('E_R', ...
                    'N', N(i_pre), 'p', P_mat(i_pre, j_post) ));
            else
                [I, J, ~] = find(MyRandomGraphGenerator('E_R_pre_post', ...
                    'N_pre', N(i_pre),'N_post', N(j_post), 'p', P_mat(i_pre, j_post) ));
            end
            K = ones(size(I))*K_mat(i_pre,j_post);
            D = rand(size(I))*1;
            writeChemicalConnection(FID_syn, Type_mat(i_pre, j_post),  i_pre, j_post,   I,J,K,D);
            clear I J K D;
        end
    end
    
end



% Explanatory (ExplVar) and response variables (RespVar) for cross-simulation data gathering and post-processing
% Record explanatory variables, also called "controlled variables"

writeExplVar(FID, 'discard_transient', discard_transient, ...
    'loop_num', loop_num, ...
    'k', kk,...
    'r', rr, ...
    'Mnum', Mnum, ...
    'EE_factor', EE_factor, ...
    'II_factor', II_factor, ...
    'I_ext_strength', I_ext_strength, ...
    'lesion_1',lesion_1, ...
    'lesion_2',lesion_2, ...
    'lesion_3',lesion_3, ...
    'lesion_4',lesion_4);


% Adding comments in raster plot
comment1 = 'p=[0.2 0.5 0.5 0.5], k = [2.4*EE_fatcor 1.4;4.5 5.7*II_factor]*k*10^-3, tau_decay_GABA=3';
comment2 = datestr(now,'dd-mmm-yyyy-HH:MM');
writeExplVar(FID, 'comment1', comment1, 'comment2', comment2);


% append this file self into .ygin for future reference
appendThisMatlabFile(FID)

disp('Matlab pre-processing done.')
end



% This function must be here!
function appendThisMatlabFile(FID)
breaker = ['>',repmat('#',1,80)];
fprintf(FID, '%s\n', breaker);
fprintf(FID, '%s\n', '> MATLAB script generating this file: ');
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

