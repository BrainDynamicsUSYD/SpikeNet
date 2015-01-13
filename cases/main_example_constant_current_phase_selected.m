function main_example_constant_current_phase_selected(varargin)
% Do it!!!
% Find it!!!
% Hunt it down!!!


% varargin is for PBS arrary job
if nargin == 0
    clc;clear all;close all;
    cd /import/yossarian1/yifan/Project1/
    addpath(genpath(cd));
    cd tmp_data
end % Basic parameters
N = [4000; 1000]; %!!!
dt = 0.1;
sec = round(10^3/dt); % 1*(10^3/dt) = 1 sec
step_tot = 100*sec; % use 10 second!

% Loop number for PBS array job
Num_pop = length(N);
discard_transient = 500; % ms
EE_factor = 0.6;
II_factor = 0.8;
rr = 0.7; % this is different from 0.6!!
Mnum = 8; %!!!
kk = 1; %2:5; % use 2 to roughly compensate synaptic saturation

lesion_1 = 1.1;  %+(-0.05:0.05:0.05) %1.1:0.1:1.4 % range [0-1]
lesion_2 = 1;    %+(-0.05:0.05:0.05)
lesion_3 = 1;    %+(-0.05:0.05:0.05)
lesion_4 = 0.6;  %+(-0.05:0.05:0.05)


loop_num = 0;

for syn_i = 1:8
    syn_filename = sprintf('cc_phase_all_visiting/all_visiting-00%d.ygin_syn', syn_i);
    
    for phi_I = 0.5:0.1:1.5
        for phi_E = 0.5:0.1:1.5
            
            for I_ext_strength = 1.5*ones(1, 1)
                
                
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
                writeBasicPara(FID, dt, step_tot, N)
                % write pop para
                writePopPara(FID, 1,  'tau_ref', 2);
                writePopPara(FID, 2,  'tau_ref', 2);
                % write synapse para
                writeSynPara(FID, 'tau_decay_GABA', 3);
                
                
                %%%%%%% write runaway killer
                min_ms = 10*1000; % 10 sec
                runaway_Hz = 20; % ??
                Hz_ms = 1000; % ms
                writeRunawayKiller(FID, 1, min_ms, runaway_Hz, Hz_ms);
                % writeRunawayKiller(FID, 2, min_ms, runaway_Hz*2, Hz_ms);
                %%%%%%%%%%%%%%%%%%%%%%%
                
                
                
                Kmat = [2.4*EE_factor  1.4;
                    4.5  5.7*II_factor]*kk*10^-3; % miuSiemens
                
                Kmat(1,:) = Kmat(1,:)*phi_E;
                Kmat(2,:) = Kmat(2,:)*phi_I;
                
                Pmat = [0.2 0.5;
                    0.5 0.5];
                
                TYPEmat = [1 2];
                
                
                
                % External current
                writeExtCurrentSettings(FID, 1, I_ext_strength, 0)
                writeExtCurrentSettings(FID, 2, I_ext_strength, 0)
                
                
                
                %%%%%%% data sampling
                writeNeuronSampling(FID, 1, [1,1,1,1,0,0,1],[100:500:4000]);
                writeNeuronSampling(FID, 2, [1,1,1,1,0,0,1],[100;600]);
                
                
                %%%%%%% random initial condition settings (int pop_ind, double p_fire)
                p_fire = 0.00*ones(size(N)); % between [0,1], 0.05
                writeInitV(FID, p_fire);
                
                %%%%%%%%%%%%%%%%%%% Chemical Connections %%%%%%%%%%%%%%%%%%%%%%%
                % type(1:AMAP, 2:GABAa, 3:NMDA)
                
                %             %
                %             for i_pre = 1:2
                %                 for j_post = 1:2
                %
                %                     if i_pre == 1 && j_post == 1
                %
                %                         %%%%%%%% Hierarchical Connection and Lesion%%%%%%%%%%%%%%%%%
                %                         % Get the guts out of
                %                         % it so I can do lesion
                %                         % study
                %                         P0 = Pmat(i_pre,j_post);
                %
                %                         % Generate vector of first-level module size
                %                         Msize = Mnum_2_Msize( Mnum, N(i_pre) );
                %
                %                         % Generate inter modular connection probability matrix
                %                         [P, CL] = inter_module_Pmatrix(Msize, P0, rr);
                %
                %                         % Do lesion here
                %                         P(CL==4) = P(CL==4)*lesion_4; % highest level connection
                %                         P(CL==3) = P(CL==3)*lesion_3;  %*lesion_left;  % second-highest level connection
                %                         P(CL==2) = P(CL==2)*lesion_2; % third-highest level connection
                %                         P(CL==1) = P(CL==1)*lesion_1;
                %
                %
                %                         disp('Lesion in the hierarchical network!');
                %
                %                         % Generate full connection matrix from P
                %                         A11 = P_2_A(P,Msize);
                %
                %                         % Display
                %                         fprintf('Hierarchical Graph: N=%d, Mnum=%d, P0=%g, r=%g\n', N, Mnum, P0, rr);
                %                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                         [I, J, ~] = find(A11);
                %                         % % save A11
                %                         % save([name,'A11.mat'], 'A11');
                %                     else
                %                         [I, J, ~] = find(MyRandomGraphGenerator('E_R_pre_post','N_pre',N(i_pre),'N_post',N(j_post),'p',Pmat(i_pre,j_post)));
                %                     end
                %
                %
                %                     K = ones(size(I))*Kmat(i_pre,j_post);
                %                     D = rand(size(I))*1;
                %                     writeChemicalConnection(FID_syn, TYPEmat(i_pre),  i_pre,j_post,   I,J,K,D); % (FID, type, i_pre, j_post, I, J, K, D)
                %                     clear I J K D;
                %                 end
                %             end
                
                
                writeSynFilename(FID, syn_filename);
                
                
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
                    'lesion_4',lesion_4,'lesion_3',lesion_3,'lesion_2',lesion_2,'lesion_1',lesion_1,...
                    'phi_E', phi_E, 'phi_I', phi_I, ...
                    'syn_i', syn_i);
                
                
                % Adding comments in raster plot
                comment1 = 'p=[0.2 0.5 0.5 0.5], k = [2.4*EE_fatcor 1.4;4.5 5.7*II_factor]*k*10^-3, tau_decay_GABA=3';
                comment2 = datestr(now,'dd-mmm-yyyy-HH:MM');
                writeExplVar(FID, 'comment1', comment1, 'comment2', comment2);
                
                
                % append this file self into .ygin for future reference
                appendThisMatlabFile(FID)
                appendThisMatlabFile(FID_syn)
                
            end
        end
        
    end
    
end


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
