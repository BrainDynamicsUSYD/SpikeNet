function main_phase_diagram(varargin)
% <<Asynchronous states in cortex>>

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
step_tot = 100*sec;

% Loop number for PBS array job
Num_pop = length(N);
loop_num = 0;
discard_transient = 0.05;

% Hierarchical structure
hierarchy_model = 1;
% lognormal_syn = 0;
% Inhib_range_change = 2.5;
%for Excit_strength_change = 1;
    %for Inhib_strength_change = 1; %1/Inhib_range_change * 1;
    
    for lesion_left = 0.5 %1.1:0.1:1.4 % range [0-1]
    % for P_random_control = 0:0.01:0.2
    for E_factor =  0.8:0.1:1.2;
        for I_factor = 0.8:0.1:1.2;
                for kk = 1 %2:5; % use 2 to roughly compensate synaptic saturation
                    %for pp = 0.2
                        for rr = [0.6]
                            for Mnum = 8 %!!!
                                % for I_ext_strength = 2:0.5:2.5 %0.5:0.5:2.5; %nA   run-away at 3.0!!!!
                                for rate_ext = 4.4*ones(1,5) %linspace(4.0,4.0,45) %4.0:0.025:4.5 %4.1:0.025:4.5; % Hz
                                    
                                    loop_num = loop_num + 1;
                                    
                                    % For PBS array job
                                    if nargin ~= 0
                                        PBS_ARRAYID = varargin{1};
                                        if loop_num ~=  PBS_ARRAYID
                                            continue;
                                        end
                                    end
                                    
                                    % seed the matlab rand function! The seed is global.
                                    % Be very careful about that!!!!!!!!!!!
                                    rand_seed = loop_num*10^5+eval(datestr(now,'SSFFF'));
                                    rng(rand_seed,'twister');

                                    % Creat ygin file
                                    % using loop_num in filename to ensure unique naming!
                                     % Otherwise overwriting may occur when using PBS.
                                    name = [ sprintf('%03g-', loop_num), datestr(now,'yyyymmddHHMM-SSFFF')];
                                    
                                    fprintf('Data file name is: /n%s/n', strcat(name,'.ygin') ); % write the file name to stdout and use "grep ygin" to extract it
                                    FID = fopen([name,'.ygin'], 'w'); % creat file
                                    FID_syn = fopen([name,'.ygin_syn'], 'w'); % creat file
                                    
                                    % write basic parameters
                                    writeBasicPara(FID, dt, step_tot, N)
                                    % write pop para
                                    writePopPara(FID, 1,  'tau_ref', 2);
                                    writePopPara(FID, 2,  'tau_ref', 2);
                                    % write synapse para
                                    writeSynPara(FID, 'tau_decay_GABA', 3);
                                     
%                                     %%%%%%% write runaway killer
%                                     runaway_steps = round(50/dt);
%                                     runaway_mean_num_ref = 0.4;
%                                     writeRunawayKiller(FID, runaway_steps, runaway_mean_num_ref);
                                    %%%%%%%%%%%%%%%%%%%%%%%%
                                    
                                    
                                    
                                    Kmat = [2.4*0.6  1.4;
                                            4.5  5.7*0.8]*kk*10^-3; % miuSiemens
                                        
                                    Kmat(1,:) = Kmat(1,:)*E_factor;
                                    Kmat(2,:) = Kmat(2,:)*I_factor;
                                    
                                    Pmat = [0.2 0.5;
                                            0.5 0.5];
                                    
                                    TYPEmat = [1 2];
                                    
                                    
                                    
                                    %%%%%%% external spikes settings (FID, pop_ind, type_ext, K_ext, Num_ext, rate_ext)
                                    Num_ext = 400;
                                    K_ext = Kmat(1,1);
                                    rate_t = zeros(1,step_tot);
                                    
                                    rate_t(:) = rate_ext;
                                    
%                                     rate_t(5*10^4:6*10^4) = 0;
%                                     rate_t(15*10^4:16*10^4) = 0;
                                    
                                    writeExtSpikeSettings(FID, 1, 1, K_ext,  Num_ext, rate_t);
                                    writeExtSpikeSettings(FID, 2, 1, K_ext,  Num_ext, rate_t);
          
                                    
                                    %%%%%%% random initial condition settings (int pop_ind, double p_fire)
                                    p_fire = 0.00*ones(size(N)); % between [0,1], 0.05
                                    writeInitV(FID, p_fire);
                                    
                                    %%%%%%%%%%%%%%%%%%% Chemical Connections %%%%%%%%%%%%%%%%%%%%%%%
                                    % type(1:AMAP, 2:GABAa, 3:NMDA)
                                    
                                    %
                                    for i_pre = 1:2
                                        for j_post = 1:2
                                            if hierarchy_model ~= 0
                                                if i_pre == 1 && j_post == 1 
                                                    % hierarchical structure
                                                    % A11 = MyHierarchyGraph('N', N(i_pre), 'Mnum', Mnum, 'P0', Pmat(i_pre,j_post), 'r', rr);
                                                    
                                                    %%%%%%%% Hierarchical Connection and Lesion%%%%%%%%%%%%%%%%%
                                                    % Get the guts out of
                                                    % it so I can do lesion
                                                    % study
                                                    P0 = Pmat(i_pre,j_post);
                                                
                                                    % Generate vector of first-level module size
                                                    Msize = Mnum_2_Msize(Mnum, N(i_pre));

                                                    % Generate inter modular connection probability matrix
                                                    [P, CL] = inter_module_Pmatrix(Msize, P0, rr);
                                                    
                                                    % Do lesion here
                                                    P(CL==4) = P(CL==4)*lesion_left; % highest level connection
                                                    %P(CL==3) = P(CL==3)*lesion_left;  %*lesion_left;  % second-highest level connection
                                                    %P(CL==2) = P(CL==2)*lesion_left;% third-highest level connection
                                                      %P(CL>1) = P_random_control;
                                                    

                                                    disp('Lesion in the hierarchical network!');
                                                    
                                                    % Generate full connection matrix from P
                                                    A11 = P_2_A(P,Msize);
                                                    
                                                    % Display
                                                    fprintf('Hierarchical Graph: N=%d, Mnum=%d, P0=%g, r=%g\n', N, Mnum, P0, rr);
                                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                    [I, J, ~] = find(A11);
                                                    % % save A11
                                                    % save([name,'A11.mat'], 'A11');
                                                else
                                                    [I, J, ~] = find(MyRandomGraphGenerator('E_R_pre_post','N_pre',N(i_pre),'N_post',N(j_post),'p',Pmat(i_pre,j_post)));
                                                end
                                            else
                                                [I, J, ~] = find(MyRandomGraphGenerator('E_R_pre_post','N_pre',N(i_pre),'N_post',N(j_post),'p',Pmat(i_pre,j_post)));
                                            end
                                            
%                                             if lognormal_syn == 1
%                                                 % lognorm synaptic weight
%                                                 miu = Kmat(i_pre,j_post); % mean
%                                                 sigma = miu*0.5; % std
%                                                 K = lognrnd(log(miu^2 / sqrt(sigma^2+miu^2)), sqrt(log(sigma^2/miu^2 + 1)), size(I));
%                                             else
                                                % Identical synaptic weight
                                                K = ones(size(I))*Kmat(i_pre,j_post);
%                                             end
                                            
                                            
                                            D = rand(size(I))*1;
                                            writeChemicalConnection(FID_syn, TYPEmat(i_pre),  i_pre,j_post,   I,J,K,D); % (FID, type, i_pre, j_post, I, J, K, D)
                                            clear I J K D;
                                        end
                                    end
                                    
                                    
                                    
                                    % Explanatory (ExplVar) and response variables (RespVar) for cross-simulation data gathering and post-processing
                                    % Record explanatory variables, also called "controlled variables"
                                    
                                    writeExplVar(FID, 'discard_transient', discard_transient, ...
                                                        'loop_num', loop_num, ...
                                                           'k', kk,... 
                                                           'rate_ext', rate_ext,... 
                                                           'r', rr, ...
                                                           'Mnum', Mnum, ...
                                                           'E_factor', E_factor, ...
                                                           'I_factor', I_factor,...
                                                           ...% 'P_random_control', P_random_control);
                                                            'lesion_left',lesion_left);

                                    
                                    % Adding comments in raster plot
                                    comment1 = 'p=[0.2 0.5 0.5 0.5], k = [2.4*EE_fatcor 1.4;4.5 5.7*II_factor]*k*10^-3, tau_decay_GABA=3';
                                    comment2 = datestr(now,'dd-mmm-yyyy-HH:MM');
                                    writeExplVar(FID, 'comment1', comment1, 'comment2', comment2);
                                    
                                    
                                    % append this file self into .ygin for future reference
                                    appendThisMatlabFile(FID)
                                    
                                end
                            end
                        end
                    end
                 end
             end
         end  
%     end
% end

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
