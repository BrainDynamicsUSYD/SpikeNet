function main_cc_ER_random(varargin)
% Do it!!!
% Find it!!!
% Hunt it down!!!



dt = 0.1;
sec = round(10^3/dt); % 1*(10^3/dt) = 1 sec
step_tot = 10*sec; % use 10 second!

% Loop number for PBS array job
loop_num = 0;


discard_transient = 100; % ms
tau_ref = 5;
delay = 2;

for EE_factor = [0.6 0.8 1.0]; % 0.6?
    for II_factor = [0.6 0.8 1.0]
        for EI_factor = [0.8 1.0 1.2]
            
            for  P0_init = [0.1] % 0.25 gives P0_actual = 0.2
                for I_ext_strength =  [0.8 1.0 1.2 1.4]
                    
                    
                    loop_num = loop_num + 1;
                    
                    % For PBS array job
                    if nargin ~= 0
                        PBS_ARRAYID = varargin{1};
                        if loop_num ~=  PBS_ARRAYID
                            continue;
                        end
                    end
                    
                    P_mat = [P0_init 0.3;
                        0.3 0.3];
                    
                    
                    
                    N_e = 4000;
                    sample_neuron = 1:500:4000;
                    
                    
                    N_i = 1000;
                    N = [N_e, N_i];
                    Num_pop = length(N);
                    
                    
                    
                    % seed the matlab rand function! The seed is global.
                    [FID, FID_syn] = new_ygin_files_and_randseed(loop_num);
                    
                    % write basic parameters
                    writeBasicPara(FID, dt, step_tot, N);
                    
                    
                    % write pop para
                    for pop_ind = 1:Num_pop
                        writePopPara(FID, pop_ind,  'tau_ref', tau_ref);
                        writeExtCurrentSettings(FID, pop_ind, I_ext_strength, 0);
                    end
                    
                    % write synapse para
                    writeSynPara(FID, 'tau_decay_GABA', 3);
                    
                    %%%%%%% write runaway killer
                    min_ms = 5*1000; % 5 sec
                    runaway_Hz = 20; % ??
                    Hz_ms = 1000; % ms
                    writeRunawayKiller(FID, 1, min_ms, runaway_Hz, Hz_ms);
                    %%%%%%%%%%%%%%%%%%%%%%%
                    
                    %%%%%%% data sampling
                    sample_pop = 1;
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
                    writeNeuronSampling(FID, sample_pop, [1,1,1,1,0,0,1], sample_neuron, ones(1, step_tot) )
                    
                    
                    %%%%%%% random initial condition settings (int pop_ind, double p_fire)
                    p_fire = 0.00*ones(size(N)); % between [0,1], 0.05
                    writeInitV(FID, p_fire);
                    
                    %%%%%%%%%%%%%%%%%%% Chemical Connections %%%%%%%%%%%%%%%%%%%%%%%
                    % type(1:AMAP, 2:GABAa, 3:NMDA)
                    
                    
                    
                    K_mat = [2.4*EE_factor  1.4;
                        4.5*EI_factor  5.7*II_factor]*10^-3; % miuSiemens
                    
                    Type_mat = ones(Num_pop);
                    Type_mat(end, :) = 2;
                    
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
                            
                            D = rand(size(I))*delay;
                            
                            writeChemicalConnection(FID_syn, Type_mat(i_pre, j_post),  i_pre, j_post,   I,J,K,D);
                            clear I J K D;
                        end
                    end
                    
                    
                    % Explanatory (ExplVar) and response variables (RespVar) for cross-simulation data gathering and post-processing
                    % Record explanatory variables, also called "controlled variables"
                    
                    writeExplVar(FID, 'discard_transient', discard_transient, ...
                        'loop_num', loop_num, ...
                        'EE_factor', EE_factor, ...
                        'II_factor', II_factor, ...
                        'I_ext_strength', I_ext_strength,...
                        'P0_init', P0_init, ...
                        'delay', delay, ...
                        'EI_factor', EI_factor);
                    
                    
                    % Adding comments in raster plot
                    comment1 = 'p=[0.2 0.5 0.5 0.5], k = [2.4*EE_fatcor 1.4;4.5 5.7*II_factor]*k*10^-3, tau_decay_GABA=3';
                    comment2 = datestr(now,'dd-mmm-yyyy-HH:MM');
                    writeExplVar(FID, 'comment1', comment1, 'comment2', comment2);
                    
                    
                    % append this file self into .ygin for future reference
                    appendThisMatlabFile(FID)
                    
                    disp('Matlab pre-processing done.')
                end
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

