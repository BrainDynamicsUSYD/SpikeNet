function main_cc_in_out_search(varargin)
% Do it!!!
% Find it!!!
% Hunt it down!!!



dt = 0.1;
sec = round(10^3/dt); % 1*(10^3/dt) = 1 sec
step_tot = 0.5*sec; % use 10 second!

% Loop number for PBS array job
loop_num = 0;
tau_ref = 4;
delay = 2;


in_out_corrcoef = 0.1;
discard_transient = 0; % ms

for I_ext_CV = 0
    for EE_CV = 0.2
        for EI_CV = 0.2
            for phi_E = 0.4
                phi_I = phi_E;
                for STD_on = [0]
                    for EE_factor = [0.6 ]; % 0.6?
                        for II_factor = [ 0.6 0.8]
                            for EI_factor = [ 1.0 ]
                                for degree_CV = [0.5 ] % 0.5?
                                    for  P0_init = [ 0.1 ]
                                        for I_ext_strength_E = [ 0.2:0.1:0.4 ]% 0.9*ones(1,10)]
                                            for I_ext_strength_I = [0.1:0.1:0.2]% 0.9*ones(1,10)]
                                                for  tau_c = [10 15]
                                                    loop_num = loop_num + 1;
                                                    
                                                    % For PBS array job
                                                    if nargin ~= 0
                                                        PBS_ARRAYID = varargin{1};
                                                        if loop_num ~=  PBS_ARRAYID
                                                            continue;
                                                        end
                                                    end
                                                    
                                                    
                                                    % sptially embedded network
                                                    
                                                    hw = 31; % half-width, (31*2+1)^2 = 3969 ~ 4000
                                                    
                                                    [ll_matrix, N_e, ~] = spatial_embed_in_out_network(hw, P0_init, degree_CV, tau_c, in_out_corrcoef);
                                                    [I_e,J_e,~] = find(ll_matrix);
                                                    
                                                    A = sparse(I_e,J_e, ones(size(I_e)));
                                                    in_degree = full(sum(A, 1));
                                                    out_degree = full(sum(A, 2));
                                                    
                                                    
                                                    
                                                    [~,ind_sorted] = sort(in_degree);
                                                    sample_neuron = ind_sorted(1:500:end);
                                                    
                                                    
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
                                                    end
                                                    
                                                    % write external currents
                                                    I_ext = I_ext_strength_E*sqrt(in_degree/mean(in_degree));
                                                    writeExtCurrentSettings(FID, 1, I_ext, I_ext*I_ext_CV);
                                                    writeExtCurrentSettings(FID, 2, I_ext_strength_I, I_ext_strength_I*I_ext_CV);
                                                    
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
                                                    
                                                    writeNeuronSampling(FID, 2, [1,1,1,1,0,0,1], [1 100], ones(1, step_tot) )
                                                    
                                                    %%%%%%% random initial condition settings (int pop_ind, double p_fire)
                                                    p_fire = -1.00*ones(size(N)); % between [0,1], 0.05
                                                    writeInitV(FID, p_fire);
                                                    
                                                    %%%%%%%%%%%%%%%%%%% Chemical Connections %%%%%%%%%%%%%%%%%%%%%%%
                                                    % type(1:AMAP, 2:GABAa, 3:NMDA)
                                                    
                                                    P_mat = [P0_init 0.3;
                                                        0.3 0.3];
                                                    
                                                    K_mat = [2.4*EE_factor  1.4;
                                                        4.5  5.7*II_factor]*10^-3; % miuSiemens
                                                    
                                                    Type_mat = ones(Num_pop);
                                                    Type_mat(end, :) = 2;
                                                    
                                                    for i_pre = 1:Num_pop
                                                        for j_post = 1:Num_pop
                                                            if i_pre == 1 && j_post == 1
                                                                I = I_e;
                                                                J = J_e;
                                                            elseif i_pre == j_post % no self-connection!!!!!!!!
                                                                [I, J, ~] = find(MyRandomGraphGenerator('E_R', ...
                                                                    'N', N(i_pre), 'p', P_mat(i_pre, j_post) ));
                                                            else
                                                                [I, J, ~] = find(MyRandomGraphGenerator('E_R_pre_post', ...
                                                                    'N_pre', N(i_pre),'N_post', N(j_post), 'p', P_mat(i_pre, j_post) ));
                                                            end
                                                            
                                                            if i_pre == 1 && j_post == 1 % E to E
                                                                K_tot =  K_mat(1,1)*length(I_e);
                                                                scale = sqrt(in_degree);
                                                                K = scale*K_tot/sum(scale);
                                                                K_mu = K./in_degree;
                                                                K_ee = zeros(size(J_e));
                                                                for j = 1:length(in_degree)
                                                                    Mu_log = K_mu(j); %/(mean(in_degree));
                                                                    Sigma_log = Mu_log*EE_CV;
                                                                    Mu_norm = log((Mu_log.^2)./sqrt(Sigma_log.^2+Mu_log.^2));
                                                                    Sigma_norm = sqrt(log(Sigma_log.^2./(Mu_log.^2)+1));
                                                                    err = 1;
                                                                    while err > 1e-2
                                                                        K_tmp = exp(randn(size(1,in_degree(j)))*Sigma_norm + Mu_norm);
                                                                        err = abs(K_tmp - Mu_log)/Mu_log;
                                                                    end
                                                                    K_ee(J_e == j) = K_tmp;
                                                                    
                                                                end
                                                                EE_input = full(sum(sparse(I_e,J_e,K_ee),1));
                                                                K = K_ee;
                                                                
                                                            elseif i_pre == 2 && j_post == 1 % I to E
                                                                K_ei = [];
                                                                for E_i = 1:N(1)
                                                                    Mu_log = EE_input(E_i)*EI_factor/(P_mat(2,1)*N(2));
                                                                    Sigma_log = Mu_log*EI_CV; % EI_CV?
                                                                    Mu_norm = log((Mu_log.^2)./sqrt(Sigma_log.^2+Mu_log.^2));
                                                                    Sigma_norm = sqrt(log(Sigma_log.^2./(Mu_log.^2)+1));
                                                                    K_ei = [K_ei exp(randn(1,sum(J==E_i))*Sigma_norm + Mu_norm)];
                                                                    K = K_ei;
                                                                end
                                                                EI_input = full(sum(sparse(I,J,K_ei),1));
                                                                
                                                            else
                                                                K = ones(size(I))*K_mat(i_pre,j_post);
                                                            end
                                                            D = rand(size(I))*delay;
                                                            writeChemicalConnection(FID_syn, Type_mat(i_pre, j_post),  i_pre, j_post,   I,J,K,D);
                                                            clear I J K D;
                                                        end
                                                    end
                                                    if STD_on == 1
                                                        % writeSTD(FID, 1, 1, 5*sec);
                                                        writeSTD(FID, 1, 1, 1);
                                                    end
                                                    % Explanatory (ExplVar) and response variables (RespVar) for cross-simulation data gathering and post-processing
                                                    % Record explanatory variables, also called "controlled variables"
                                                    
                                                    writeExplVar(FID, 'discard_transient', discard_transient, ...
                                                        'loop_num', loop_num, ...
                                                        'delay', delay, ...
                                                        'phi_I', phi_I, ...
                                                        'phi_E', phi_E, ...
                                                        'EE_factor', EE_factor, ...
                                                        'EI_factor', EI_factor, ...
                                                        'II_factor', II_factor, ...
                                                        'I_ext_strength_E', I_ext_strength_E,...
                                                        'I_ext_strength_I', I_ext_strength_I,...
                                                        'P0_init', P0_init, ...
                                                        'degree_CV', degree_CV,...
                                                        'EE_CV',EE_CV,...
                                                        'EI_CV',EI_CV,...
                                                        'in_out_corrcoef', in_out_corrcoef, ...
                                                        'tau_c', tau_c, ...
                                                        'STD_on', STD_on, ...
                                                        'I_ext_CV', I_ext_CV);
                                                    
                                                    
                                                    % Adding comments in raster plot
                                                    comment1 = 'p=[0.2 0.5 0.5 0.5], k = [2.4*EE_fatcor 1.4;4.5 5.7*II_factor]*k*10^-3, tau_decay_GABA=3';
                                                    comment2 = datestr(now,'dd-mmm-yyyy-HH:MM');
                                                    writeExplVar(FID, 'comment1', comment1, 'comment2', comment2);
                                                    
                                                    
                                                    % save in_degree and sample neuron data based on in_degree
                                                    save([sprintf('%04g-', loop_num), datestr(now,'yyyymmddHHMM-SSFFF'),...
                                                        '_in_degree.mat'], 'in_degree', 'out_degree', 'EE_input','EI_input');
                                                    
                                                    % append this file self into .ygin for future reference
                                                    appendThisMatlabFile(FID)
                                                    
                                                    disp('Matlab pre-processing done.')
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
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


