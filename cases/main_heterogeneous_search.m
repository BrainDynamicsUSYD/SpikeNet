
function main_heterogeneous_search(varargin)
% Do it!!!
% Find it!!!
% Hunt it down!!!



dt = 0.1;
sec = round(10^3/dt); % 1*(10^3/dt) = 1 sec
step_tot = 0.3*sec; % use 10 second!

% Loop number for PBS array job
loop_num = 0;
tau_ref = 4;
delay = 2;

% parameter
in_out_corrcoef = 0.2;
cn_scale_wire = 4;
cn_scale_weight = 4;
discard_transient = 0; % ms
STD_on = 0;
P0_init = 0.1;
degree_CV = 0.1 ; % 0.2 works
iter_num = 5;
mu_p = -0.2;  
s_p = 0.5;

%  K_ee_mean is about 0.5, need 1000 in-coming connections.
%  this is not good.
%  what can I do??? 
%  ref: A Lognormal Recurrent Network Model for Burst Generation during Hippocampal Sharp Waves

for g_EI_over_EE = 0.25:0.05:0.5
    for g_IE = 4:1:7 % 6
        for g_II = 1.5:0.5:2.5 % 1.5
            for I_ext_CV = 0
                for I_ext_strength_E = 0 %[ 0.05:0.05:0.2 ]% 0.9*ones(1,10)]
                    for I_ext_strength_I = 0 %[0.1]% 0.9*ones(1,10)]
                        for  tau_c = [10 15]
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
                            
                            K_mat = [NaN  g_IE;
                                     NaN  g_II]*10^-3; % miuSiemens
                            P_mat = [P0_init 0.3;
                                        0.3  0.3];
                            
                            % sptially embedded network
                            hw = 31; % half-width, (31*2+1)^2 = 3969 ~ 4000, hw=49 gives 9801
                            N_e = (hw*2+1)^2;
                            
                            N_i = 1000;
                            N = [N_e, N_i];
                            Num_pop = length(N);
                            Type_mat = ones(Num_pop);
                            Type_mat(end, :) = 2;
                            
                            
                            % generate in- and out-degree accroding to (1) distance-dependent rule
                            % (2) degree distributios and (3) common neighbour rule
                            deg_mean = N_e*P0_init;
                            deg_std = degree_CV*deg_mean;
                            [ degree_in_0, degree_out_0 ] = logn_in_out_degree( N_e, [deg_mean deg_mean], [deg_std,  deg_std], in_out_corrcoef  );% log-normal
                            [ I_ee, J_ee, dist_IJ, iter_hist ] = generate_IJ_2D( degree_in_0, degree_out_0, tau_c, cn_scale_wire, iter_num );
                            in_degree = full(sum(sparse(I_ee,J_ee,ones(size(I_ee))), 1)); % row vector
                            out_degree = full(sum(sparse(I_ee,J_ee,ones(size(I_ee))), 2));
                            
                            % Generate K according to
                            
                            g_pool_generator_hld = @(N)g_pool_generator(N, mu_p, s_p);
                            K_scale = sqrt(in_degree);
                            K_cell = inverse_pool( in_degree, K_scale, g_pool_generator_hld);
                             K_ee = NaN;
                            if ~isnan(K_cell{1})
                                 K_ee = zeros(size(J_ee)); for j = 1:N_e;  K_ee(J_ee==j) = K_cell{j}'; end, clear K_cell; % reformat K
                            end
                            % shuffle K accordind to common neighbour rule
                            if ~isnan( K_ee)
                                [  K_ee ] = shuffle_K_common_neighbour(  K_ee, I_ee, J_ee, cn_scale_weight );
                            end
                            K_ee_mean = mean(K_ee);
                            
                            
                            [~,ind_sorted] = sort(in_degree);
                            sample_neuron = ind_sorted(1:500:end);

                            % write basic parameters
                            writeBasicPara(FID, dt, step_tot, N);
                            
                            
                            % write pop para
                            for pop_ind = 1:Num_pop
                                writePopPara(FID, pop_ind,  'tau_ref', tau_ref);
                            end
                            
                            % write external currents
                            I_ext_E = I_ext_strength_E*sqrt(in_degree/mean(in_degree));
                            I_ext_I = I_ext_strength_I*ones(1,N_i);
                            writeExtCurrentSettings(FID, 1, I_ext_E, I_ext_E*I_ext_CV);
                            writeExtCurrentSettings(FID, 2, I_ext_I, I_ext_I*I_ext_CV);
                            
                            % write synapse para
                            writeSynPara(FID, 'tau_decay_GABA', 3);
                            
                            %%%%%%% write runaway killer
                            min_ms = 1000; % 5 sec
                            runaway_Hz = 200; % ??
                            Hz_ms = 500; % ms
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
                            p_fire = [0.1 0.05]; % between [0,1], 0.05
                            writeInitV(FID, p_fire);
                            
                            %%%%%%%%%%%%%%%%%%% Chemical Connections %%%%%%%%%%%%%%%%%%%%%%%
                            % type(1:AMAP, 2:GABAa, 3:NMDA)
                            
                            
                            
                            for i_pre = 1:Num_pop
                                for j_post = 1:Num_pop
                                    if i_pre == 1 && j_post == 1
                                        I = I_ee;
                                        J = J_ee;
                                    elseif i_pre == j_post % no self-connection!!!!!!!!
                                        [I, J, ~] = find(MyRandomGraphGenerator('E_R', ...
                                            'N', N(i_pre), 'p', P_mat(i_pre, j_post) ));
                                    else
                                        [I, J, ~] = find(MyRandomGraphGenerator('E_R_pre_post', ...
                                            'N_pre', N(i_pre),'N_post', N(j_post), 'p', P_mat(i_pre, j_post) ));
                                    end
                                    
                                    if i_pre == 1 && j_post == 1 % E to E
                                        EE_input = full(sum(sparse(I_ee,J_ee,K_ee),1));
                                        K = K_ee;
                                        
                                    elseif i_pre == 2 && j_post == 1 % I to E
                                        K = zeros(size(J));
                                        for i_E = 1:N(1)
                                            K(J==i_E) = EE_input(i_E)/sum(J==i_E)*g_EI_over_EE;
                                        end
 
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
                                'g_EI_over_EE', g_EI_over_EE, ...
                                'g_IE', g_IE, ...
                                'g_II', g_II, ...
                                'I_ext_strength_E', I_ext_strength_E,...
                                'I_ext_strength_I', I_ext_strength_I,...
                                'P0_init', P0_init, ...
                                'degree_CV', degree_CV,...
                                'in_out_corrcoef', in_out_corrcoef, ...
                                'tau_c', tau_c, ...
                                'STD_on', STD_on, ...
                                'I_ext_CV', I_ext_CV, ...
                                'cn_scale_wire', cn_scale_wire, ...
                                'cn_scale_weight', cn_scale_weight, ...
                                'mu_p', mu_p,...
                                's_p', s_p);
                            
                            
                            % Adding comments in raster plot
                            comment1 = 'p=[p0_init 0.3 0.3 0.3], k = ?, tau_decay_GABA=3';
                            comment2 = datestr(now,'dd-mmm-yyyy-HH:MM');
                            writeExplVar(FID, 'comment1', comment1, 'comment2', comment2);
                            
                            % save in_degree and sample neuron data based on in_degree
                            save([sprintf('%04g-', loop_num), datestr(now,'yyyymmddHHMM-SSFFF'),...
                                '_in_degree.mat'], 'in_degree', 'out_degree', 'dist_IJ', 'iter_hist','K_ee_mean');
                            
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


