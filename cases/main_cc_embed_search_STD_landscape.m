function main_cc_embed_search_STD_landscape(varargin)
% Do it!!!
% Find it!!!
% Hunt it down!!!



dt = 0.1;
sec = round(10^3/dt); % 1*(10^3/dt) = 1 sec
step_tot = 40*sec; % use 10 second!

% Loop number for PBS array job
loop_num = 0;
tau_ref = 4;
delay = 2;

in_deg_scale_exp = -0.5;


xc_m = (1:2:9)*7;
yc_m = (1:2:9)*7;
[xc_m, yc_m] = meshgrid(xc_m, yc_m);
xc_m = xc_m(:);
yc_m = yc_m(:);

hw = 31; % half-width, (31*2+1)^2 = 3969 ~ 4000
[Lattice, ~] = lattice_nD(2, hw);

for phi_E = 1
    phi_I = phi_E;
    for STD_on = [1 0]
        discard_transient = 10; % ms
        for EE_factor = [0.5 ]; % 0.6?
            for II_factor = [0.7 ]
                for EI_factor = [ 0.8 ]
                    for degree_CV = [0.1 ] % 0.5?
                        for  P0_init = [ 0.08 ] % 0.25 gives P0_actual = 0.2
                            for I_ext_strength = [ 0.8 0.8]% 0.9*ones(1,10)]
                                for  tau_c = [10 ]
                                    for i_centre = 1:length(xc_m)
                                        for  rate_ext = [0.2 0.4]
                                            
                                            loop_num = loop_num + 1;
                                            
                                            % For PBS array job
                                            if nargin ~= 0
                                                PBS_ARRAYID = varargin{1};
                                                if loop_num ~=  PBS_ARRAYID
                                                    continue;
                                                end
                                            end
                                            
                                            % seed the matlab rand function! The seed is global.
                                            [FID, ~] = new_ygin_files_and_randseed(loop_num,0);
                                            
                                            writeSynFilename(FID, 'sample.ygin_syn')
                                            
                                            
                                            xc = xc_m(i_centre);
                                            yc = yc_m(i_centre);
                                            dist = lattice_nD_find_dist(Lattice, hw, [xc, yc]);
                                            neurons_trans_input = dist < tau_c*0.6;
                                            
                                            N_ext = 1000;
                                            g_ext = 2*10^-3;
                                            rate_ext_t = zeros(1,step_tot);
                                            rate_ext_t(0.5*sec:0.6*sec) = rate_ext;
                                            writeExtSpikeSettings(FID, 1, 1, g_ext,  N_ext, rate_ext_t, neurons_trans_input);
                                            
                                            % sptially embedded network
                                            
                                            
                                            
                                            [ll_matrix, N_e, ~] = spatial_embed_network(hw, P0_init, degree_CV, tau_c);
                                            [I_e,J_e,~] = find(ll_matrix);
                                            
                                            A = sparse(I_e,J_e, ones(size(I_e)));
                                            in_degree = full(sum(A, 1));
                                            %                                     in_deg_ratio = in_degree/mean(in_degree);
                                            %                                     in_deg_ratio = in_deg_ratio(J_e);
                                            %
                                            %                                     in_deg_scaling = in_deg_ratio(:).^in_deg_scale_exp; % use -0.5 maybe?
                                            
                                            %                                     % save in_degree and sample neuron data based on in_degree
                                            %                                     save([sprintf('%04g-', loop_num), datestr(now,'yyyymmddHHMM-SSFFF'),...
                                            %                                         '_in_degree.mat'], 'in_degree');
                                            [~,ind_sorted] = sort(in_degree);
                                            sample_neuron = ind_sorted(1:1000:end);
                                            
                                            
                                            N_i = 1000;
                                            N = [N_e, N_i];
                                            Num_pop = length(N);
                                            
                                            
                                            
                                            % write basic parameters
                                            writeBasicPara(FID, dt, step_tot, N);
                                            
                                            
                                            % write pop para
                                            for pop_ind = 1:Num_pop
                                                writePopPara(FID, pop_ind,  'tau_ref', tau_ref);
                                                writeExtCurrentSettings(FID, pop_ind, I_ext_strength*ones(1,N(pop_ind)), 0*ones(1,N(pop_ind)));
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
                                            writeNeuronSampling(FID, sample_pop, [1,1,1,1,0,0,1, 1], sample_neuron, ones(1, step_tot) )
                                            
                                            
                                            %%%%%%% random initial condition settings (int pop_ind, double p_fire)
                                            p_fire = -2.00*ones(size(N)); % between [0,1], 0.05
                                            writeInitV(FID, p_fire);
                                            
                                            %%%%%%%%%%%%%%%%%%% Chemical Connections %%%%%%%%%%%%%%%%%%%%%%%
                                            % type(1:AMAP, 2:GABAa, 3:NMDA)
                                            
                                            %                                     P_mat = [0.1 0.3;
                                            %                                         0.3 0.3];
                                            %
                                            %                                     K_mat = [2.4*EE_factor  1.4;
                                            %                                         4.5*EI_factor  5.7*II_factor]*10^-3; % miuSiemens
                                            %
                                            %                                     K_mat(1,:) = K_mat(1,:)*phi_E;
                                            %                                     K_mat(2,:) = K_mat(2,:)*phi_I;
                                            %
                                            %                                     Type_mat = ones(Num_pop);
                                            %                                     Type_mat(end, :) = 2;
                                            %
                                            %                                     for i_pre = 1:Num_pop
                                            %                                         for j_post = 1:Num_pop
                                            %                                             if i_pre == 1 && j_post == 1
                                            %                                                 I = I_e;
                                            %                                                 J = J_e;
                                            %                                             elseif i_pre == j_post % no self-connection!!!!!!!!
                                            %                                                 [I, J, ~] = find(MyRandomGraphGenerator('E_R', ...
                                            %                                                     'N', N(i_pre), 'p', P_mat(i_pre, j_post) ));
                                            %                                             else
                                            %                                                 [I, J, ~] = find(MyRandomGraphGenerator('E_R_pre_post', ...
                                            %                                                     'N_pre', N(i_pre),'N_post', N(j_post), 'p', P_mat(i_pre, j_post) ));
                                            %                                             end
                                            %
                                            %                                             if i_pre == 1 && j_post == 1
                                            %                                                 K = ones(size(I))*K_mat(i_pre,j_post).*in_deg_scaling;
                                            %                                             else
                                            %                                                 K = ones(size(I))*K_mat(i_pre,j_post);
                                            %                                             end
                                            %                                             D = rand(size(I))*delay;
                                            %                                             writeChemicalConnection(FID_syn, Type_mat(i_pre, j_post),  i_pre, j_post,   I,J,K,D);
                                            %                                             clear I J K D;
                                            %                                         end
                                            %                                     end
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
                                                'I_ext_strength', I_ext_strength,...
                                                'P0_init', P0_init, ...
                                                'degree_CV', degree_CV,...
                                                'in_deg_scale_exp', in_deg_scale_exp, ...
                                                'tau_c', tau_c, ...
                                                'STD_on', STD_on, ...
                                                'xc', xc, 'yc', yc, ...
                                                'rate_ext', rate_ext);
                                            
                                            
                                            % Adding comments in raster plot
                                            comment1 = 'tau_decay_GABA=3';
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

