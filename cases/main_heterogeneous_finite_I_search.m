function main_heterogeneous_finite_I_search(varargin)
% Do it!!!
% Find it!!!
% Hunt it down!!!



dt = 0.1;
sec = round(10^3/dt); % 1*(10^3/dt) = 1 sec

step_tot = 200*sec; % use 10 second!
discard_transient = 0; % ms

% Loop number for PBS array job
loop_num = 0;
tau_ref = 4;
delay = 4;


P0_init = 0.05;

P_mat = [P0_init 0.1;
            0.1  0.2];

% sptially embedded network
hw = 44; % half-width, (31*2+1)^2 = 3969 ~ 4000, hw=44 gives 7921
N_e = (hw*2+1)^2; %

N_i = 2000;
N = [N_e, N_i];
Num_pop = length(N);
Type_mat = ones(Num_pop);
Type_mat(end, :) = 2;


% parameter
for SpikeFreqAapt = [0 1]
    for in_out_r = [0.2 ];
        for cn_scale_wire = [2 ];
            for cn_scale_weight = [2 ];
                iter_num = 5;
                
                
                STD_on = 0;
                
                
                N_ext = 1000;
                g_ext = 2*10^-3;
                
                [ fit_g_2_EPSP_2, ~ ] = g_EPSP_conversion( );
                
                
                for deg_hybrid = [0.4 ]
                    degree_CV = 0.2; % 0.2 works
                    
                    for g_mu = [4]*10^-3;
                        
                        
                        EPSP_mu = fit_g_2_EPSP_2(g_mu);
                        EPSP_sigma = 1;
                        
                        
                        
                        for inh_STDP = [0 ];
                            
                            
                            %  K_ee_mean is about 0.5, need 1000 in-coming connections.
                            %  this is not good.
                            %  what can I do???
                            %  ref: A Lognormal Recurrent Network Model for Burst Generation during Hippocampal Sharp Waves
                            
                            
                            for g_EI = [ 11 12 ]*10^-3
                                for g_IE = [5]*10^-3
                                    for g_II = [25]*10^-3
                                        
                                        for rate_ext = [1.6:0.1:2.5];
                                            for  tau_c_E = [8 10]
                                                for tau_c_I = [15 20]
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
                                                        g_EI  g_II]; % miuSiemens
                                                    
                                                    if SpikeFreqAapt == 1
                                                        writeSpikeFreqAdpt(FID, 1);
                                                    end
                                                    
                                                    % write basic parameters
                                                    writeBasicPara(FID, dt, step_tot, N);
                                                    
                                                    
                                                    % write pop para
                                                    for pop_ind = 1:Num_pop
                                                        writePopPara(FID, pop_ind,  'tau_ref', tau_ref);
                                                    end
                                                    
                                                    % write external currents
                                                    writeExtSpikeSettings(FID, 1, 1, g_ext,  N_ext, rate_ext*ones(1, step_tot),  1, N(1) );
                                                    writeExtSpikeSettings(FID, 2, 1, g_ext,  N_ext, rate_ext*ones(1, step_tot),  1, N(2) );
                                                    
                                                    % write synapse para
                                                    writeSynPara(FID, 'tau_decay_GABA', 3);
                                                    
                                                    % inhibitory STDP
                                                    if inh_STDP == 1
                                                        writeInhSTDP(FID, 2, 1, 1*sec);
                                                    end
                                                    
                                                    if STD_on == 1
                                                        % writeSTD(FID, 1, 1, 5*sec);
                                                        writeSTD(FID, 1, 1, 1);
                                                    end
                                                    
                                                    % write runaway killer
                                                    min_ms = 500; % 5 sec
                                                    runaway_Hz = 40; % ??
                                                    Hz_ms = 200; % ms
                                                    writeRunawayKiller(FID, 1, min_ms, runaway_Hz, Hz_ms);
                                                    
                                                    % random initial condition settings (int pop_ind, double p_fire)
                                                    p_fire = [0.1 0.00]; % between [0,1], 0.05
                                                    writeInitV(FID, p_fire);
                                                    
                                                    
                                                    
                                                    %%%%%%%%%%%%%%%%%%% Chemical Connections %%%%%%%%%%%%%%%%%%%%%%%
                                                    % type(1:AMAP, 2:GABAa, 3:NMDA)
                                                    
                                                    Lattice_I = quasi_lattice_2D( N(2) , hw);
                                                    %%%%%%%%%%%%%%%%%%%%%%
                                                    % generate in- and out-degree accroding to (1) distance-dependent rule
                                                    % (2) degree distributios and (3) common neighbour rule
                                                    deg_mean = N_e*P0_init;
                                                    deg_std_logn = degree_CV*deg_mean;
                                                    
                                                    [ deg_in_0, deg_out_0 ] = hybrid_degree( N_e, deg_mean, deg_std_logn, in_out_r, deg_hybrid );
                                                    
                                                    [ I_ee, J_ee, dist_IJ, iter_hist, Lattice_E ] = generate_IJ_2D( deg_in_0, deg_out_0, tau_c_E, cn_scale_wire, iter_num );
                                                    in_degree = full(sum(sparse(I_ee,J_ee,ones(size(I_ee))), 1)); % row vector
                                                    out_degree = full(sum(sparse(I_ee,J_ee,ones(size(I_ee))), 2));
                                                    
                                                    % Generate K according to
                                                    mu_p = log((EPSP_mu^2)/sqrt(EPSP_sigma^2+EPSP_mu^2));
                                                    s_p = sqrt(log(EPSP_sigma^2/(EPSP_mu^2)+1));
                                                    mu_p = mu_p + s_p^2;
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
                                                    EE_input = full(sum(sparse(I_ee,J_ee,K_ee),1));
                                                    
                                                    D = rand(size(I_ee))*delay;
                                                    writeChemicalConnection(FID_syn, Type_mat(1, 1),  1, 1,   I_ee,J_ee,K_ee,D);
                                                    clear I J K D;
                                                    
                                                    [~,ind_sorted] = sort(in_degree);
                                                    sample_neuron = ind_sorted(1:250:end);
                                                    
                                                    %%%%%%%%%%%%%%%%%%%%%%
                                                    [ I,J ] = Lattice2Lattice( Lattice_I, Lattice_E, hw, tau_c_I, P_mat(2,1) );
                                                    D = rand(size(I))*delay;
                                                    K = zeros(size(J));
                                                    for i_E = 1:N(1)
                                                        mu_K_tmp = EE_input(i_E)/sum(J==i_E)*(g_EI/g_mu);
                                                        K(J==i_E) = abs(randn([1 sum(J==i_E)])*(mu_K_tmp/4) + mu_K_tmp); % this is a bit too arbitary!
                                                    end
                                                    % K = ones(size(I))*K_mat(2,1);
                                                    writeChemicalConnection(FID_syn, Type_mat(2, 1),  2, 1,   I,J,K,D);
                                                    clear I J K D;
                                                    
                                                    %%%%%%%%%%%%%%%%%%%%%%
                                                    [ I,J ] = Lattice2Lattice( Lattice_E, Lattice_I, hw, tau_c_E, P_mat(1,2) );
                                                    D = rand(size(I))*delay;
                                                    K = ones(size(I))*K_mat(1,2);
                                                    writeChemicalConnection(FID_syn, Type_mat(1, 2),  1, 2,   I,J,K,D);
                                                    clear I J K D;
                                                    
                                                    %%%%%%%%%%%%%%%%%%%%%%
                                                    [ I,J ] = Lattice2Lattice( Lattice_I, Lattice_I, hw, tau_c_I, P_mat(2,2) );
                                                    D = rand(size(I))*delay;
                                                    K = ones(size(I))*K_mat(2,2);
                                                    writeChemicalConnection(FID_syn, Type_mat(2, 2),  2, 2,   I,J,K,D);
                                                    clear I J K D;
                                                    
                                                    
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
                                                    writeNeuronSampling(FID, sample_pop, [1,1,1,1,0,0,1, 0], sample_neuron, ones(1, step_tot) )
                                                    writeNeuronSampling(FID, 2, [1,1,1,1,0,0,1,0], [1 100], ones(1, step_tot) )
                                                    
                                                    
                                                    % Explanatory (ExplVar) and response variables (RespVar) for cross-simulation data gathering and post-processing
                                                    % Record explanatory variables, also called "controlled variables"
                                                    
                                                    writeExplVar(FID, 'discard_transient', discard_transient, ...
                                                        'loop_num', loop_num, ...
                                                        'delay', delay, ...
                                                        'g_mu', g_mu, ...
                                                        'g_EI', g_EI, ...
                                                        'g_IE', g_IE, ...
                                                        'g_II', g_II, ...
                                                        'g_ext', g_ext,...
                                                        'SpikeFreqAapt', SpikeFreqAapt, ...
                                                        'rate_ext', rate_ext,...
                                                        'P0_init', P0_init, ...
                                                        'degree_CV', degree_CV,...
                                                        'in_out_r', in_out_r, ...
                                                        'tau_c_E', tau_c_E, ...
                                                        'tau_c_I', tau_c_I, ...
                                                        'STD_on', STD_on, ...
                                                        'cn_scale_wire', cn_scale_wire, ...
                                                        'cn_scale_weight', cn_scale_weight, ...
                                                        'mu_p', mu_p,...
                                                        's_p', s_p, ...
                                                        'inh_STDP', inh_STDP, ...
                                                        'deg_hybrid', deg_hybrid);
                                                    
                                                    
                                                    % Adding comments in raster plot
                                                    comment1 = ' '; %'p=[p0_init 0.3 0.3 0.3], k = ?, tau_decay_GABA=3';
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

