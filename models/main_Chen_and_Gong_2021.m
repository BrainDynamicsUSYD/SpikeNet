function main_Chen_Gong_2021(varargin)

dt = 0.1;
sec = round(10^3/dt); % 1*(10^3/dt) = 1 sec

step_tot = 10*sec;%7.1*sec; % use 10 second!
discard_transient = 0; % ms

% Loop number for PBS array job
loop_num = 0;
tau_ref = 4;
delay = 4;

%>>>>>>>>> stimuli's parameters >>>>>>>>>>>>

qisig = 0.6;
qi_x_offset = pi;
qi_pos = [qi_x_offset qi_x_offset; 0 0];


dg_K = 0.003;
for g_balance = 0.98
    qi_contrast = 1;
    qi_c_gradient = 0;
    qi_contrasts = [qi_contrast; qi_contrast + qi_c_gradient];
    repeats = 5;
    for P0_init = 0.08*ones(1,repeats)
        
        P_mat = [P0_init 0.1;
            0.1  0.2]*2;
        
        % sptially embedded network
        hw = 31; % half-width, (31*2+1)^2 = 3969 ~ 4000, hw=44 gives 7921
        N_e = (hw*2+1)^2; %
        
        N_i = 1000;
        N = [N_e, N_i];
        Num_pop = length(N);
        Type_mat = ones(Num_pop);
        Type_mat(end, :) = 2;
        in_out_r = [0.13 ];% correlation between in-degree and out-degree
        
        % parameter
        for SpikeFreqAapt = [ 1] % introduce potassium current
            
            for LFP_range_sigma = [8]; % 8
                for cn_scale_wire = [2 ];
                    for cn_scale_weight = [2 ];
                        iter_num = 5;
                        
                        
                        STD_on = 0;
                        
                        
                        N_ext = 1000;
                        g_ext = 2*10^-3;
                        
                        [ fit_g_2_EPSP_2, ~ ] = g_EPSP_conversion( );%????
                        
                        
                        for deg_hybrid = [0.4 ]
                            degree_CV = 0.2; % 0.2 works
                            
                            for g_mu = [4]*10^-3;                                
                                
                                EPSP_mu = fit_g_2_EPSP_2(g_mu);
                                EPSP_sigma = 1;  
                                
                                for inh_STDP = [0 ]; 
                                    for g_EI = [ 13.5 ]*10^-3 %11 12
                                        for g_IE = [5 ]*10^-3
                                            for g_II = [25]*10^-3
                                                
                                                for rate_ext_I = [1];
                                                    for rate_ext_E = [0.85 ]; %
                                                        for  tau_c_EE = [8]
                                                            tau_c_IE = 10;
                                                            for tau_c_I = [20]
                                                                loop_num = loop_num + 1;
                                                                
                                                                % For PBS array job
                                                                if nargin ~= 0
                                                                    PBS_ARRAYID = varargin{1};
                                                                    if loop_num ~=  PBS_ARRAYID
                                                                        continue;
                                                                    end
                                                                end                                                                
                                                                
                                                                [FID] = new_ygin_files_and_randseedHDF5(loop_num);
                                                                % write basic parameters
                                                                writeBasicParaHDF5(FID, dt, step_tot, N);
                                                                
                                                                K_mat = [NaN  g_IE;
                                                                    g_EI  g_II]; % miuSiemens
                                                                
                                                                if SpikeFreqAapt == 1
                                                                    writeSpikeFreqAdptHDF5(FID, 1,dg_K);
                                                                end
                                                                
                                                                % write pop para
                                                                for pop_ind = 1:Num_pop
                                                                    writePopParaHDF5(FID, pop_ind,  'tau_ref', tau_ref);
                                                                end
                                                                
                                                                % write external spikes of base line    
                                                                rate_ext_on_base1= ones(1,step_tot); % duration of stimulus
                                                                writeExtSpikeTinvSettingsHDF5(FID, 1, 1, g_ext,  N_ext, rate_ext_E*ones(1,N_e), rate_ext_on_base1);
                                                                %                                                                     writeExtCurrentSettingsHDF5(FID, 2, 1, 0)%, modify)
                                                                writeExtSpikeSettingsHDF5(FID, 2, 1, g_ext,  N_ext, rate_ext_I*ones(1, step_tot),  ones(1, N(2)) );
                                                                % add sti as external spikes
                                                                qix = (-31:31)*2*pi/63; %spatial coordinate in (-pi,pi)
                                                                rate_ext_E_scale = 1;%0.5; %global scaling factor to input (default 0.85 kHz)
                                                                rate_ext_on = zeros(1,step_tot);
                                                                sti_start_step = 4e4;
                                                                sti_end_step = step_tot;
                                                                rate_ext_on(sti_start_step:sti_end_step) = 1; % duration of stimulus
                                                                rate_ext_n = 0;
                                                                for input_pos = 1:size(qi_pos,1)
                                                                    [qiY,qiX]=ndgrid( wrapToPi(qix- qi_pos(input_pos,1))/qisig, wrapToPi(qix-qi_pos(input_pos,2))/qisig );
                                                                    qiR2 = qiY.*qiY + qiX.*qiX;
                                                                    rate_ext_n = rate_ext_n + qi_contrasts(input_pos)*exp(-0.5*qiR2(:)'); %add gaussian inputs
                                                                end
                                                                rate_ext_n =  rate_ext_E_scale*rate_ext_n;
                                                                writeExtSpikeTinvSettingsHDF5(FID, 1, 1, g_ext,  N_ext, rate_ext_n, rate_ext_on); %run this again to add more inputs
                                                                disp('Custom input has been written.')
                                                                % write synapse para
                                                                writeSynParaHDF5(FID, 'tau_decay_AMPA',5.8,'tau_decay_GABA', 6.5);
                                                                
                                                                % inhibitory STDP
                                                                if inh_STDP == 1
                                                                    writeInhSTDPHDF5(FID, 2, 1, 1*sec);
                                                                end
                                                                
                                                                if STD_on == 1
                                                                    % writeSTD(FID, 1, 1, 5*sec);
                                                                    writeSTDHDF5(FID, 1, 1, 1);
                                                                end
                                                                
                                                                % write runaway killer
                                                                min_ms = 500; % 5 sec
                                                                runaway_Hz = 100; % ??
                                                                Hz_ms = 200; % ms
                                                                writeRunawayKillerHDF5(FID, 1, min_ms, runaway_Hz, Hz_ms);
                                                                
                                                                % random initial condition settings (int pop_ind, double p_fire)
                                                                p_fire = [0.1 0.00]; % between [0,1], 0.05
                                                                writeInitVHDF5(FID, p_fire); 
                                                                
                                                                %%%%%%%%%%%%%%%%%%% Chemical Connections %%%%%%%%%%%%%%%%%%%%%%%
                                                                % type(1:AMAP, 2:GABAa, 3:NMDA)
                                                                
                                                                Lattice_I = quasi_lattice_2D( N(2) , hw);
                                                                %%%%%%%%%%%%%%%%%%%%%%
                                                                % generate in- and out-degree accroding to (1) distance-dependent rule
                                                                % (2) degree distributios and (3) common neighbour rule
                                                                deg_mean = N_e*P0_init;
                                                                deg_std_logn = degree_CV*deg_mean;
                                                                
                                                                [ deg_in_0, deg_out_0 ] = hybrid_degree( N_e, deg_mean, deg_std_logn, in_out_r, deg_hybrid );
                                                                
                                                                [ I_ee, J_ee, dist_IJ, iter_hist, Lattice_E ] = generate_IJ_2D( deg_in_0, deg_out_0, tau_c_EE, cn_scale_wire, iter_num );
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
                                                                writeChemicalConnectionHDF5(FID, Type_mat(1, 1),  1, 1,   I_ee,J_ee,K_ee,D);
                                                                clear I J K D;
                                                                
                                                                [~,ind_sorted] = sort(in_degree);
                                                                %                                                                 sample_neuron = ind_sorted(1:500:end);
                                                                
                                                                %%%%%%%%%%%%%%%%%%%%%%
                                                                [ I,J ] = Lattice2Lattice( Lattice_I, Lattice_E, hw, tau_c_I, P_mat(2,1) );
                                                                D = rand(size(I))*delay;
                                                                K = zeros(size(J));
                                                                for i_E = 1:N(1)
                                                                    mu_K_tmp = EE_input(i_E)/sum(J==i_E)*(g_EI/g_mu)*g_balance;
                                                                    K(J==i_E) = abs(randn([1 sum(J==i_E)])*(mu_K_tmp/4) + mu_K_tmp); % this is a bit too arbitary!
                                                                end
                                                                % K = ones(size(I))*K_mat(2,1);
                                                                writeChemicalConnectionHDF5(FID, Type_mat(2, 1),  2, 1,   I,J,K,D);
                                                                clear I J K D;
                                                                
                                                                %%%%%%%%%%%%%%%%%%%%%%
                                                                [ I,J ] = Lattice2Lattice( Lattice_E, Lattice_I, hw, tau_c_IE, P_mat(1,2) );
                                                                D = rand(size(I))*delay;
                                                                K = ones(size(I))*K_mat(1,2);
                                                                writeChemicalConnectionHDF5(FID, Type_mat(1, 2),  1, 2,   I,J,K,D);
                                                                clear I J K D;
                                                                
                                                                %%%%%%%%%%%%%%%%%%%%%%
                                                                [ I,J ] = Lattice2Lattice( Lattice_I, Lattice_I, hw, tau_c_I, P_mat(2,2) );
                                                                D = rand(size(I))*delay;
                                                                K = ones(size(I))*K_mat(2,2);
                                                                writeChemicalConnectionHDF5(FID, Type_mat(2, 2),  2, 2,   I,J,K,D);
                                                                clear I J K D;
                                                                
                                                                
                                                                %%%%%%% data sampling
                                                                sample_pop = 1;
                                                                writePopStatsRecordHDF5(FID, sample_pop);
                                                                for pop_ind_pre = 1:Num_pop
                                                                    pop_ind_post = sample_pop;
                                                                    if pop_ind_pre == Num_pop
                                                                        syn_type = 2;
                                                                    else
                                                                        syn_type = 1;
                                                                    end
                                                                    %writeSynSampling(FID, pop_ind_pre, pop_ind_pos228t, syn_type, sample_neurons, sample_steps)
                                                                    writeSynStatsRecordHDF5(FID, pop_ind_pre, pop_ind_post, syn_type)
                                                                end
                                                                %traditional
                                                                sample_neuron = 1:N(1);
                                                                sample_steps = zeros(1,step_tot);
                                                                sample_steps(1:50:end) = 1;                                                       
                                                             
                                                                % Add LFP sampling
                                                                [Lattice, ~] = lattice_nD(2, hw);
                                                                LFP_neurons = [];
                                                                
                                                                LFP_centre_x = linspace(-hw,hw,41); % E16:9  E100:21 E400:41 E1600:81
                                                                LFP_centre_y = linspace(-hw,hw,41);
                                                                LFP_centre_x = LFP_centre_x(2:2:40); % E16(2:2:8)  E100(2:2:20) E400(2:2:40) E1600(2:2:80)
                                                                LFP_centre_y = LFP_centre_y(2:2:40);
                                                                [LFP_centre_x, LFP_centre_y] = meshgrid(LFP_centre_x, LFP_centre_y);
                                                                LFP_centre_x = LFP_centre_x(:);
                                                                LFP_centre_y = LFP_centre_y(:);
                                                                %LFP_centre_x = [0,31]; % one in centere, another in corner
                                                                %LFP_centre_y = [0,31];
                                                                for cc = 1:length(LFP_centre_x)
                                                                    dist = lattice_nD_find_dist(Lattice, hw, LFP_centre_x(cc) , LFP_centre_y(cc));
                                                                    gaus_tmp = 1/(LFP_range_sigma*sqrt(2*pi))*exp(-0.5*(dist/LFP_range_sigma).^2) .* double(dist <= LFP_range_sigma*2.5);
                                                                    LFP_neurons = [LFP_neurons; transpose(gaus_tmp(:))]; %#ok<AGROW>
                                                                end
                                                                
                                                                writeLFPRecordHDF5(FID, 1, LFP_neurons);
                                                                
                                                                % Explanatory (ExplVar) and response variables (RespVar) for cross-simulation data gathering and post-processing
                                                                % Record explanatory variables, also called "controlled variables"
                                                                
                                                                writeExplVarHDF5(FID, 'discard_transient', discard_transient, ...
                                                                    'loop_num', loop_num, ...
                                                                    'delay', delay, ...
                                                                    'g_mu', g_mu, ...
                                                                    'g_EI', g_EI, ...
                                                                    'g_IE', g_IE, ...
                                                                    'g_II', g_II, ...
                                                                    'g_ext', g_ext,...
                                                                    'SpikeFreqAapt', SpikeFreqAapt, ...
                                                                    'rate_ext_E', rate_ext_E,...
                                                                    'rate_ext_I', rate_ext_I,...
                                                                    'P0_init', P0_init, ...
                                                                    'degree_CV', degree_CV,...
                                                                    'in_out_r', in_out_r, ...
                                                                    'tau_c_EE', tau_c_EE, ...
                                                                    'tau_c_IE', tau_c_IE, ...
                                                                    'tau_c_I', tau_c_I, ...
                                                                    'STD_on', STD_on, ...
                                                                    'cn_scale_wire', cn_scale_wire, ...
                                                                    'cn_scale_weight', cn_scale_weight, ...
                                                                    'mu_p', mu_p,...
                                                                    's_p', s_p, ...
                                                                    'inh_STDP', inh_STDP, ...
                                                                    'deg_hybrid', deg_hybrid,...
                                                                    'LFP_range_sigma', LFP_range_sigma,...
                                                                    'sti_start_step', sti_start_step,...
                                                                    'sti_end_step',sti_end_step,...
                                                                    'g_balance',g_balance);
                                                                
                                                                
                                                                %                                                     % Adding comments in raster plot
                                                                %                                                     comment1 = ' '; %'p=[p0_init 0.3 0.3 0.3], k = ?, tau_decay_GABA=3';
                                                                %                                                     comment2 = datestr(now,'dd-mmm-yyyy-HH:MM');
                                                                %                                                     writeExplVarHDF5(FID, 'comment1', comment1, 'comment2', comment2);
                                                                
                                                                % save in_degree and sample neuron data based on in_degree
                                                                save([FID(1:end-6),...
                                                                    '_config_data.mat'], 'in_degree', 'out_degree', 'dist_IJ', 'iter_hist','K_ee_mean', 'LFP_centre_x', 'LFP_centre_y', 'Lattice_I', 'Lattice_E');
                                                                
                                                                % append this file self into .ygin for future reference
                                                                appendThisMatlabFileHDF5(FID)
                                                                
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
    end
end
end




% This function must be here!
function appendThisMatlabFileHDF5(FID)

% need to coopy and past the following code into the file to be appended!
text = fileread([mfilename('fullpath'),'.m']);
hdf5write(FID,'/config/MATLAB/config.m',text,'WriteMode','append');

end