function main_GU_et_al_2018(varargin)
% The simulation configuration file for GU et al., 2018
% All parameter values are for the default system (the working point)


%%%%% simulation time
dt = 0.1; % simulation time step
sec = round(10^3/dt); % number of time steps for 1 second
step_tot = 10*sec; % total number of simulation time steps



%%%%% neuron and synapse parameters
tau_ref = 4; % (ms) refractory period
delay_max = 4; % (ms)  the range of conduction delay: [0, delay_max]



%%%%% connectivity 
P_mat = [0.16 0.2; 0.2  0.4]; % P_mat(i,j) gives the mean connection probility between pre-synaptic population i and post-synaptic population j
hw = 31; % half-width of the square grid for excitatory population, hw = 31 gives (31*2+1)^2 = 3969 ~ 4000 excitatory neurons
N_e = (hw*2+1)^2; % size of the excitatory population
N_i = 1000; % size of the inhibitory population

N = [N_e, N_i];
Num_pop = length(N);
Type_mat = ones(Num_pop);
Type_mat(end, :) = 2;
in_out_r =  0.13; % correlation coefficent between the in-degree and out-degree of excitatory neurons
dg_K = 0.01; % (muS) unit increase in potassium conductance (strength of spike frequency adaptation)


a_Gamma = 2; % common neighbor coefficent

iter_num = 5;




g_ext = 2*10^-3; % (muS) strength of external (excitatory) connections



q = 0.4; % hybrid parameter of degree distribution
degree_CV = 0.2; % 0.2 works

g_EE_mu = 4*10^-3; % (muS) mean value of E-to-E connection strength

[ fit_g_2_EPSP, ~ ] = g_EPSP_conversion(); % convertion function
EPSP_mu = fit_g_2_EPSP(g_EE_mu); % (mV) mean value of E-to-E connection strength converted to unit EPSP from resting potenital
EPSP_sigma = 1; % (mV) SD of E-to-E connection strength






g_EI_mu = 13.5*10^-3; % (muS) mean value of I-to-E connection strength
g_IE = 5*10^-3; % (muS) value of E-to-I connection strength
g_II = 25*10^-3; % (muS) value of I-to-I connection strength

rate_ext_I = 1; % (kHz) firing rate of  external input into excitatory neurons
rate_ext_E = 0.85; % (kHz) firing  rate of  external input into inihibitory neurons
tau_c_EE = [8]
tau_c_IE = 10;
tau_c_I = [20]


% seed the matlab rand function
[FID] = new_ygin_files_and_randseedHDF5(loop_num);


% write basic parameters
writeBasicParaHDF5(FID, dt, step_tot, N);

K_mat = [NaN  g_IE;
    g_EI_mu  g_II]; % miuSiemens


writeSpikeFreqAdptHDF5(FID, 1, dg_K);


% write pop para
for pop_ind = 1:Num_pop
    writePopParaHDF5(FID, pop_ind,  'tau_ref', tau_ref);
end

% write external currents
writeExtSpikeSettingsHDF5(FID, 1, 1, g_ext,  N_ext, rate_ext_E*ones(1, step_tot),  ones(1, N(1)) );
writeExtSpikeSettingsHDF5(FID, 2, 1, g_ext,  N_ext, rate_ext_I*ones(1, step_tot),  ones(1, N(2)) );

% write synapse para
writeSynParaHDF5(FID, 'tau_decay_GABA', 3);


% random initial condition settings
writeInitVHDF5(FID, [0.1 0.00]);



%%%%%%%%%%%%%%%%%%% Chemical Connections %%%%%%%%%%%%%%%%%%%%%%%
% type(1:AMAP, 2:GABAa, 3:NMDA)

Lattice_I = quasi_lattice_2D( N(2) , hw);
%%%%%%%%%%%%%%%%%%%%%%
% generate in- and out-degree accroding to (1) distance-dependent rule
% (2) degree distributios and (3) common neighbour rule
deg_mean = N_e*P_mat(1,1);
deg_std_logn = degree_CV*deg_mean;

[ deg_in_0, deg_out_0 ] = hybrid_degree( N_e, deg_mean, deg_std_logn, in_out_r, q );

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


EE_input = full(sum(sparse(I_ee,J_ee,K_ee),1));

D = rand(size(I_ee))*delay_max;
writeChemicalConnectionHDF5(FID, Type_mat(1, 1),  1, 1,   I_ee,J_ee,K_ee,D);
clear I J K D;



%%%%%%%%%%%%%%%%%%%%%%
[ I,J ] = Lattice2Lattice( Lattice_I, Lattice_E, hw, tau_c_I, P_mat(2,1) );
D = rand(size(I))*delay_max;
K = zeros(size(J));
for i_E = 1:N(1)
    mu_K_tmp = EE_input(i_E)/sum(J==i_E)*(g_EI_mu/g_EE_mu);
    K(J==i_E) = abs(randn([1 sum(J==i_E)])*(mu_K_tmp/4) + mu_K_tmp); % this is a bit too arbitary!
end
% K = ones(size(I))*K_mat(2,1);
writeChemicalConnectionHDF5(FID, Type_mat(2, 1),  2, 1,   I,J,K,D);
clear I J K D;

%%%%%%%%%%%%%%%%%%%%%%
[ I,J ] = Lattice2Lattice( Lattice_E, Lattice_I, hw, tau_c_IE, P_mat(1,2) );
D = rand(size(I))*delay_max;
K = ones(size(I))*K_mat(1,2);
writeChemicalConnectionHDF5(FID, Type_mat(1, 2),  1, 2,   I,J,K,D);
clear I J K D;

%%%%%%%%%%%%%%%%%%%%%%
[ I,J ] = Lattice2Lattice( Lattice_I, Lattice_I, hw, tau_c_I, P_mat(2,2) );
D = rand(size(I))*delay_max;
K = ones(size(I))*K_mat(2,2);
writeChemicalConnectionHDF5(FID, Type_mat(2, 2),  2, 2,   I,J,K,D);
clear I J K D;


%%%%%%% data sampling
writePopStatsRecordHDF5(FID, 1);
writePopStatsRecordHDF5(FID, 2);
for pop_ind_pre = 1:Num_pop
    pop_ind_post = 1;
    if pop_ind_pre == Num_pop
        syn_type = 2;
    else
        syn_type = 1;
    end
    writeSynStatsRecordHDF5(FID, pop_ind_pre, pop_ind_post, syn_type)
end

writeNeuronSamplingHDF5(FID, 2, [1,1,1,1,0,0,1,0], [1 100], ones(1, step_tot) )


sample_neuron = randperm(N(1),20);
writeNeuronSamplingHDF5(FID, 1, [1,1,1,1,0,0,1,1], sample_neuron, ones(1, step_tot) )



writeExplVarHDF5(FID, 'discard_transient', 0);


% append this file self into .ygin for future reference
appendThisMatlabFileHDF5(FID)

disp('Matlab pre-processing done.')


end


% This function must be here!
function appendThisMatlabFileHDF5(FID)

% need to coopy and past the following code into the file to be appended!
text = fileread([mfilename('fullpath'),'.m']);
hdf5write(FID,'/config/MATLAB/config.m',text,'WriteMode','append');

end