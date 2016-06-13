function main_demo()
% Use coherent units (msec+mV+nF+miuS+nA) unless otherwise stated


% Seed the Matlab random number generator
seed = 1;
% Open new (uniquely time-stamped) input files for the SpikeNet C++
% simulator.
% FID is the main input file, FID_syn is the input file with the synaptic
% connectivity definitions (could be very large).
[FID, FID_syn] = new_ygin_files_and_randseed(seed);



% Time step (ms)
dt = 0.1; % 0.1 is usually small enough
% Total number of simulation steps
step_tot = 10^2;
% Create two neuron populations with 50 neurons in the 1st population and
% 10 neurons in the 2nd population.
N = [40, 10];
% Write the above basic parameters to the input file
writeBasicPara(FID, dt, step_tot, N);


% Define non-parameters for the neuorn model
% For demo purpose, the values used as following are the still the
% default ones.
% If default values are to be used, there is no need to re-define them.
Cm = 0.25; % (nF) membrane capacitance
tau_ref = 2.0; % (ms) absolute refractory time
V_rt = -60.0;   % (mV) reset potential
V_lk = -70.0;   % (mV) leaky reversal potential
V_th = -50.0;   % (mV) firing threshold
g_lk = 0.0167;   % (miuS) leaky conductance (note that Cm/gL=15 ms)
for pop_ind = 1:2
    writePopPara(FID, pop_ind,'Cm',Cm,'tau_ref',tau_ref,'V_rt',V_rt,'V_lk',V_lk,'V_th',V_th,'g_lk',g_lk);
end

% Add spike-frequency adaptation to the 1st population
pop = 1;
writeSpikeFreqAdpt(FID, pop);


% Define the initial condition
p_fire = [0.1 0.1]; % initial firing probabilities for both populations
r_V0 = [0.5 0.5]; %set initial V distribution to be [V_rt, V_rt + (V_th-V_rt)*r_V0] 
writeInitCond(FID, r_V0, p_fire)

% Add external Poissonian spike trains into the 1st population
pop = 1;
type_ext = 1; % 1 for AMPA-like syanpses
g_ext = 2*10^-3; % synaptic coupling strength (miuS)
N_ext = 1000; % number of independent external connections onto each neuron
rate_ext = 2*ones(1, step_tot); % Poisson rate defiend for each time step (Hz) 
N_start = 1;
N_end = 10; % the external spikes trains are only added to neuron No.1-No.10.
writeExtSpikeSettings(FID, pop, type_ext, g_ext,  N_ext, rate_ext,  N_start, N_end );


% Add external currents (Gaussian white noise) to the 2nd population
pop = 2;
I_ext_mean = 0.5*ones(1,N(2)); % defined for each neuron (nA)
I_ext_std = 0.2*ones(1,N(2)); % defined for each neuron (nA)
writeExtCurrentSettings(FID, pop, I_ext_mean, I_ext_std)





% Define runaway killer
% The computational cost of the simulation is directly proportional to the
% avearage firing rate of the network. Very often, your parameters may lead
% to biologically unrealistically high firing rate or even runaway
% activity (the highest rate possible ~1/tau_ref). To save the computational
% resources for more interesting simulation cases, it is advised to define
% the runawary killer.
% Note that the simulator will still output all the relevant data before
% the kill.
min_ms = 500; % the mininum simulation duration that should be guaranteed
runaway_Hz = 40; % the threshold above which the simulation should be killed
Hz_ms = 200; % the window length over which the firing rate is averaged
pop = 1; % the population to be monitored by the runaway killer
writeRunawayKiller(FID, pop, min_ms, runaway_Hz, Hz_ms);


% Record the basic statistics for the 1st neuron population
pop = 1;
writePopStatsRecord(FID, pop);

% Define AMAP-like excitatory synaptic connection wtihin the 1st population
pop_pre = 1;
pop_post = 1;
syn_type = 1; % 1 for AMPA-like synapse
g_11 = 10*10^-3; % synaptic coupling strength (miuS)
A = rand(N(1), N(1)) < 0.5; % random adjaceny matrix with a wiring probability of 0.5
[I_11, J_11] = find(A);
K_11 = g_11*ones(size(I_11)); % identical coupling strength
D_11 = rand(size(I_11))*5; % uniformly random conduction delay between 0 and 5 ms
writeChemicalConnection(FID_syn, syn_type,  pop_pre, pop_post, I_11, J_11, K_11, D_11);


% Record the basic statistics for the excitatory synapses within the 1st
% neuron population
pop_pre = 1;
pop_post = 1;
syn_type = 1; % 1 for AMPA-like synapse
writeSynStatsRecord(FID, pop_pre, pop_post, syn_type);

% Sample detailed time serious data from the 1st population population
pop = 1;
sample_neuron = 1:10:N(1); % sample every 10th neuron
sample_steps = 1:2:step_tot; % sample every 2nd time step
sample_data_type = [1,1,1,1,0,0,1,0]; % [V,I_leak,I_AMPA,I_GABA,I_NMDA,I_GJ,I_ext, I_K]
writeNeuronSampling(FID, pop, sample_data_type, sample_neuron, sample_steps )

% Record explanatory variables that are nacessary for post-processing
discard_transient = 0; % transient period data to be discarded (ms)
writeExplVar(FID, 'discard_transient', discard_transient, ...
    'g_11', g_11, ...
    'g_ext', g_ext, ...
    'N_ext', N_ext, ...
    'rate_ext', rate_ext);


% Additional comments to be passed to the post-processing stage
comment1 = 'This is a tutorial case.';
comment2 = 'If you have any question regarding SpikeNet, please contact Yifan Gu (yigu8115@gmail.com).';
writeExplVar(FID, 'comment1', comment1, 'comment2', comment2);


% parameter

STD_on = 0;
[ fit_g_2_EPSP_2, ~ ] = g_EPSP_conversion( );
EPSP_mu = fit_g_2_EPSP_2(g_mu);
EPSP_sigma = 1;



% for inh_STDP = [0 ];

% write synapse para
writeSynPara(FID, 'tau_decay_GABA', 3);

% inhibitory STDP
writeInhSTDP(FID, 2, 1, 1*sec);

% writeSTD(FID, 1, 1, 5*sec);
writeSTD(FID, 1, 1, 1);







%%%%%%%%%%%%%%%%%%% Chemical Connections %%%%%%%%%%%%%%%%%%%%%%%
% type(1:AMAP, 2:GABAa, 3:NMDA)





% append this matlab file to the input file for future reference
appendThisMatlabFile(FID)

end



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


