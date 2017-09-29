function james_adamreprod()
% Use coherent units (msec+mV+nF+miuS+nA) unless otherwise stated
% The following configuration reproduced the dynamics in the paper by Adam
% Kean and Pulin Gong, Propagating Waves Can Explain Irregular Neural
% Dynamics.
% For weight conversion, please see conductance_model_comparison()

%%%% Seed the Matlab random number generator
seed = 1;

%%%% Open new (uniquely time-stamped) input files for the SpikeNet C++
% simulator.
% FID is the main input file, FID_syn is the input file with the 
% synaptic connectivity definitions (could be very large).
[FID] = new_ygin_files_and_randseedHDF5(seed);
% If no FID_syn is needed, use FID = new_ygin_files_and_randseed(seed,0)

% Use Adam 2016 synapse model instead of the default model
model_choice = 2;
writeSynapseModelChoiceHDF5(FID, model_choice)

%%%% Define some basic parameters
% Time step (ms)
dt = 0.05; % 0.1 is usually small enough
% Total number of simulation steps
step_tot = 10*10^3;

%Define grids [no rows, no columns, grid step size]
Grid(1,:)=[50,50,1];
Grid(2,:)=[25,25,2];


% Drange=[10,10;
%     15, 15]; % maximum connection range
% sigma=[12,12;
%     1000,1000]; % connection streength spatial Gaussain 
% W=[0.23,0.23;
%     0.3,0.3;]; % max connection strength (muS)
% 
% % weight conversion
% W(1,:) = W(1,:);
% W(2,:) = W(2,:);
% 
% pbc=1; % periodic boundary conditions
% 
% F=[15,2]*10^-3; % muS

Drange=[15,15;
    15, 15]; % maximum connection range
sigma=[12,12;
    90,90]; % connection streength spatial Gaussain 
W=[0.24,0.24;
    0.79,0.79;]; % max connection strength (muS)

pbc=1; % periodic boundary conditions

F=[15,0]*10^-3; % muS



SynapseType=[1 1;
    2 2];

N = [Grid(1,1)*Grid(1,2), Grid(2,1)*Grid(2,2)];
% Write the above basic parameters to the input file
writeBasicParaHDF5(FID, dt, step_tot, N);


%%%% Define non-parameters for the neuorn model
% For demo purpose, the values used as following are the still the
% default ones.
% If default values are to be used, there is no need to re-define them.
Cm = 1; % (nF) membrane capacitance
tau_ref = 5.0; % (ms) absolute refractory time
V_rt = -70.0;   % (mV) reset potential
V_lk = -70.0;   % (mV) leaky reversal potential
V_th = -55.0;   % (mV) firing threshold
g_lk = 0.050;   % (muS) leaky conductance 
for pop = 1:length(N)
    writePopParaHDF5(FID, pop,'Cm',Cm,'tau_ref',tau_ref,'V_rt',V_rt,...
        'V_lk',V_lk,'V_th',V_th,'g_lk',g_lk);
end

%%%% Add spike-frequency adaptation to the 1st population



%%%% Define the initial condition
p_fire = 0*ones(1,length(N)); % initial firing probabilities for both populations
% set initial V distribution to be [V_rt, V_rt + (V_th-V_rt)*r_V0] 
r_V0 = 1*ones(1,length(N));
writeInitCondHDF5(FID, r_V0, p_fire)


%%%% Define runaway killer
% The computational cost of the simulation is directly proportional to 
% the avearage firing rate of the network. Very often, your parameters 
% may lead to biologically unrealistically high firing rate or even 
% runaway activity (the highest rate possible ~1/tau_ref). To save the 
% computational resources for more interesting simulation cases, it is 
% advised to define the runawary killer.
% Note that the simulator will still output all the relevant data before
% the kill.
min_ms = 500; % the min simulation duration that should be guaranteed
runaway_Hz = 40; % the threshold above which the simu should be killed
Hz_ms = 200; % the window length over which the firing rate is averaged
pop = 1; % the population to be monitored by the runaway killer
writeRunawayKillerHDF5(FID, pop, min_ms, runaway_Hz, Hz_ms);


%%%% Record the basic statistics for the 1st neuron population
pop = 1;
writePopStatsRecordHDF5(FID, pop);

%%%%%%%%%%%%%%%%%%% Chemical Connections %%%%%%%%%%%%%%%%%%%%%%%
% type(1:AMAP, 2:GABAa, 3:NMDA)

% Define connectivity between populations
for pop_pre=1:length(N)
    for pop_post=1:length(N)
        syn_type = SynapseType(pop_pre,pop_post); % 1 for AMPA-like synapse
        % Define COnnectivity matrix
        DM=GridDM(Grid(pop_pre,:),Grid(pop_post,:),pbc);
        A=(DM<=Drange(pop_pre,pop_post)).*W(pop_pre,pop_post).*exp(-DM.^2/sigma(pop_pre,pop_post));
        if pop_pre==pop_post
           A(1:N(pop_pre)+1:N(pop_pre)*N(pop_pre)) = 0; %remove self conns 
        end
        [I, J, K] = find(A);
        D = rand(size(I))*0; % uniformly random conduction delay ms
        writeChemicalConnectionHDF5(FID, syn_type,  pop_pre, pop_post, ...
            I, J, K, D);
    end
end

%External Conductances
for pop =1:length(N)
    %uniform 
    F_ext=F(pop)*ones(1,N(pop));
    writeExtConductanceSettingsHDF5(FID, pop, F_ext, 0*ones(1,N(pop)) );
end
%Gaussian
% pop=1
% [x,y]=meshgrid(1:Grid(pop,1),1:Grid(pop,2))
% sigmag=2;
% g_Gaus_max=g_ext_mean(pop);
% g_Gaus=g_Gaus_max*exp(-((x-(Grid(pop,1)-1)/2)^2+(y-((Grid(pop,2)-1)/2))^2)/sigmaG)
% g_Gaus=reshape(g_Gauss,N(pop),1)    
% writeExtConductanceSettings(FID, pop, g_Gaus, 0);

%%%%%%%%%%%%%%%%%%% Chemical Connections Done %%%%%%%%%%%%%%%%%%%%%%%


%%%% Use an existing synapse connectivity definition file instead of
% generating a new one:
% writeSynFilename(FID, 'path/to/existing/XXX.ygin_syn')
% In this case, you should use 
% "FID = new_ygin_files_and_randseed(seed, 0)"
% to avoid creating an empty ygin_syn file.



% %%%% Define non-default synapse parameters
% writeSynPara(FID, 'tau_decay_GABA', 7, 'Dt_trans_AMPA' ,0.5,...
% 'tau_decay_AMPA', 2.0, 'Dt_trans_GABA', 0.5)

%%%% synapse parameters conversion
writeSynParaHDF5(FID, 'tau_decay_GABA', 7, 'Dt_trans_GABA', 0.5,...
'tau_decay_AMPA', 2.0, 'Dt_trans_AMPA' , 0.5)

% "help writeSynPara" to see all the parameter names and default values


%%%% Record the basic statistics for the excitatory synapses within 
% the 1st neuron population
pop_pre = 1;
pop_post = 1;
syn_type = 1; % 1 for AMPA-like synapse
writeSynStatsRecordHDF5(FID, pop_pre, pop_post, syn_type);

% %%%% Sample detailed time serious data from the 1st population
% pop = 1;
% sample_neuron = 1:500:N(1); % sample every 10th neuron
% sample_steps = 1:step_tot; % sample every 2nd time step
% sample_data_type = logical([1,1,1,1,0,0,1,0]); 
% % The logical vector above corresponds to 
% % [V,I_leak,I_AMPA,I_GABA,I_NMDA,I_GJ,I_ext, I_K]
% writeNeuronSampling(FID, pop, sample_data_type, ...
%     sample_neuron, sample_steps)

%%%% Record explanatory variables that are nacessary for post-processing
discard_transient = 0; % transient period data to be discarded (ms)
writeExplVarHDF5(FID, 'discard_transient', discard_transient, ...
    'F', F);


%%%% append this matlab file to the input file for future reference
appendThisMatlabFileHDF5(FID)

end



% This function must be here!
function appendThisMatlabFileHDF5(FID)

% need to coopy and past the following code into the file to be appended!
text = fileread([mfilename('fullpath'),'.m']);
hdf5write(FID,'/config/MATLAB/config.m',text,'WriteMode','append');

end


function DM=GridDM(InGrid, OutGrid, pbc)
% Calculates the distance between all pairs of points from one grid to
% another.
% Grid is a vector of [ no rows, no columns, grid step size]
NoIn=InGrid(1)*InGrid(2);
NoOut=OutGrid(2)*OutGrid(2);

DM=zeros(NoIn,NoOut);

for i=0:NoIn-1
    % Convert 1D index to 2D spat coords
    yi=floor(i/InGrid(1));
    xi=i-yi*InGrid(1);
    
    % Scale for grid step size
    yi=yi*InGrid(3);
    xi=xi*InGrid(3);
    
    for j=0:NoOut-1
        % Convert 1D index to 2D spat coords
        yj=floor(j/OutGrid(1));
        xj=j-yj*OutGrid(1);

        % Scale for grid step size
        yj=yj*OutGrid(3);
        xj=xj*OutGrid(3);
        
        delx=abs(xi-xj);
        dely=abs(yi-yj);
        
        if pbc
            %If grid has periodic boundaries then check to see which
            % distance is shortest (i.e via the periodic boundary or not
            delx=min([delx,abs(InGrid(2)*InGrid(3)-delx)]);
            dely=min([dely,abs(InGrid(1)*InGrid(3)-dely)]);
        end
        
        DM(i+1,j+1)=sqrt(delx^2+dely^2);
        
    end 
end
end

