function generate_learning_config()
tic
% Use units msec+mV+nF+miuS+nA unless otherwise stated

format long

%%%% Seed the Matlab random number generator
seed = 1;

%%%% Open new (uniquely time-stamped) input files for the SpikeNet C++
% simulator.
% FID is the main input file, FID_syn is the input file with the 
% synaptic connectivity definitions (could be very large).
FID = new_ygin_files_and_randseedHDF5(seed);
% If no FID_syn is needed, use FID = new_ygin_files_and_randseed(seed,0)

% Choose synapse model
model_choice = 2;
writeSynapseModelChoiceHDF5(FID, model_choice);

%%%% Define some basic parameters
% Time step (ms)
dt = 0.1; 
% Total number of simulation steps
step_tot = 10^6;

% Choose neuron type for each population (0 = excitatory, 1 = inhibitory)
neu_types=[0,0,1];

%Define 2D neuron grids [no rows, no columns, grid step size]
gsize=100; % the defines the size of the excitatory grid of neurons
Inhscale=2; % inhibitory neurons are spaced at twice separation of excitatory neurons
Stimgrid=128; % Change this to suit your stimulus grid side length!!!
Grid(1,:)=[Stimgrid,Stimgrid,(gsize-1)/(Stimgrid-1),(gsize-1)/(Stimgrid-1)]; % this specifies the dimensions of the stimulus input. 
Grid(2,:)=[gsize,gsize,1,1]; % excitatory population
Grid(3,:)=[gsize/Inhscale,gsize/Inhscale,Inhscale,Inhscale]; % inhibitory population

npops=size(Grid,1); % number of populations

% Number of neurons in each population
N=zeros(npops,1);
for i=1:npops
    N(i) = Grid(i,1)*Grid(i,2);
end

%%%% Define parameters for the neuron model

Cm = 1; % (nF) membrane capacitance
tau_ref = 5.0; % (ms) refractory time
V_rt = -70.0;   % (mV) reset potential
V_lk = -70.0;   % (mV) leaky reversal potential
V_th = -55.0;   % (mV) firing threshold
g_lk = 0.05;   % (muS) leaky conductance 
for pop = 1:length(N)
    writePopParaHDF5(FID, pop,'Cm',Cm,'tau_ref',tau_ref,'V_rt',V_rt,...
        'V_lk',V_lk,'V_th',V_th,'g_lk',g_lk);
end

% Setup parameters of the learning/plasticity

tau=1/((g_lk/Cm)*dt) % neuronal time constant in units of timesteps
tau_A=7*tau
Vspike=1*abs(V_th-V_lk); % voltage change from resting to threshold

dropout=[0,0.75,0.5];    % noise applied to learning for spikes in each population

gamma=10^-3*[0,0.683,2.73;
                            0, 1.67,6.67;
                            0, 2.78,1.11];


epsilon=10^-6*[0, 1.25, 1.25;
                            0, 0.512,0.512;
                            0,1.13,1.13];

%%%% Now define some connectivity parameters

% the maximum connection range of neurons between populations, in grid units
Drange=[15 15 15;
                 15 15 15;
                 15 15 15]; 
             
% Initial connection strengths between populations
W= 10^-2*[0, 2.2, 2.2;
                    0, 1.67*10^-3, 4.17*10^-4;
                    0, 2.2*10^-2, 5.5*10^-3]; 

% no periodic boundaries
pbc=0; 

% Define synapse types between populations 1 = excitatory, 2 = inhibitory
SynapseType=[0 1 1;
             0 1 1;
             0 2 2]; 

% Write the above basic parameters to the input file
writeBasicParaHDF5(FID, dt, step_tot, N);


%%%% Define the initial condition
p_fire = 0*ones(1,length(N)); % initial firing probabilities for populations
% set initial V distribution to be [V_rt, V_rt + (V_th-V_rt)*r_V0] 
r_V0 = 1*ones(1,length(N));
writeInitCondHDF5(FID, r_V0, p_fire)

% Define the stimulus input filename
% the stimulus input file should be and HDF5 file containing the following datasets
% t : a list of integer times (in microseconds) for each spike in the stimulus
% x: a list of integer x grid coordinates for each spike in the stimulus
% y: a list of integer y grid coordinates for each spike in the stimulus
% max_x: an integer specifying the size of the stimulus grid in the x coordinate
% max_y: an integer specifying the size of the stimulus grid in the y coordinate
writeSpikeFileInputPopHDF5(FID,'!!!ADD_YOUR_FILENAME_HERE!!!.h5',1);
% include the full path if not in the same directory as the simulator
% We have used DVS recordings. Some are available here 
% from Prof. Tobi Delbruck:  https://www.ini.uzh.ch/~tobi/dvs/
% These will need to be preprocessed into an HDF5 file as described above


% Define connectivity between populations
for pop_pre=1:length(N)
    [xpre,ypre]=meshgrid(1:Grid(pop_pre,3):1+Grid(pop_pre,3)*(Grid(pop_pre,1)-1),1:Grid(pop_pre,4):1+Grid(pop_pre,4)*(Grid(pop_pre,2)-1));
    xpre=reshape(xpre,N(pop_pre),1);
    ypre=reshape(ypre,N(pop_pre),1);
    for pop_post=1:length(N)
        if W(pop_pre,pop_post) ~= 0
            [xpost,ypost]=meshgrid(1:Grid(pop_post,3):1+Grid(pop_post,3)*(Grid(pop_post,1)-1),1:Grid(pop_post,4):1+Grid(pop_post,4)*(Grid(pop_post,2)-1));

            xpost=reshape(xpost,N(pop_post),1);
            ypost=reshape(ypost,N(pop_post),1);
            
            maxGridx=gsize;
            maxGridy=gsize;

                len=round(2*pi*Drange(pop_pre,pop_post)^2/min([N(pop_pre),N(pop_post)])*length(xpre)*length(xpost));
                si=-1*ones(1,len);
                sj=-1*ones(1,len);;
                sD=-1*ones(1,len);;
                ctr=1;
                for i=1:length(xpre)
                    i/length(xpre);
                    for j=1:length(xpost)
                        tempx=abs(xpre(i)-xpost(j));
                        tempy=abs(ypre(i)-ypost(j));
                        if pbc
                            tempx=min(tempx,abs(maxGridx-tempx));
                            tempy=min(tempy,abs(maxGridy-tempy));
                        end
                        tempd=sqrt(tempx^2+tempy^2);
                        if (tempd<Drange(pop_pre,pop_post))&&~((pop_pre==pop_post)&&(i==j))
                            si(ctr)=i;
                            sj(ctr)=j;
                            sD(ctr)=tempd;
                            ctr=ctr+1;
                        end
                    end
                end
                ind=find(si>0);
                I=si(ind);
                J=sj(ind);
                K=W(pop_pre,pop_post)*ones(size(sj(ind)));

            
            
            D = rand(size(I))*0; % uniformly random conduction delay ms
            writeChemicalConnectionHDF5(FID, SynapseType(pop_pre,pop_post),  pop_pre, pop_post, ...
                I, J, K, D);

                writeSynJHLearnHDF5(FID, pop_pre, pop_post, tau_A,epsilon(pop_pre,pop_post),0,0, 1,tau, gamma(pop_pre,pop_post),0,neu_types(pop_pre),neu_types(pop_post), 1);
                writeJHLearnDirectionHDF5(FID, pop_pre, pop_post,0);

        end
    end
end
% Population learning noise
for pop =1:length(N)
    writePopJHLearnHDF5(FID, pop, dropout(pop), dropout(pop));
end

%%%%%%%%%%%%%%%%%%% Chemical Connections Done %%%%%%%%%%%%%%%%%%%%%%%

%%%% synapse parameters conversion
writeSynParaHDF5(FID, 'tau_decay_GABA', 7, 'Dt_trans_AMPA' , 0.5,...
    'tau_decay_AMPA', 2.0, 'Dt_trans_GABA', 0.5)


%%%% Sample detailed time series data
for pop = 1:npops
    sample_neuron =round(linspace(1,N(pop),N(pop))); % sample every ?th neuron
    sample_steps = zeros(1, step_tot); sample_steps(max(1,step_tot-1e4):20:end) = 1; % sample every ? time step
    sample_data_type = [1,0,1,1,0,0,0,0,1];
    % The logical vector above corresponds to
    % [V,I_leak,I_AMPA,I_GABA,I_NMDA,I_GJ,I_ext, I_K,Q]
    writeNeuronSamplingHDF5(FID, pop, sample_data_type, ...
        sample_neuron, sample_steps)
end

%%%% append this matlab file to the input file for future reference
appendThisMatlabFileHDF5(FID)
toc
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

function appendThisMatlabFileHDF5(FID)

text = fileread([mfilename('fullpath'),'.m']);
hdf5write(FID,['/config/MATLAB/config.m'],text,'WriteMode','append');

end
