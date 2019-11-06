function main_Chen_and_Gong_2019(varargin)

% Loop number for PBS array job
loop_num = 0;
% use /SpikeNet/ForExtCurrent/ImageProcess_DOG.m to generate the external
% stimulus
d_E=dir('*E.h5');
d_I=dir('*I.h5');
for dump = 1:length(d_E)  
    
    ImageName_E=d_E(dump).name;
    ImageName_I=d_I(dump).name;     
    
    % For PBS array job
    loop_num = loop_num + 1;
    if nargin ~= 0
        PBS_ARRAYID = varargin{1};
        if loop_num ~=  PBS_ARRAYID
            continue;
        end
    end
    
    % Use coherent units (msec+mV+nF+miuS+nA) unless otherwise stated
      
    % Open new (uniquely time-stamped) input files for the SpikeNet C++
    % simulator.
    % FID is the main input file, FID_syn is the input file with the
    % synaptic connectivity definitions (could be very large).
    [FID] = new_ygin_files_and_randseedHDF5(loop_num);
    
    % Use Adam 2016 synapse model instead of the default model
    model_choice = 2;
    writeSynapseModelChoiceHDF5(FID, model_choice);
    
    %%%% Define some basic parameters
    % Time step (ms)
    dt = 0.1; % 0.1 is usually small enough
    % Total number of simulation steps
    step_tot = 1*10^5;

    %Define grids [no rows, no columns, grid step size]
    gsize=250;
    Grid = sparse(gsize,gsize);
    Grid(2:2:gsize,2:2:gsize) = true; % 0 is Exc, 1 is Inh.
    
    npops=2;
    N=zeros(npops,1);
    N(1) =3/4*gsize^2;
    N(2) =1/4*gsize^2;
    %%%% Define non-parameters for the neuorn model
    % For demo purpose, the values used as following are the still the
    % default ones.
    % If default values are to be used, there is no need to re-define them.
    Cm = 1; % (nF) membrane capacitance
    tau_ref = 5.0; % (ms) absolute refractory time
    V_rt = -75.6250;   % (mV) reset potential
    V_lk = -70;   % (mV) leaky reversal potential
    V_th = -40;   % (mV) firing threshold
    g_lk = 0.050;   % (muS) leaky conductance
    V_ex= 0; % (mV) exitatory reversal potential
    V_in= -80;% (mV) inhibitatory reversal potential 
    
    % for exponential leaky integrate and fire model    
    ELIF_VT=-60.6250; % mv for ELIF
    ELIF_delT=6.5625; % mv    
    
    for pop = 1:length(N)
        writePopParaHDF5(FID, pop,'Cm',Cm,'tau_ref',tau_ref,'V_rt',V_rt,...
            'V_lk',V_lk,'V_th',V_th,'g_lk',g_lk);
        writeELIFNeuronModelHDF5(FID,pop,ELIF_VT,ELIF_delT) % USE Exponential LIF as in Brandons work
    end
    
    % external conductance
    F=10^-2; % muS    
    F_ext=F*ones(1, N(1));
    writeExtConductanceSettingsHDF5(FID, 1, F_ext, F*ones(1,N(1)) );
    writeExtConductanceSettingsHDF5(FID, 2, F*ones(1,N(2)), F*ones(1,N(2)) );
    
    % coupling range
    Drange=[ 45 45;
        45 45];
    
    sigma=[ 18 18 ;
        9e9 9e9];    
    
    % coupling strength
    alpha = 1.65;
    WI=0.035*alpha;
    WE=0.13*alpha;    
    
    W=[WE WE ;
        WI WI]; 
    
    SynapseType=[1 1 ;
        2 2];
    
    % use periodic boundary conditions    
    pbc=1; 
    
    % Write the above basic parameters to the input file
    writeBasicParaHDF5(FID, dt, step_tot, N);    
    
    %%%% Define the initial condition
    p_fire = 0*1e-2*ones(1,length(N)); % initial firing probabilities for populations
    % set initial V distribution to be [V_rt, V_rt + (V_th-V_rt)*r_V0]
    r_V0 = 1*ones(1,length(N));
    writeInitCondHDF5(FID, r_V0, p_fire)
    
    %%%% use initial conditions in files
    %     file_name = 2;
    %     outsetVE = load(sprintf('%04d_init_V_1.mat',file_name));
    %     outsetV_E = outsetVE.V;
    %     outsetVI = load(sprintf('%04d_init_V_2.mat',file_name));
    %     outsetV_I = outsetVI.V;
    %     writeExtInitVHDF5(FID,1,outsetV_E);
    %     writeExtInitVHDF5(FID,2,outsetV_I);
    
    %%%% Add stimulus
    
    writeExtCurrentPopHDF5(FID,ImageName_E,1);
    writeExtCurrentPopHDF5(FID,ImageName_I,2);
    
    
    %%%%%%%%%%%%%%%%%%% Chemical Connections %%%%%%%%%%%%%%%%%%%%%%%
    % type(1:AMAP, 2:GABAa, 3:NMDA)
    
    % Define connectivity between populations
    for pop_pre=1:length(N)
        [xpre,ypre] = find(Grid==pop_pre-1);
        for pop_post=1:length(N)
            if W(pop_pre,pop_post) ~= 0
                syn_type = SynapseType(pop_pre,pop_post); % 1 for AMPA-like synapse
                
                [xpost,ypost] = find(Grid==pop_post-1);
                if pbc
                    xdist=pdist2(xpre,xpost);
                    ydist=pdist2(ypre,ypost);
                    maxGridx=size(Grid,1);
                    maxGridy=size(Grid,1);
                    xdist=min(xdist,abs(maxGridx-xdist));
                    ydist=min(ydist,abs(maxGridy-ydist));
                    xdist2 = xdist.*xdist;
                    clear xdist
                    ydist2 = ydist.*ydist;
                    clear ydist
                    sqrsum = xdist2+ydist2;
                    clear xdist2 ydist2
                    DM=sqrt(sqrsum);
                    clear sqrsum                    
                else
                    DM=pdist2([xpre,ypre],[xpost,ypost]);  % FAST
                end
                if pop_pre == 2
                    A = (DM<=Drange(pop_pre,pop_post)).*W(pop_pre,pop_post);
                else
                    A=(DM<=Drange(pop_pre,pop_post)).*W(pop_pre,pop_post).*exp(-DM.^2/sigma(pop_pre,pop_post));
                end

                if pop_pre==pop_post
                    A(1:N(pop_pre)+1:N(pop_pre)*N(pop_pre)) = 0; %remove self conns
                end
                [I, J, K] = find(A);
                D = rand(size(I))*0; % uniformly random conduction delay ms
                writeChemicalConnectionHDF5(FID, syn_type,  pop_pre, pop_post, ...
                    I, J, K, D);
                clear I J K A D DM
            end
        end
    end

    
    %%%%%%%%%%%%%%%%%%% Chemical Connections Done %%%%%%%%%%%%%%%%%%%%%%%    
    
   
    %%%% synapse parameters conversion
    writeSynParaHDF5(FID, 'tau_decay_GABA', 3, 'Dt_trans_AMPA' , 0.3,...
        'tau_decay_AMPA', 2.0, 'Dt_trans_GABA', 0.3,'V_ex',V_ex,'V_in',V_in)
    
    % "help writeSynPara" to see all the parameter names and default values    
    
    %%% Sample detailed time series data
        for pop = 1:npops
            sample_neuron = 1:N(pop); % sample every ?th neuron
            sample_steps = zeros(1, step_tot); sample_steps(16000:10:24000) = 1; % sample every ? time step
            %sample_steps(1) = 1;
            sample_data_type = [1,0,0,0,0,0,0,0];
            % The logical vector above corresponds to
            % [V,I_leak,I_AMPA,I_GABA,I_NMDA,I_GJ,I_ext, I_K]
            writeNeuronSamplingHDF5(FID, pop, sample_data_type, ...
                sample_neuron, sample_steps)
        end
    %%%% Record explanatory variables that are nacessary for post-processing
    discard_transient = 0; % transient period data to be discarded (ms)
    writeExplVarHDF5(FID, 'discard_transient', discard_transient, ...
        'loop_num', loop_num, ...
        'F', F);        
    
  
    
    %%%% append this matlab file to the input file for future reference
    appendThisMatlabFileHDF5(FID)
    
end
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
