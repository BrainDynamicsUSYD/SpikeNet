function writeExtSpikeSettingsHDF5(FID, pop_ind, type_ext, K_ext,  Num_ext, rate_ext_t, neurons, seed)
% write external spike settings
%        FID: file id for writing data
%    pop_ind: index of neuron population to receive external spikes
%   type_ext: type of external chemical connection (1:AMPA, 2:GABA, 3:NMDA)
%      K_ext: strength for external chemical connection
%    Num_ext: number of external neurons connected to each neuron in pop_ind
% rate_ext_t: spiking rate for each time step
%    neurons: logical vector, true if the neuron receives the noisy input
%       seed: (optional) manually set RNG seed, a positive integer
%
% Note that each external neuron is independent Poissonian neuron

pop_ind = pop_ind - 1;
type_ext = type_ext - 1;
try
    n_syns = h5read(FID,'/config/syns/n_syns');
catch ME
    if (strcmp(ME.identifier,'MATLAB:imagesci:h5read:libraryError'))
        n_syns = 0;
        h5create(FID,'/config/syns/n_syns', 1);
        h5write(FID,'/config/syns/n_syns', n_syns);
    end
end
n_syns = n_syns + 1;
h5write(FID,'/config/syns/n_syns', n_syns);
n = n_syns - 1; % c++ zero-based index

% fprintf(FID, '%s\n', '# external spikes // (pop_ind, type_ext, K_ext:miuSiemens,  Num_ext, ia, ib;  rate_ext(t):Hz)');
hdf5write(FID,['/config/syns/syn',num2str(n),'/INIT005/pop_ind'],pop_ind,'WriteMode','append'); 
hdf5write(FID,['/config/syns/syn',num2str(n),'/INIT005/type_ext'],type_ext,'WriteMode','append'); 
hdf5write(FID,['/config/syns/syn',num2str(n),'/INIT005/K_ext'],K_ext,'WriteMode','append'); 
hdf5write(FID,['/config/syns/syn',num2str(n),'/INIT005/Num_ext'],Num_ext,'WriteMode','append'); 
hdf5write(FID,['/config/syns/syn',num2str(n),'/INIT005/neurons'],neurons,'WriteMode','append'); 
hdf5write(FID,['/config/syns/syn',num2str(n),'/INIT005/rate_ext_t'],rate_ext_t,'WriteMode','append'); 
fprintf('C++ external spikes: n_syn = %d, pop = %d, type = %d\n'  , n, pop_ind, type_ext);

% set seed
if nargin == 8
    if mod(seed,1) ~= 0 || seed < 0
        disp('The seed should be a positive integer!\n')
    else
        hdf5write(FID,['/config/syns/syn',num2str(n),'/SEED002/seed'],seed,'WriteMode','append');
        fprintf('\t RNG seed set manually.\n')
    end
end

end

