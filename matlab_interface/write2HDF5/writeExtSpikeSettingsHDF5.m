function writeExtSpikeSettingsHDF5(FID, pop_ind, type_ext, K_ext,  Num_ext, rate_ext_t, ia, ib)
% write external spike settings
%      FID: file id for writing data
%  pop_ind: index of neuron population to receive external spikes
% type_ext: type of external chemical connection (1:AMPA, 2:GABA, 3:NMDA)
%    K_ext: strength for external chemical connection
%  Num_ext: number of external neurons connected to each neuron in pop_ind
% rate_ext_t: spiking rate for each time step
%       ia: start neuron index
%       i:  end neuron index
% Note that each external neuron is independent Poissonian neuron

pop_ind = pop_ind - 1;
type_ext = type_ext - 1;


ia = ia - 1;
ib = ib - 1;


% fprintf(FID, '%s\n', '# external spikes // (pop_ind, type_ext, K_ext:miuSiemens,  Num_ext, ia, ib;  rate_ext(t):Hz)');
hdf5write(FID,['/config/pop',num2str(pop_ind),'/INIT005/type_ext'],type_ext,'WriteMode','append'); 
hdf5write(FID,['/config/pop',num2str(pop_ind),'/INIT005/K_ext'],K_ext,'WriteMode','append'); 
hdf5write(FID,['/config/pop',num2str(pop_ind),'/INIT005/Num_ext'],Num_ext,'WriteMode','append'); 
hdf5write(FID,['/config/pop',num2str(pop_ind),'/INIT005/ia'],ia,'WriteMode','append'); 
hdf5write(FID,['/config/pop',num2str(pop_ind),'/INIT005/ib'],ib,'WriteMode','append'); 
hdf5write(FID,['/config/pop',num2str(pop_ind),'/INIT005/rate_ext_t'],rate_ext_t,'WriteMode','append'); 

end

