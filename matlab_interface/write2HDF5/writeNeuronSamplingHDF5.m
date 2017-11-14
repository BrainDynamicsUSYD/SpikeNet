function writeNeuronSamplingHDF5(FID, pop_ind, data_type, sample_ind, time_index)
% write data sampling settings for individual neurons
%        FID: file id for writing data
%    pop_ind: neuron population index
%  data_type: logical vector for [V,I_leak,I_AMPA,I_GABA,I_NMDA,I_GJ,I_ext, I_K,rhat]
%		note rhat requires JH_Learning turned on and is optional
% sample_ind: vector of neuron indices to be sampled in pop_ind
% time_index: logical vector to define time points to be sampled
%
% For example, if V and I_GABA are needed, use data_type = [1,0,0,1,0,0,0,0,0]
% 
% Note that time_index should have the length of total simulation steps.
% For example, if step_tot = 10 and the last half time points are to be
% sampled, use time_index = [0,0,0,0,0,1,1,1,1,1,0]
% 4,000 neurons x 10,000 time points will be 300MB data


% if sum((time_index ~= 0) & (time_index ~= 1)) ~= 0
%     error('sample_time_index must be logical vectors: ');
% end

% for C/C++ index convetion
pop_ind = pop_ind-1;
sample_ind = sample_ind-1;
if sum((data_type ~= 0) & (data_type ~= 1)) ~= 0
    error('Data_type must be logical vectors: ');
end

if length(data_type) < 8
    error('data_type must have a length of at least 8.')
end

% sampling frequency check
dt = hdf5read(FID,'/config/Net/INIT002/dt');
freq_s = 1000 / (mode(diff(find(time_index)))*dt);
if freq_s < 1000
    warning('You are assuming the Nyquist frequency of the system is lower than %g Hz. /n', freq_s);
end

% write
% fprintf(FID, '%s\n', '# neuronal membrane potential and currents sampling setting // pop_ind;sample_ind');
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SAMP001/data_type/V'],data_type(1),'WriteMode','append'); 
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SAMP001/data_type/I_leak'],data_type(2),'WriteMode','append'); 
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SAMP001/data_type/I_AMPA'],data_type(3),'WriteMode','append'); 
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SAMP001/data_type/I_GABA'],data_type(4),'WriteMode','append'); 
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SAMP001/data_type/I_NMDA'],data_type(5),'WriteMode','append'); 
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SAMP001/data_type/I_GJ'],data_type(6),'WriteMode','append'); 
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SAMP001/data_type/I_ext'],data_type(7),'WriteMode','append'); 
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SAMP001/data_type/I_K'],data_type(8),'WriteMode','append'); 
if length(data_type)==9
	hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SAMP001/data_type/rhat'],data_type(9),'WriteMode','append')
end


hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SAMP001/neurons'],sample_ind,'WriteMode','append');
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SAMP001/time_points'],time_index,'WriteMode','append');
end
