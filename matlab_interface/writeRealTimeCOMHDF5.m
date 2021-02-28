function writeRealTimeCOMHDF5(FID, pop_ind, data_type, time_index)
% write real-time center-of-mass recording settings for specific time points
%        FID: file id for writing data
%    pop_ind: neuron population index
%  data_type: logical vector for [V,I]
% time_index: logical vector to define time points to be sampled
%
% For example, if V and I are needed, use data_type = [1,1]
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
if sum((data_type ~= 0) & (data_type ~= 1)) ~= 0
    error('Data_type must be logical vectors: ');
end

if length(data_type) < 2
    error('data_type must have a length of at least 2.')
end

% sampling frequency check
dt = hdf5read(FID,'/config/Net/INIT002/dt');
freq_s = 1000 / (mode(diff(find(time_index)))*dt);
if freq_s < 1000
    warning('You are assuming the Nyquist frequency of the system is lower than %g Hz. /n', freq_s);
end

% write
% fprintf(FID, '%s\n', '# neuronal membrane potential and currents sampling setting // pop_ind;sample_ind');
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SAMP006/data_type/V_flag'],data_type(1),'WriteMode','append'); 
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SAMP006/data_type/I_flag'],data_type(2),'WriteMode','append'); 
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/SAMP006/time_points'],time_index,'WriteMode','append');
end
