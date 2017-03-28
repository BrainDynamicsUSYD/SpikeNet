function writeSynSamplingHDF5(FID, pop_ind_pre, pop_ind_post, syn_type, sample_ind, time_index)
% write data sampling settings for individual neurons
%        FID: file id for writing data
% sample_ind: vector of post-synaptic neuron indices to be sampled in pop_ind
% time_index: logical vector to define time points to be sampled
% 
% Note that time_index should have the length of total simulation steps.
% For example, if step_tot = 10 and the last half time points are to be
% sampled, use time_index = [0,0,0,0,0,1,1,1,1,1]
% 4,000 neurons x 10,000 time points will be 300MB data


if sum((time_index ~= 0) & (time_index ~= 1)) ~= 0
    error('sample_time_index must be logical vectors: ');
end

% for C/C++ index convetion
pop_ind_pre = pop_ind_pre-1;
pop_ind_post = pop_ind_post-1;
syn_type = syn_type-1;
sample_ind = sample_ind-1;

% sampling frequency check
dt = hdf5read(FID,'/config/Net/INIT002/dt');
freq_s = 1000 / (mode(diff(find(time_index)))*dt);
if freq_s < 1000
    warning('You are assuming the Nyquist frequency of the system is lower than %g Hz. /n', freq_s);
end

n_syns = h5read(FID,'/config/syns/n_syns');
n_match = NaN; 
for n = (1:n_syns)-1
    try
        type = hdf5read(FID,['/config/syns/syn',num2str(n),'/INIT006/type']);
        i_pre = hdf5read(FID,['/config/syns/syn',num2str(n),'/INIT006/i_pre']);
        j_post = hdf5read(FID,['/config/syns/syn',num2str(n),'/INIT006/j_post']);
        
        if pop_ind_pre == i_pre && pop_ind_post == j_post &&  syn_type == type
            n_match = n;
        end
    catch ME
    end
end

if isnan(n_match)
    error('Cannot find syn with identical pop_ind_pre, pop_ind_post and syn_type!')
else
    hdf5write(FID,['/config/syns/syn',num2str(n_match),'/SAMP002/neurons'], sample_ind,'WriteMode','append');
    hdf5write(FID,['/config/syns/syn',num2str(n_match),'/SAMP002/time_points'], time_index,'WriteMode','append');
end

end

