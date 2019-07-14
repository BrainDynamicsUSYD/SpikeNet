function writeExtCurrentTimeVariantSettingsHDF5(FID, pop_ind, mean, std, mean_TV, std_TV, modify)
% write external current settings with time-variant factors
%       FID: file id for writing data
%   pop_ind: neuron population index
%      mean: mean value for Gaussian currents (nA) for each neuron
%       std: std for Gaussrian currents for each neuron
%   mean_TV: a time-variant factor multiplied to the mean value for Gaussian currents (nA) for each neuron
%    std_TV: a time-variant factor multiplied to the std value for Gaussrian currents for each neuron
%    modify: whether this is for modifying existing HDF5 file content or not
%
% Explanation: For neuron i, at time step t, the external current will be a Gaussian random number with a mean of mean[i]*mean_TV[t] and a std of std[i]*std_TV[t]	

if length(mean) == 1 || length(std) == 1
    warning('INIT004 has been updated. MEAN and STD must be specified for each neuron.')
end

pop_ind = pop_ind - 1;

if nargin == 6
    modify = 0;
end

if modify == 0
    h5create(FID,['/config/pops/pop',num2str(pop_ind),'/INIT014/mean'],length(mean));
    h5create(FID,['/config/pops/pop',num2str(pop_ind),'/INIT014/std'],length(std));
    h5create(FID,['/config/pops/pop',num2str(pop_ind),'/INIT014/mean_TV'],length(mean_TV));
    h5create(FID,['/config/pops/pop',num2str(pop_ind),'/INIT014/std_TV'],length(std_TV));
end

h5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT014/mean'],mean);
h5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT014/std'],std);
h5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT014/mean_TV'],mean_TV);
h5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT014/std_TV'],std_TV);

end
