function writeExtCurrentTimeVariantSettingsMultiGroupHDF5(FID, pop_ind, mean, std, mean_TV, std_TV, TV_group)
% write external current settings with time-variant factors
%       FID: file id for writing data
%   pop_ind: neuron population index
%      mean: mean value for Gaussian currents (nA) for each neuron
%       std: std for Gaussrian currents for each neuron
%   mean_TV: a time-variant factor multiplied to the mean value for 
%            Gaussian currents (nA) for each neuron, 
%            size: number_of_groups-by-number_of_time_steps
%    std_TV: a time-variant factor multiplied to the std value for
%            Gaussrian currents for each neuron
%            size: number_of_groups-by-number_of_time_steps
%  TV_group: the group number of each neuron (for choosing the row in mean_TV and std_TV)
%
% Explanation: For neuron i, at time step t, the external current will be a 
% Gaussian random number with a mean of mean[i]*mean_TV[TV_group[i],t] and a std of std[i]*std_TV[TV_group[i],t]	

pop_ind = pop_ind - 1;
if min(TV_group) == 0
    error('Group number must start from 1 (matlab indexing).')
end
TV_group = TV_group - 1;

[n_group, step_tot] = size(mean_TV);
if n_group > step_tot % if the size is not right
    error('For multiple groups, use multiple rows.')
else
    mean_TV = mean_TV'; % columns
    std_TV = std_TV'; % columns
    hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT0019/mean_TV'],mean_TV,'WriteMode','append');
    hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT0019/std_TV'],std_TV,'WriteMode','append');
    hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT0019/mean'],mean,'WriteMode','append');
    hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT0019/std'],std,'WriteMode','append');
    hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT0019/TV_group'],TV_group,'WriteMode','append');

end
