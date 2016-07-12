function writeExtCurrentSettings(FID, pop_ind, mean, std)
% write external current settings
%     FID: file id for writing data
% pop_ind: neuron population index
%    mean: mean value for Gaussian currents (nA) for each neuron
%     std: std for Gaussrian currents for each neuron

if length(mean) == 1 || length(std) == 1
    warning('INIT004 has been updated. MEAN and STD must be specified for each neuron.')
end

pop_ind = pop_ind - 1;

hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT004/mean'],mean,'WriteMode','append'); 
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT004/std'],std,'WriteMode','append'); 

end

