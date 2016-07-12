function writeExtConductanceSettingsHDF5(FID, pop_ind, mean, std)
% write external Conductance settings
%     FID: file id for writing data
% pop_ind: neuron population index
%    mean: mean value for Gaussian conductance (uS) for each neuron
%     std: std for Gaussian conductance (uS) for each neuron

if length(mean) == 1 || length(std) == 1
    warning('INIT012: MEAN and STD must be specified for each neuron.')
end

pop_ind = pop_ind - 1;

hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT012/mean'],mean,'WriteMode','append');
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT012/std'],std,'WriteMode','append');
end

