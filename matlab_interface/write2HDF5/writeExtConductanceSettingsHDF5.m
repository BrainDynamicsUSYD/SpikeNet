function writeExtConductanceSettingsHDF5(FID, pop_ind, mean, std, modify)
% write external Conductance settings
%     FID: file id for writing data
% pop_ind: neuron population index
%    mean: mean value for Gaussian conductance (uS) for each neuron
%     std: std for Gaussian conductance (uS) for each neuron

if length(mean) == 1 || length(std) == 1
    warning('INIT012: MEAN and STD must be specified for each neuron.')
end

pop_ind = pop_ind - 1;

if nargin == 4
     modify = 0;
end

if  modify == 0
    h5create(FID,['/config/pops/pop',num2str(pop_ind),'/INIT012/mean'],length(mean));
    h5create(FID,['/config/pops/pop',num2str(pop_ind),'/INIT012/std'],length(std));
end
h5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT012/mean'],mean);
h5write(FID,['/config/pops/pop',num2str(pop_ind),'/INIT012/std'],std);
end

