function writePoissonPopHDF5(FID, pop_ind, rate)
% This function sets the population pop_ind to be a population 
% where each neuron generates a Poisson spike train 
% rate is assumed to be in units of spikes per second


pop_ind = pop_ind - 1; % from matlab to c++ index
hdf5write(FID,['/config/pops/pop',num2str(pop_ind),'/poisson_pop/rate'],rate,'WriteMode','append'); 

end

