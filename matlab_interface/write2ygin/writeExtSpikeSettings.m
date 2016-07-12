function writeExtSpikeSettings(FID, pop_ind, type_ext, K_ext,  Num_ext, rate_ext_t, neurons)
% write external spike settings
%        FID: file id for writing data
%    pop_ind: index of neuron population to receive external spikes
%   type_ext: type of external chemical connection (1:AMPA, 2:GABA, 3:NMDA)
%      K_ext: strength for external chemical connection
%    Num_ext: number of external neurons connected to each neuron in pop_ind
% rate_ext_t: spiking rate for each time step
%    neurons: logical vector, true if the neuron receives the noisy input
% Note that each external neuron is independent Poissonian neuron

pop_ind = pop_ind - 1;
type_ext = type_ext - 1;



% fprintf(FID, '%s\n', '# external spikes // (pop_ind, type_ext, K_ext:miuSiemens,  Num_ext, ia, ib;  rate_ext(t):Hz)');
fprintf(FID, '%s\n', '> INIT005');
fprintf(FID, '%.9f,', [pop_ind, type_ext, K_ext,  Num_ext]); fprintf(FID,'\n');
fprintf(FID, '%d,', neurons); fprintf(FID,'\n');
fprintf(FID, '%.9f,', rate_ext_t); fprintf(FID,'\n\n');
end

