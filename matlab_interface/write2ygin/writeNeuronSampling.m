function writeNeuronSampling(FID, pop_ind, data_type, sample_ind, time_index)
% write data sampling settings for individual neurons
%        FID: file id for writing data
%    pop_ind: neuron population index
%  data_type: logical vector for [V,I_leak,I_AMPA,I_GABA,I_NMDA,I_GJ,I_ext, I_K]
% sample_ind: vector of neuron indices to be sampled in pop_ind
% time_index: logical vector to define time points to be sampled
%
% For example, if V and I_GABA are needed, use data_type = [1,0,0,1,0,0,0,0]
% 
% Note that time_index should have the length of total simulation steps.
% For example, if step_tot = 10 and the last half time points are to be
% sampled, use time_index = [0,0,0,0,0,1,1,1,1,1]
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

if length(data_type) ~= 8
    error('data_type must have a length of 8.')
end

% write
% fprintf(FID, '%s\n', '# neuronal membrane potential and currents sampling setting // pop_ind;sample_ind');
fprintf(FID, '%s\n', '> SAMP001');
fprintf(FID, '%d,\n', pop_ind);
fprintf(FID, '%d,', data_type); fprintf(FID,'\n');
fprintf(FID, '%d,', sample_ind); fprintf(FID,'\n');
fprintf(FID, '%d,', time_index); fprintf(FID,'\n\n');