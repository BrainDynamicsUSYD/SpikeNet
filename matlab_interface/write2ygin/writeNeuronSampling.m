function writeNeuronSampling(FID, pop_ind, data_type, sample_ind)
% write data sampling settings for individual neurons
%        FID: file id for writing data
%    pop_ind: neuron population index
%  data_type: logical vector for [V,I_leak,I_AMPA,I_GABA,I_NMDA,I_GJ,I_ext]
% sample_ind: vector of neuron indices to be sampled in pop_ind
%
% For example, if V and I_GABA are needed, use data_type = [1,0,0,1,0,0,0]
% Entire time sequences of the data will be recorded for individual neurons.


% for C/C++ index convetion
pop_ind = pop_ind-1;
sample_ind = sample_ind-1;
if sum((data_type ~= 0) & (data_type ~= 1)) ~= 0
    error('Data_type must be logical vectors: ');
end
% write
% fprintf(FID, '%s\n', '# neuronal membrane potential and currents sampling setting // pop_ind;sample_ind');
fprintf(FID, '%s\n', '> SAMP001');
fprintf(FID, '%d,\n', pop_ind);
fprintf(FID, '%d,', data_type); fprintf(FID,'\n');
fprintf(FID, '%d,', sample_ind); fprintf(FID,'\n');
fprintf(FID, '\n');
