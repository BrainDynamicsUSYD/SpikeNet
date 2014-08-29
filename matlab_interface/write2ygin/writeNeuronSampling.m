function writeNeuronSampling(FID, pop_ind, data_type, sample_ind)
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
