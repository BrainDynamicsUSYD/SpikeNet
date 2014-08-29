function writePopSampling(FID, pop_ind, sample_time_index)
% sample index should be logical vector with the length of time vector
% 4,000 neurons x 10,000 time steps will be 300MB data

if nnz(sample_time_index) > 10^4
    disp('Warning: too many time steps to be sampled!')
end
if sum((sample_time_index ~= 0) & (sample_time_index ~= 1)) ~= 0
    error('sample_time_index must be logical vectors: ');
end

pop_ind = pop_ind - 1;
% fprintf(FID, '%s\n', '# populational potential sampling // pop_ind; sample_time_index');
fprintf(FID, '%s\n', '> SAMP002');
fprintf(FID, '%d,', pop_ind); fprintf(FID,'\n');
fprintf(FID, '%d,', sample_time_index); fprintf(FID,'\n\n');
end



