function writePopSampling(FID, pop_ind, data_type, time_index)
% write data sampling settings for entire population
%         FID: file id for writing data
%     pop_ind: neuron population index
%   data_type: logical vector for [V,I_leak,I_AMPA,I_GABA,I_NMDA,I_GJ,I_ext]
%  time_index: logical vector to define time points to be sampled
%
% For example, if V and I_GABA are needed, use data_type = [1,0,0,1,0,0,0]
% 
% Note that time_index should have the length of total simulation steps.
% For example, if step_tot = 10 and the last half time points are to be
% sampled, use time_index = [0,0,0,0,0,1,1,1,1,1]
% 4,000 neurons x 10,000 time points will be 300MB data

if nnz(time_index) > 10^4
    warning('Too many time steps to be sampled!');
end
if sum((time_index ~= 0) & (time_index ~= 1)) ~= 0
    error('sample_time_index must be logical vectors: ');
end

pop_ind = pop_ind - 1;
% fprintf(FID, '%s\n', '# populational potential sampling // pop_ind; sample_time_index');
fprintf(FID, '%s\n', '> SAMP002');
fprintf(FID, '%d,', pop_ind); fprintf(FID,'\n');
fprintf(FID, '%d,', data_type); fprintf(FID,'\n');
fprintf(FID, '%d,', time_index); fprintf(FID,'\n\n');
end
