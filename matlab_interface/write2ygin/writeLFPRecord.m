function writeLFPRecord(FID, pop_ind, LFP_neurons)
%  write data recording for local field potential
%         FID: file id for writing data
%     pop_ind: 
% LFP_neurons: logical vector that specifies which neurons contribute the
%              LFP measure. For multiple LFP measures, use multiple rows.


% for C/C++ index convetion
pop_ind = pop_ind-1;

% write
[n_LFP, N] = size(LFP_neurons);
if n_LFP > N % if the size is not right
    warning('For multiple LFP measures, use multiple rows.')
else
    fprintf(FID, '%s\n', '> SAMP005');
    fprintf(FID, '%d,%d', pop_ind, n_LFP); fprintf(FID,'\n');
    for i = 1:n_LFP
        fprintf(FID, '%d,', LFP_neurons(i,:)); fprintf(FID,'\n');
    end
    fprintf(FID,'\n');
end

end