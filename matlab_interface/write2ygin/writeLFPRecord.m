function writeLFPRecord(FID, pop_ind, LFP_neurons)
%  write data recording for local field potential
%         FID: file id for writing data
%     pop_ind: 
% LFP_neurons: logical vector that specifies which neurons contribute the
%              LFP measure


% for C/C++ index convetion
pop_ind = pop_ind-1;

% write
fprintf(FID, '%s\n', '> SAMP005');
fprintf(FID, '%d,', pop_ind); fprintf(FID,'\n');
fprintf(FID, '%d,', LFP_neurons); fprintf(FID,'\n\n');
end