function writeSpikeFreqAdpt(FID, pop_ind)
% write spike-frequency adaptation setting
% ref: Alessandro Treves, 1993, Mean-field analysis of neuronal spike dynamics
%          FID: file id for writing data
%      pop_ind:

fprintf(FID, '%s\n', '> INIT010');
fprintf(FID, '%d,', pop_ind-1);
fprintf(FID,'\n\n');
end

