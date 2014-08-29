function writeExtSpikeSettings(FID, pop_ind, type_ext, K_ext,  Num_ext, rate_ext)
pop_ind = pop_ind - 1;
type_ext = type_ext - 1;
% fprintf(FID, '%s\n', '# external spikes // (pop_ind, type_ext, K_ext:miuSiemens,  Num_ext;  rate_ext(t):Hz)');
fprintf(FID, '%s\n', '> INIT005');
fprintf(FID, '%.9f,', [pop_ind, type_ext, K_ext,  Num_ext]); fprintf(FID,'\n');
fprintf(FID, '%.9f,', rate_ext); fprintf(FID,'\n\n');
end

