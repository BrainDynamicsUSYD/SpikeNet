function writeSynStatsRecord(FID, pop_ind_pre, pop_ind_post, syn_type)



% for C/C++ index convetion
pop_ind_pre = pop_ind_pre-1;
pop_ind_post = pop_ind_post-1;
syn_type = syn_type-1;

% write
fprintf(FID, '%s\n', '> SAMP004');
fprintf(FID, '%d, %d, %d,\n', pop_ind_pre, pop_ind_post, syn_type);

end