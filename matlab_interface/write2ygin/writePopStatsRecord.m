function writePopStatsRecord(FID, pop_ind)

% for C/C++ index convetion
pop_ind = pop_ind-1;

% write
fprintf(FID, '%s\n', '> SAMP003');
fprintf(FID, '%d,\n', pop_ind);

end