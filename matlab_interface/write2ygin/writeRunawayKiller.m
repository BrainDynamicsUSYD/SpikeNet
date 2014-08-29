function writeRunawayKiller(FID, steps, mean_num_ref)
if ~(mean_num_ref >= 0 && mean_num_ref <= 1)
    disp('mean_num_ref must be between [0,1]!');
end
%fprintf(FID, '%s\n', '# runaway killer setting //runaway_steps, runaway_mean_num_ref(0~1),');
fprintf(FID, '%s\n', '> KILL001');
fprintf(FID, '%d,%f,\n\n', steps, mean_num_ref);
end