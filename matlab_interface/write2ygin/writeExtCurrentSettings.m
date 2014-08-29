function writeExtCurrentSettings(FID, pop_ind, mean, std)
pop_ind = pop_ind - 1;
% fprintf(FID, '%s\n', '# external current // (pop_ind, mean, std) in nA');
fprintf(FID, '%s\n', '> INIT004');
fprintf(FID, '%.6f,', [pop_ind, mean, std]);fprintf(FID,'\n\n');
end

