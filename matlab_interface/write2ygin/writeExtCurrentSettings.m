function writeExtCurrentSettings(FID, pop_ind, mean, std)
% write external current settings
%     FID: file id for writing data
% pop_ind: neuron population index
%    mean: mean value for Gaussian currents (nA) for each neuron
%     std: std for Gaussrian currents for each neuron

if length(mean) == 1 || length(std) == 1
    warning('INIT004 has been updated. MEAN and STD must be specified for each neuron.')
end

pop_ind = pop_ind - 1;
% fprintf(FID, '%s\n', '# external current // (pop_ind, mean, std) in nA');
fprintf(FID, '%s\n', '> INIT004');
fprintf(FID, '%d,\n', pop_ind);
fprintf(FID, '%.9f,', mean);fprintf(FID,'\n');
fprintf(FID, '%.9f,', std);fprintf(FID,'\n\n');
end

