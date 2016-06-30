function writeExtConductanceSettings(FID, pop_ind, mean, std)
% write external Conductance settings
%     FID: file id for writing data
% pop_ind: neuron population index
%    mean: mean value for Gaussian conductance (uS) for each neuron
%     std: std for Gaussian conductance (uS) for each neuron

if length(mean) == 1 || length(std) == 1
    warning('INIT012: MEAN and STD must be specified for each neuron.')
end

pop_ind = pop_ind - 1;
% fprintf(FID, '%s\n', '# external current // (pop_ind, mean, std) in uS');
fprintf(FID, '%s\n', '> INIT012');
fprintf(FID, '%d,\n', pop_ind);
fprintf(FID, '%.9f,', mean);fprintf(FID,'\n');
fprintf(FID, '%.9f,', std);fprintf(FID,'\n\n');
end

