function writeRunawayKiller(FID, runaway_steps, mean_num_ref)
% write runaway activity killer
%            FID: file id for writing data
%  runaway_steps: see mean_num_ref
%   mean_num_ref: mean number of refractory neurons in a population at one
%                 time step averaged over the lastest runaway_steps.
%
% So if any neuron population's activity exceeds the threshold define by
% (runaway_steps,mean_num_ref), the simulation will be killed.
% Note that only those populations with (by defualt) >100 neurons will be
% considered.

if ~(mean_num_ref >= 0 && mean_num_ref <= 1)
    disp('mean_num_ref must be between [0,1]!');
end
%fprintf(FID, '%s\n', '# runaway killer setting //runaway_steps, runaway_mean_num_ref(0~1),');
fprintf(FID, '%s\n', '> KILL001');
fprintf(FID, '%d,%f,\n\n', runaway_steps, mean_num_ref);
end