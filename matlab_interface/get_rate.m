function [R] = get_rate(R)
fprintf('\t Getting firing rate...\n');
% dump fields
dt = R.dt;
step_tot = R.step_tot;
N = R.N;
Num_pop = R.Num_pop;
num_spikes = R.num_spikes;
spike_hist = R.spike_hist;
%
rate = cell(Num_pop,1); % Initialise individual neuron average firing rate
Hz_overall = 0; % Network-wide average firing rate in Hz over simulation
for pop_ind = 1:Num_pop
    if nnz(num_spikes{pop_ind}) > 0
        rate{pop_ind} = full(sum(spike_hist{pop_ind},2))/(dt/1000*step_tot);
        Hz_overall = Hz_overall + sum(rate{pop_ind})/sum(N);
    end
end
% Record results
R.Analysis.Hz_overall = Hz_overall;
R.Analysis.rate = rate;
end
