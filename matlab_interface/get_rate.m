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

% kernel for instantaneous rate estimation
CC_kernel_width = 50; % ms, kernel length
choice = 'gaussian_Hz';
kernel = spike_train_kernel_YG(CC_kernel_width, dt, choice);
    
rate = cell(Num_pop,1); % Initialise individual neuron average firing rate
Hz_t = cell(Num_pop,1); % population instantaneous rate estimation
Hz_overall = 0; % Network-wide average firing rate in Hz over simulation
for pop = 1:Num_pop
    if nnz(num_spikes{pop}) > 0
        rate{pop} = full(sum(spike_hist{pop},2))/(dt/1000*step_tot);
        Hz_t{pop} = SpikeTrainConvolve( sum(spike_hist{pop}, 1)/N(pop), kernel );
        Hz_overall = Hz_overall + sum(rate{pop})/sum(N);
    end
end
% Record results
R.Analysis.Hz_overall = Hz_overall;
R.Analysis.Hz_t = Hz_t;
R.Analysis.rate = rate;
end
