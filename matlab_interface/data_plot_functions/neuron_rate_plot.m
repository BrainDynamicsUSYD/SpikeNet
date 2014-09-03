function neuron_rate_plot(R, pop_ind, sample_ind)

sigma_gaussian = 50; % ms, which is width???

% Dump fields
dt = R.dt;
step_tot = R.step_tot;
T = (1:step_tot)*dt;
ind = R.neuron_sample.neuron_ind{pop_ind}(sample_ind); % need to record neuron_ind!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
spike_hist = R.spike_hist{pop_ind}(ind, :);

% Gaussian filter
kernel = spike_train_kernel_YG(sigma_gaussian, dt, 'gaussian');
            
% rate plot
rate = SpikeTrainConvolve(spike_hist, kernel);
plot(T, rate);
set(gca,'xtick',[],'box','off');


end

