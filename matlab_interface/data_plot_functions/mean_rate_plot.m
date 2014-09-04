function mean_rate_plot(R, pop_ind, neuron_ind, seg, seg_size)
 
sigma_gaussian = 50; % ms, which is width???

% Input check and default values
if nargin < 5
    seg_size = 4*10^4; % 2*10^4 for 2-pop, segmentation size for each plot
end
if nargin < 4
    seg = 1;
end


% Dump fields
dt = R.reduced_dt;
step_tot = R.reduced_step_tot;

% Segmetation
seg_num = ceil(step_tot/seg_size);
if seg < seg_num
    seg_ind = ((seg-1)*seg_size+1):(seg*seg_size);
else
    seg_ind = ((seg-1)*seg_size+1):(step_tot);
end


% Dump fields
T = seg_ind*dt;
spike_hist = R.reduced_spike_hist{pop_ind}(neuron_ind,seg_ind);

% Gaussian filter
kernel = spike_train_kernel_YG(sigma_gaussian, dt, 'gaussian');

% mean rate
cluster_rate = SpikeTrainConvolve(sum(spike_hist, 1)/length(neuron_ind), kernel);
            
% rate plot
plot(T, cluster_rate);

ylabel('Hz'); xlabel('t (ms)');
set(gca,'box','off', 'TickDir','out');


end

