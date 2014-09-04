function neuron_rate_plot(R, pop_ind, sample_ind, seg, seg_size)
 
sigma_gaussian = 50; % ms, which is width???

% Input check and default values
if nargin < 5
    seg_size = 4*10^4; % 2*10^4 for 2-pop, segmentation size for each plot
end
if nargin < 4
    seg = 1;
end


% Dump fields
dt = R.dt;
step_tot = R.step_tot;

% Segmetation
seg_num = ceil(step_tot/seg_size);
if seg < seg_num
    seg_ind = ((seg-1)*seg_size+1):(seg*seg_size);
else
    seg_ind = ((seg-1)*seg_size+1):(step_tot);
end


% Dump fields
T = seg_ind*dt;
ind = R.neuron_sample.neuron_ind{pop_ind}(sample_ind); % need to record neuron_ind!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
spike_hist = R.spike_hist{pop_ind}(ind,  seg_ind);

% Gaussian filter
kernel = spike_train_kernel_YG(sigma_gaussian, dt, 'gaussian');
            
% rate plot
rate = SpikeTrainConvolve(spike_hist, kernel);
plot(T, rate);
ylabel('Hz'); xlabel('t (ms)');
set(gca,'box','off','TickDir','out');


end
