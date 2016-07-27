function mean_rate_plot(R, pop_ind, neuron_ind, seg, varargin)
 
sigma_gaussian = 50; % ms, which is width???

% Input check and default values
seg_size = 4*10^4; % 2*10^4 for 2-pop, segmentation size for each plot

if nargin < 4
    seg = 1;
end


text_fontsize = 12;
for i = 1:(length(varargin)/2)
    eval([varargin{i*2-1}, '=', num2str(varargin{i*2}) ]);
end



% Dump fields
dt = R.reduced.dt;
step_tot = R.reduced.step_tot;

% Segmetation
seg_ind = get_seg(step_tot, seg_size, seg);

% Dump fields
T = seg_ind*dt/1000;
spike_hist = R.reduced.spike_hist{pop_ind}(neuron_ind,seg_ind);

% Gaussian filter
kernel = spike_train_kernel_YG(sigma_gaussian, dt, 'gaussian');
kernel = spike_train_kernel_YG(1000, dt, 'square');

% mean rate
cluster_rate = SpikeTrainConvolve(sum(spike_hist, 1)/length(neuron_ind), kernel);
            
% rate plot
plot(T, cluster_rate);

ylabel('Hz','fontsize',text_fontsize); xlabel('t (sec)','fontsize',text_fontsize);
set(gca,'box','off','ticklength', [0 0],'fontsize',text_fontsize) %, 'TickDir','out');


end

