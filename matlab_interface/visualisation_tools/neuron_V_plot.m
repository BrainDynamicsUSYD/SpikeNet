function neuron_V_plot(R, pop_ind, sample_ind, seg, varargin)
% axes_matrix(1) = subplot(6, 8, 1:7 );hold on;


% Dump fields
dt = R.dt;
step_tot = R.step_tot;

% Input check and default values
seg_size = 4*10^4; % 2*10^4 for 2-pop, segmentation size for each plot
if nargin < 4
    seg = 1;
end


text_fontsize = 12;
for i = 1:(length(varargin)/2)
    eval([varargin{i*2-1}, '=', num2str(varargin{i*2}) ]);
end

% Segmetation
seg_ind = get_seg(step_tot, seg_size, seg);


T = seg_ind*dt;
V = R.neuron_sample.V{pop_ind}(sample_ind,seg_ind);
V_th = R.PopPara{pop_ind}.V_th;

% plot potential
plot(T, V);

% show spkies
V_spike_gap = 2; % mV
spike_length = 3; % mV
X_a = T([(V(1:end-1) < V_th) & (V(2:end) >= V_th), false]);
X_b = X_a;
Y_a = ones(size(X_a))*V_th + V_spike_gap;
Y_b = Y_a + spike_length;
line([X_a;X_b], [Y_a;Y_b],'Color', 'k');

% % show threshold
% plot(T, V_th*ones(size(T)), 'r--');
ylim([  min( -70, min(V) ) V_th+10  ]); % [reset spike_peak]

% use scale bar instead axis label and ticks
% scalebar( [10,10], {'10 ms','10 mV'}); 
% scalebar( [1000,10], {'1 sec','10 mV'}); 
scalebar( [0,10], { '','10 mV'}, 'text_fontsize',text_fontsize); 
        
end
