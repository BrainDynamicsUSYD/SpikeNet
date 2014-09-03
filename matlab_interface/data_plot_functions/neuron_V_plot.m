function neuron_V_plot(R, pop_ind, sample_ind)
% axes_matrix(1) = subplot(6, 8, 1:7 );hold on;


% Dump fields
dt = R.dt;
step_tot = R.step_tot;
T = (1:step_tot)*dt;
V = R.neuron_sample.V{pop_ind}(sample_ind,:);
V_th = R.PopPara{pop_ind}.V_th;

% plot potential
plot(T, V);

% show spkies
V_spike_gap = 5; % mV
spike_length = 10; % mV
X_a = T([(V(1:end-1) < V_th) & (V(2:end) >= V_th), false]);
X_b = X_a;
Y_a = ones(size(X_a))*V_th + V_spike_gap;
Y_b = Y_a + spike_length;
line([X_a;X_b], [Y_a;Y_b],'Color', 'k');

% % show threshold
% plot(T, V_th*ones(size(T)), 'r--');
ylim([  min( -70, min(V) ) V_th+10  ]); % [reset spike_peak]

% use scale bar instead axis label and ticks
scalebar( [10,10], {'10 ms','10 mV'}); 
        
end
