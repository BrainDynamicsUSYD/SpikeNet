
function neuron_V_barplot(R, pop_ind, sample_ind, YLim, varargin)
% subplot(6,8,8);hold on;  

bin_num = 100;

if nargin == 4
    set(gca,'YLim',YLim);
end

text_fontsize = 12;
for i = 1:(length(varargin)/2)
    eval([varargin{i*2-1}, '=', num2str(varargin{i*2}) ]);
end

% Dump fields
dt = R.dt;
V = R.neuron_sample.V{pop_ind}(sample_ind,:);
V_th = R.PopPara{pop_ind}.V_th;
tau_ref = R.PopPara{pop_ind}.tau_ref;
ref_steps = round(tau_ref/dt);

% plot
V_remove = RemoveSpikeRef(V, V_th, ref_steps); % remove spike and refractory data
bin_edge = linspace(min(V_remove), max(V_remove), bin_num+1);
bin_count = histc(V_remove,bin_edge);
bin_count = bin_count/sum(bin_count); % normalization
barh(bin_edge, bin_count, 'histc');


set(gca,'xtick',[],'box','off', 'TickDir','out', 'XColor','w');
set(gcf, 'InvertHardCopy', 'off'); % prevent the x-axis line to reappear when printed

ylim([-70,-50]);
ylabel('mV','fontsize', text_fontsize);

end

function V = RemoveSpikeRef(V, theta, ref_steps)
% logical vector of spiking event
spike = [(V(1:end-1) < theta) & (V(2:end) >= theta), false];
% logical vector of being refractory
ref = spike;
for i = 1:ref_steps
    ref = ref | [false(1,i) spike(1:end-i)];
end
% delete refractoriness from V
V(ref) = [];
end
