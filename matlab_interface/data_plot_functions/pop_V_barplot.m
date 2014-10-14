function pop_V_barplot(R, pop_ind, sample_t_ind, neuron_ind)
% axes_matrix(1) = subplot(6, 8, 1:7 );hold on;

bin_num = 100;

% Dump fields
if nargin < 4
    neuron_ind = 1:R.N(pop_ind);
end
V = R.pop_sample.V{pop_ind}(:,sample_t_ind);
% t_ind = R.pop_sample.t_ind{pop_ind}(sample_t_ind);
% num_ref = R.num_ref{pop_ind}(t_ind-1:t_ind+1);
% V_rt = R.PopPara{pop_ind}.V_rt;
% V = Remove_V_rt(V, V_rt, num_ref);% remove V_rt (reset potential duing refractory period)

V = V(neuron_ind); % extract V for neuron_ind
% V( isnan(V) ) = [];

% plot
bin_edge = linspace(min(V), max(V), bin_num+1);
bin_count = histc(V,bin_edge);
bin_count = bin_count/sum(bin_count); % normalization
barh(bin_edge, bin_count, 'histc');


set(gca,'xtick',[],'box','off', 'TickDir','out', 'XColor','w');
set(gcf, 'InvertHardCopy', 'off'); % prevent the x-axis line to reappear when printed
ylabel('mV');
ylim([-70,-50]);
axis off;

end


% function V = Remove_V_rt(V, V_rt, num_ref)
% ref = find(V == V_rt);
% ref = ref(1:num_ref);
% V(ref) = NaN;
% end