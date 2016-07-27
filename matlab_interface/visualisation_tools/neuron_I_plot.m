function neuron_I_plot(R, pop_ind, sample_ind, current_type, seg, seg_size)
% current_type can be any non-empty subset of 
% {'I_leak','I_AMPA','I_GABA','I_NMDA','I_GJ','I_ext'}

marker = {'r','g','b','y','k','c','m'}; % for 6 different currents


% Input check and default values
if nargin < 6
    seg_size = 4*10^4; % 2*10^4 for 2-pop, segmentation size for each plot
end
if nargin < 5
    seg = 1;
end


% Dump fields
dt = R.dt;
step_tot = R.step_tot;

% Segmetation
seg_ind = get_seg(step_tot, seg_size, seg);

T = seg_ind*dt;

% plot
hold on;
for type = 1:length(current_type)
        I = R.neuron_sample.(current_type{type}){pop_ind}(sample_ind,seg_ind);
        plot(T,I,marker{type});
end
lh = legend(current_type);
set(lh,'Interpreter','none'); % plain text instead of Tex

set(gca,'box','off', 'TickDir','out')

end