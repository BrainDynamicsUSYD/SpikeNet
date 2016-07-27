function num_spikes_plot(R, pop_ind, seg, seg_size)

% Input check and default values
if nargin < 4
    seg_size = 4*10^4; % 2*10^4 for 2-pop, segmentation size for each plot
end
if nargin < 3
    seg = 1;
end

% Default line color
Color = [255 30 30]/255;

% Dump fields
dt = R.reduced.dt;
step_tot = R.reduced.step_tot;
N = R.N;

%
dt = dt/1000; % sec

% Segmetation
seg_ind = get_seg(step_tot, seg_size, seg);

% Dump fields
num_spikes = R.reduced.num_spikes{pop_ind}(seg_ind);
T = seg_ind*dt;

%  Plot number of spikes
line([T; T], [zeros(1, length(T)); num_spikes/N(pop_ind)*100], 'Color', Color);


ylabel('% Firing');
set(gca, 'TickDir','out');
xlabel('t (sec)');

end