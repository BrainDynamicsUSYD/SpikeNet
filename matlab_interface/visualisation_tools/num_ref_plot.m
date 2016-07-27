function num_ref_plot(R, pop_ind, seg, seg_size)

% Input check and default values
if nargin < 4
    seg_size = 4*10^4; % 2*10^4 for 2-pop, segmentation size for each plot
end
if nargin < 3
    seg = 1;
end

% Default line color
Color = [35 163 200]/255;

% Dump fields
dt = R.reduced.dt;
step_tot = R.reduced.step_tot;
N = R.N;

% Segmetation
seg_ind = get_seg(step_tot, seg_size, seg);

% Dump fields
T = seg_ind*dt;
num_ref = R.reduced.num_ref{pop_ind}(seg_ind);

% Plot number of refractory neurons
line([T; T], [zeros(1, length(T)); num_ref/N(pop_ind)*100], 'Color', Color);


ylabel('% Refractory');
xlabel('t (ms)');
set(gca, 'TickDir','out');

end

