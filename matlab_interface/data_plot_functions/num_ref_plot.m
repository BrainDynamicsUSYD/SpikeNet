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
reduced_dt = R.reduced_dt;
reduced_step_tot = R.reduced_step_tot;
N = R.N;

% Segmetation
seg_num = ceil(reduced_step_tot/seg_size);
if seg < seg_num
    seg_ind = ((seg-1)*seg_size+1):(seg*seg_size);
else
    seg_ind = ((seg-1)*seg_size+1):(reduced_step_tot);
end

% Dump fields
reduced_T = seg_ind*reduced_dt;
reduced_num_ref = R.reduced_num_ref{pop_ind}(seg_ind);

% Plot number of refractory neurons
line([reduced_T; reduced_T], [zeros(1, length(reduced_T)); reduced_num_ref/N(pop_ind)*100], 'Color', Color);
if pop_ind == 1
    ylabel('% Refractory');
end
xlabel(['t (ms), pop No.', num2str(pop_ind)]);

end

