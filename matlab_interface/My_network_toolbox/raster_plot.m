function raster_plot(R, pop_ind, seg, seg_size, sample_size)

% Input check and default values
if nargin < 5
    sample_size = 500; % sample neurons for raster plot
end
if nargin < 4
    seg_size = 4*10^4; % 2*10^4 for 2-pop, segmentation size for each plot
end
if nargin < 3
    seg = 1;
end

% Dump fields
reduced_dt = R.reduced_dt;
reduced_step_tot = R.reduced_step_tot;
N = R.N;
rate_sorted = R.Analysis.rate_sorted;

% Segmetation
seg_num = ceil(reduced_step_tot/seg_size);
if seg < seg_num
    seg_ind = ((seg-1)*seg_size+1):(seg*seg_size);
else
    seg_ind = ((seg-1)*seg_size+1):(reduced_step_tot);
end

% Dump fields
reduced_num_spikes = R.reduced_num_spikes{pop_ind}(seg_ind);
reduced_spike_hist = R.reduced_spike_hist{pop_ind}(:,seg_ind);
reduced_T = seg_ind*reduced_dt;

% Plot raster plot
if nnz(reduced_num_spikes) > 0
    % down-sampling
    if N(pop_ind) >= sample_size
        ind_sample = sort(randperm(N(pop_ind),sample_size));
    else
        ind_sample = 1:1:N(pop_ind);
    end
    
    [Y,X,~] = find(reduced_spike_hist(ind_sample,:));
    line(([X(:)'; X(:)']+seg_ind(1)-1)*reduced_dt,[Y(:)'-1;Y(:)'],'Color','k');
    ylim([0,length(ind_sample)]);
    if pop_ind == 1
        if rate_sorted == 1
            ylabel('Rate-sorted sample neuron index');
        else
            ylabel('Sample neuron index');
        end
    end
    set(gca, 'xtick', []);
    xlim([reduced_T(1), reduced_T(1)+(length(seg_ind)-1)*reduced_dt]); % make sure all the plots have the same axis scale
end

end

