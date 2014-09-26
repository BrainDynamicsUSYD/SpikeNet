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
dt = R.reduced.dt;
step_tot = R.reduced.step_tot;
N = R.N;
% rate_sorted = R.Analysis.rate_sorted;

% Segmetation
seg_num = ceil(step_tot/seg_size);
if seg < seg_num
    seg_ind = ((seg-1)*seg_size+1):(seg*seg_size);
else
    seg_ind = ((seg-1)*seg_size+1):(step_tot);
end

% Dump fields
num_spikes = R.reduced.num_spikes{pop_ind}(seg_ind);
spike_hist = R.reduced.spike_hist{pop_ind}(:,seg_ind);
T = seg_ind*dt;

% Plot raster plot
if nnz(num_spikes) > 0
    % down-sampling
    if N(pop_ind) >= sample_size
        ind_sample = sort(randperm(N(pop_ind),sample_size));
    else
        ind_sample = 1:1:N(pop_ind);
    end
    
    [Y,X,~] = find(spike_hist(ind_sample,:));
    line(([X(:)'; X(:)']+seg_ind(1)-1)*dt,[Y(:)'-1;Y(:)'],'Color','k');
    ylim([0,length(ind_sample)]);
    ylabel('Neurons')
%     if pop_ind == 1
%         if rate_sorted == 1
%             ylabel('Rate-sorted sample neuron index');
%         else
%             ylabel('Sample neuron index');
%         end
%     end
    
    xlim([T(1), T(1)+(length(seg_ind)-1)*dt]); % make sure all the plots have the same axis scale
    
    % Keep tick lables while remove tick marks
    % set(gca, 'xtick', [], 'Ticklength', [0 0], 'TickDir','out');
    axis off; box off;
    
        
end

end

