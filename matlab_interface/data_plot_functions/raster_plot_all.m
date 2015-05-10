function raster_plot_all(R, pop_sample, seg, varargin)
% pop_sample: [400 100] means sample 400 neurons from pop1 and 100 from
% pop2 and plot them in one raster plot

% Input check and default values
if nargin < 3
    seg = 1;
end

seg_size = 4*10^4; % 2*10^4 for 2-pop, segmentation size for each plot

text_fontsize = 12;
for i = 1:(length(varargin)/2)
    eval([varargin{i*2-1}, '=', num2str(varargin{i*2}), ';' ]);
end



% Dump fields
dt = R.reduced.dt;
step_tot = R.reduced.step_tot;
N = R.N;
% rate_sorted = R.Analysis.rate_sorted;

% 
dt = dt/1000; %second

% Segmetation
seg_ind = get_seg(step_tot, seg_size, seg);

% Dump fields
T = seg_ind*dt;

neurons_sampled = 0;    
for pop_ind = 1:length(N)
    spike_hist = R.reduced.spike_hist{pop_ind}(:,seg_ind);
    % sampling
    if N(pop_ind) >= pop_sample(pop_ind)
        ind_sample = sort(randperm( N(pop_ind), pop_sample(pop_ind) ));
    else
        ind_sample = 1:1:N(pop_ind);
    end
    
    for i = 1:length(ind_sample)
        [Y,X,~] = find(spike_hist(ind_sample(i),:));
        xdata = ( [X(:)'; X(:)']+seg_ind(1)-1 )*dt; % sec
        ydata = ([Y(:)'-1;Y(:)']+i-1) + neurons_sampled ;
        hold on;
        line(xdata, ydata, 'Color', 'k', 'linewidth',0.1);
    end
    
    neurons_sampled = neurons_sampled + length(ind_sample); % bookkeeping
end

    
    ylim([0,neurons_sampled]);
    ylabel('Neurons','fontsize',text_fontsize)

    xlim([T(1)-dt, T(1)+(length(seg_ind)-1)*dt]); % make sure all the plots have the same axis scale
    % Keep tick lables while remove tick marks
    % set(gca, 'xtick', [], 'Ticklength', [0 0], 'TickDir','out');
    set(gca,'ytick',[],'ydir','reverse');
    % set(gca,'xtick',[]);
    box on;

end






