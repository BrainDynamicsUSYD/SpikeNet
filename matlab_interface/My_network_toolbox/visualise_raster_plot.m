function visualise_raster_plot( N, dt, spike_hist_cell, num_spikes_cell, num_ref_cell, A0 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


Num_pop = length(num_spikes_cell);
step_tot = length(num_ref_cell{1});
T = (1:step_tot)*dt;

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Get Firing rate and ISI
spike_count = cell(Num_pop,1); % Initialise individual neuron average firing rate
last_spike = cell(Num_pop,1);% Initialise populational Inter-spike-interval (ISI) statistics
ISI = zeros(Num_pop,step_tot);
% Get the data
for pop_ind = 1:Num_pop
    spike_count{pop_ind} = zeros(N(pop_ind),1);
    last_spike{pop_ind} = zeros(N(pop_ind),1);
    hist_ind_current = 0;
    for i = 1:step_tot
        num_spikes_current = num_spikes_cell{pop_ind}(i);
        if num_spikes_current > 0
            Y = spike_hist_cell{pop_ind}(hist_ind_current+1:hist_ind_current+num_spikes_current);
            hist_ind_current = hist_ind_current+num_spikes_current;
            % count the number of spikes for each neuron
            spike_count{pop_ind}(Y) = spike_count{pop_ind}(Y)+1; 
            % ISI
            interval_current = (i-last_spike{pop_ind}(Y)).*double(last_spike{pop_ind}(Y)>0);
            interval_current(interval_current==0) = [];
            ISI(pop_ind,interval_current) = ISI(pop_ind,interval_current)+1;
            last_spike{pop_ind}(Y) = i;
        end
    end
end
%-------------------------------------------------------------------------%
% Plot them
h_ISI = figure('NumberTitle','Off','Name','Average firing rate and ISI distribution', 'position', [300, 100, 768*2, 900]);
for pop_ind = 1:Num_pop
    %---------------------------------------------------------------------%
    % Firing rate (Hz) distribution
    [rate,rate_sort_ind] = sort(spike_count{pop_ind}/T(end)*1000);
    axes_temp = subplot(2,Num_pop,pop_ind);
    xlim([1, N(pop_ind)]);
    if max(rate) > 0 
        ylim([0 max(rate)]);
    end
    xlabel('Sorted neuron index w.r.t. average firing rate');
    ylabel('Average firing rate (Hz)');
    hold on;
    line(1:N(pop_ind), rate, 'Color', 'k', 'LineWidth', 1.5);
    mean_rate = mean(rate);
    line([1, N(pop_ind)], [mean_rate, mean_rate],'Color','k','LineStyle','--','LineWidth', 1.5);% populational average
    
    % Insert histogram
    hist_axes_temp = subplot(2,Num_pop+1,2*(Num_pop+1)); % use some irrelevant axes since it will be deleted anyway
    bin_num = 25;
    hist(rate,bin_num);set(gca, 'ytick',[]);xlabel('Hz');
    hold on;
    insert_plot(h_ISI,axes_temp,hist_axes_temp,0.1,0.7,0.25);

    % Overlap in-degree (in-strength) distribution of pop. No.1 for
    % comparison
    if pop_ind == 1 && nargin == 6
        in_strength = sum(A0,1); % in-strength
        in_strength = in_strength(rate_sort_ind);
        ax2 = axes('Position',get(axes_temp,'Position'),...
           'XAxisLocation','bottom',...
           'YAxisLocation','right',...
           'Color','none',...
           'xtick',[],...
           'XColor','k','YColor','b');
       linkaxes([axes_temp ax2],'x'); hold on;
       % line(1:N(pop_ind),in_strength,'Parent',ax2);
       patchline(1:N(pop_ind),in_strength,'Parent',ax2,'EdgeColor', 'b', 'EdgeAlpha', 0.1, 'LineWidth',0.1);
       
       ylabel('In-strength of Pop. No.1 internal connection');
       ylim([0 max(in_strength)]);
    end

    %---------------------------------------------------------------------%
    % ISI (ms) distribution
    subplot(2,Num_pop,pop_ind+Num_pop);
    line([T; T],[zeros(1,length(T));ISI(pop_ind,:)],'Color','r');
    
    if nnz(ISI(pop_ind,:)) >= 2
        xlim([T(find(ISI(pop_ind,:),1,'first')), T(find(ISI(pop_ind,:),1,'last'))]);
    end
    
    xlabel(['Inter-spike-interval (ms) of Pop. No.',num2str(pop_ind)]);
    ylabel('Number of counts');
end


%--------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Raster plot
figure('NumberTitle','Off','Name','Raster plot','position', [300, 50, 768*2.1, 900]);
axes_matrix = zeros(3,Num_pop);
for pop_ind = 1:Num_pop
    axes_matrix(1,pop_ind) = subplot(6,Num_pop,(0:3)*Num_pop+pop_ind);hold on;
    hist_ind_current = 0;
    % Get new neuron index according to their average firing rate (a bit tricky!!)
    [~,rate_sort_ind] = sort(spike_count{pop_ind});
    [~,sort_back] = sort(rate_sort_ind,'ascend');
    New_ind = 1:N(pop_ind);
    New_ind = New_ind(sort_back);
    
    for i = 1:step_tot
        num_spikes_current = num_spikes_cell{pop_ind}(i);
        if num_spikes_current > 0
            X = ones(1,num_spikes_current)*i*dt;
            Y = spike_hist_cell{pop_ind}(hist_ind_current+1:hist_ind_current+num_spikes_current);
            Y_sorted = New_ind(Y);
            line([X; X],[Y_sorted-0.5; Y_sorted+0.5],'Color','w'); 
            hist_ind_current = hist_ind_current+num_spikes_current;
        end
    end
    ylim([0,N(pop_ind)]);
    if pop_ind == 1
        ylabel('Sorted neuron index w.r.t. average firing rate');
    end
    set(gca, 'xtick', []);
    xlim([0, step_tot*dt]);
    
    % Plot number of spikes
    axes_matrix(2,pop_ind) = subplot(6,Num_pop,4*Num_pop+pop_ind);hold on;
    line([T; T], [zeros(1, step_tot); num_spikes_cell{pop_ind}/N(pop_ind)], 'Color', [255 30 30]/255);
    if pop_ind == 1
        ylabel('% Firing');
    end
    set(gca, 'xtick', []);
    
    % Plot number of refractory neurons
    axes_matrix(3,pop_ind) = subplot(6,Num_pop,5*Num_pop+pop_ind);hold on;
    line([T; T], [zeros(1, step_tot); num_ref_cell{pop_ind}/N(pop_ind)], 'Color', [35 163 200]/255);
    if pop_ind == 1
        ylabel('% Refractory');
    end
    xlabel(['t (ms), pop No.', num2str(pop_ind)]);
    
end

% Background color
set(axes_matrix(:)', 'Color', [0, 0, 0]);

% Link axes to synchronise them when zooming
for pop_ind = 1:Num_pop
    linkaxes(axes_matrix(:,pop_ind),'x');
end

% Set initial y-axes to be same, since only one property can be linked
for i = 2:3
    ylimData = [inf -inf]; % initialise ylim
    % find common ylim
    for pop_ind = 1:Num_pop
        ylimData_new = ylim(axes_matrix(i,pop_ind));
        ylimData(1) = min(ylimData(1),ylimData_new(1));
        ylimData(2) = max(ylimData(2),ylimData_new(2));
    end
    % set common ylim
    for pop_ind = 1:Num_pop
        ylim(axes_matrix(i,pop_ind),ylimData);
    end
end


end

