function HistogramsYG( Result_cell, save_figure )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
disp('HistogramsYG...');
tic;

if nargin == 0
    Result_cell = CollectRYG();
end

if nargin <= 1
    save_figure = 1; % default
end

if save_figure == 1
    figure_visibility = 'off'; % 'on', 'off'
else
    figure_visibility = 'on';
end

Result_num = length(Result_cell);

for r_num = 1:Result_num
    R_temp = Result_cell{r_num};
    
    % Dump fields
    Num_pop = R_temp.Num_pop;
    ISI_dist = R_temp.Analysis.ISI_dist;
    CV_ISI = R_temp.Analysis.CV_ISI;
    CV_ISI_overall = R_temp.Analysis.CV_ISI_overall;
    CC_network_sample = R_temp.Analysis.CC_network_sample;
    CC_network_sample_mean = R_temp.Analysis.CC_network_sample_mean;
    CC_kernel_type = R_temp.Analysis.CC_kernel_type;
    CC_kernel_width = R_temp.Analysis.CC_kernel_width;
    
    %-------------------------------------------------------------------------%
    % Plot
    h_hist = figure('NumberTitle','Off','Name','Histograms','units','normalized','position',[0 0 1 1], ...
        'visible', figure_visibility, 'Color','w', 'PaperPositionMode', 'default');
    axes_matrix = zeros(3,Num_pop);
    for pop_ind = 1:Num_pop
        if isempty(ISI_dist{pop_ind}) % if empty, skip plotting
            break;
        end
        %%%%%%%%%%% semi-log ISI histogram using bar
        axes_matrix(1,pop_ind) = subplot(3,Num_pop,0*Num_pop+pop_ind);hold on;
        bin_num = 100;
        bin_edge = logspace(floor(log10(min(ISI_dist{pop_ind}))),ceil(log10(max(ISI_dist{pop_ind}))), bin_num);
        bin_count = histc(ISI_dist{pop_ind},bin_edge);
        bar(bin_edge, bin_count, 'histc');
        set(gca, 'Xscale', 'log');
        set(findobj(gca,'Type','line'),'Marker','none'); % prevent the asterisks from being displayed on unequal bins in plots created with BAR
        
        xlimA = 10^floor(log10(min(ISI_dist{pop_ind})));
        xlimB = 10^ceil(log10(max(ISI_dist{pop_ind})));
        if ~isnan(xlimA) && ~isnan(xlimB) 
            if xlimA < xlimB
                xlim([xlimA xlimB]);
            end
        end
        
        xlabel(['ISI (ms), pop ' num2str(pop_ind)]);
        set(gca,'ytick',[],'Ycolor','w','box','off')
        
        %%%%%%%%%%% CV-ISI
        axes_matrix(2,pop_ind) = subplot(3,Num_pop,1*Num_pop+pop_ind);hold on;
        
        bin_num = 100;
        bin_edge = logspace(floor(log10(min(CV_ISI{pop_ind}))),ceil(log10(max(CV_ISI{pop_ind}))), bin_num);
        bin_count = histc(CV_ISI{pop_ind},bin_edge);
        bar(bin_edge, bin_count, 'histc');
        set(gca, 'Xscale', 'log');
        set(findobj(gca,'Type','line'),'Marker','none'); % prevent the asterisks from being displayed on unequal bins in plots created with BAR
        
        xlimA = 10^floor(log10(min(CV_ISI{pop_ind})));
        xlimB = 10^ceil(log10(max(CV_ISI{pop_ind})));
        if ~isnan(xlimA) && ~isnan(xlimB)
            if xlimA < xlimB
                xlim([xlimA xlimB]);
            end
        end
        
        label_temp = sprintf('Coefficient of variation of ISI, pop %g \n network mean value %g', pop_ind, CV_ISI_overall);
        xlabel(label_temp);
        set(gca,'ytick',[],'Ycolor','w','box','off')
        
        %%%%%%%%%%% spike train CC sample
        axes_matrix(3,pop_ind) = subplot(3,Num_pop,2*Num_pop+(1:Num_pop));hold on;
        if pop_ind == 1 % only need to plot once
           network_corr_plot(R_temp);
        end
    end

    % save figure
    if save_figure == 1
        fprintf('\t Saving figure...');
        print(h_hist, '-pdf', strcat('data/', R_temp.name{1}, '_hist')); 
        delete(h_hist);
        fprintf('Saving done.\n');
    else
        next = input('\t Next figure?');
        delete(h_hist);
    end
end
toc;

end








