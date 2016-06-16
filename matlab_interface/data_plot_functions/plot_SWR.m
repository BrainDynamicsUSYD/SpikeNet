function plot_SWR( Result_cell, save_figure )
%   -1 for no figure, 0 for displaying figure, 1 for saving figure
%
%   If the firing history is too long, data will be segmented into several
%   history fractions and plotted separately.
disp('plot_SWR...');


% if nargin == 0
%     Result_cell = CollectRYG();
% end
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
    R = Result_cell{r_num};
    dt = R.dt;
    step_tot = R.step_tot;
    %LFP = R.LFP{1};
    SWR = R.SWR;
    comments = R.comments;
    
    seg_size = 4*10^4; % 2*10^4 for 2-pop, segmentation size for each plot
    seg_num = ceil(step_tot/seg_size);
    for seg = 1:seg_num
        
        seg_ind = get_seg(step_tot, seg_size, seg);
        
        h_SWR = figure('NumberTitle','Off','Name',strcat('Raster plot:', R.stamp),'units','normalized','position',[0 0 1 1], ...
            'visible', figure_visibility, 'Color','w', 'PaperPositionMode', 'default');
        
        t = seg_ind*dt*1e-3; % second
        ax1 = subplot(3,1,1);
        plot(t, SWR.LFP_broad(seg_ind));
        ylabel('Broadband LFP')
        ax2 = subplot(3,1,2);
        plot(t, SWR.LFP_ripple(seg_ind));
        ylabel('Rippleband LFP')
        xlabel('t (sec)');
        
        % Link axes to synchronise them when zooming
        linkaxes([ax1 ax2],'x');

        

        % Write comments
        subplot(3,1,3, 'visible','off')
        text(0.5, 0.5, comments, ...
            'VerticalAlignment', 'top', ...
            'HorizontalAlignment', 'center',...
            'FontSize',10,'FontWeight','normal', 'interpreter', 'none'); % ...'interpreter', 'none'... to show underscore

        
        % save figure
        if save_figure == 1
            fprintf('\t Saving figure...');
            print(h_SWR, '-dpdf', strcat( R.stamp, '_SWR_',sprintf('%02g', seg)));
            % delete(h_SWR);
            close gcf; clear gcf;
            fprintf('done.\n');
        else
            next = input('\t Next figure?');
            delete(h_raster);
        end
    end
end


end










