function plot_SWR( R, save_figure )
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

dt = R.dt;
step_tot = R.step_tot;

comments = R.comments;

seg_size = 2*10^4; % 2*10^4 for 2-pop, segmentation size for each plot
seg_num = ceil(step_tot/seg_size);
[no, ~] = size(R.LFP.LFP_broad);
for i = 1:no
    for seg = 1:seg_num
        
        seg_ind = get_seg(step_tot, seg_size, seg);
        
        
        
        h_SWR = figure('NumberTitle','Off','Name',strcat('Raster plot:', R.stamp),'units','normalized','position',[0 0 1 1], ...
            'visible', figure_visibility, 'Color','w', 'PaperPositionMode', 'default');
        
        t = (seg_ind-1)*dt*1e-3; % second
        %
        ax1 = subplot(8,1,1);
        LFP_tmp = R.LFP.LFP{1}(i,seg_ind);
        plot(t, LFP_tmp);
        hold on;
        plot(t, ones(size(t))*mean(LFP_tmp),'r');
        plot(t, ones(size(t))*(mean(LFP_tmp)+1*std(LFP_tmp)),'r');
        plot(t, ones(size(t))*(mean(LFP_tmp)+3*std(LFP_tmp)),'r'); 
        % plot(t, ones(size(t))*(mean(LFP_tmp)+5*std(LFP_tmp)),'r');
        ylabel('Unfiltered LFP')
        %
        ax2 = subplot(8,1,2);
        rip_tmp = R.LFP.LFP_ripple(i,seg_ind);
        plot(t,  rip_tmp, 'b');
        ylabel('Rippleband LFP');
        hold on;
        
%         rms_tmp = R.LFP.LFP_ripple_rms(i,seg_ind);
%         rms_scaled = rms_tmp/max(rms_tmp)*max(rip_tmp);
%         plot(t, rms_scaled,'g'); %ylabel('Ripple Hilbert')
        
        hil_tmp = R.LFP.LFP_ripple_hilbert(i,seg_ind);
        hil_scaled = hil_tmp/max(hil_tmp)*max(rip_tmp);
        plot(t, hil_scaled,'r'); %ylabel('Ripple Hilbert')
        
        plot(t, ones(size(t))*R.LFP.ripple_event.hil_mean_baseline(i,end)/max(hil_tmp)*max(rip_tmp) );
        plot(t, ones(size(t))*(R.LFP.ripple_event.hil_mean_baseline(i,end) + R.LFP.ripple_event.no_std*R.LFP.ripple_event.hil_std_baseline(i,end))/max(hil_tmp)*max(rip_tmp) );
      
        
        %
        ax3 = subplot(8,1,3); % Scaleogram with pseudo-Frquency
        freqrange = [R.LFP.lowFreq R.LFP.hiFreq];
        imagesc('XData',t,'YData',R.LFP.wavelet.pseudoFreq,'CData',transpose(R.LFP.wavelet.coeffs{i}(seg_ind,:)));
        ylim(freqrange)
        ylabel('Hz')
        
        
        % raster plot
        ax4 = subplot(8,1,4:7);
        reduced = R.reduced;
        R_LFP.N(1) = sum(R.LFP.LFP_neurons{1}(i,:));
        reduced.spike_hist{1} = reduced.spike_hist{1}(logical((R.LFP.LFP_neurons{1}(i,:))), :);
        R_LFP.reduced = reduced; clear reduced;
        raster_plot(R_LFP, 1, seg, [], 'seg_size', seg_size*(R.dt/R.reduced.dt))
        xlabel('t (sec)');
        % highlight detected ripple events
        if ~isempty(R.LFP.ripple_event.ripple_du_steps{i});
            ripple_start = R.LFP.ripple_event.ripple_start_steps{i};
            ripple_du = R.LFP.ripple_event.ripple_du_steps{i};
            hold on;
            y_lim = get(gca,'ylim');
            for i_r = 1:length(ripple_du)
                tA = ripple_start(i_r);
                tB = tA + ripple_du(i_r) - 1;
                plot([tA*dt tA*dt]*1e-3,y_lim,'r');
                plot([tB*dt tB*dt]*1e-3,y_lim,'r');
            end
        end
        
        
        % Link axes to synchronise them when zooming
        linkaxes([ax1 ax2 ax3 ax4],'x');
        
        
        
        
        % Write comments
        subplot(8,1,8, 'visible','off')
        text(0.5, 0.5, comments, ...
            'VerticalAlignment', 'top', ...
            'HorizontalAlignment', 'center',...
            'FontSize',10,'FontWeight','normal', 'interpreter', 'none'); % ...'interpreter', 'none'... to show underscore
        
        
        % save figure
        if save_figure == 1
            fprintf('\t Saving figure...');
            print(h_SWR, '-dpdf', strcat( R.stamp, '_SWR_',sprintf('%02g_', i), sprintf('%02g', seg)));
            % delete(h_SWR);
            close gcf; clear gcf;
            fprintf('done.\n');
        else
            next = input('\t Next figure?');
            delete(gcf);
        end
    end
end


end










