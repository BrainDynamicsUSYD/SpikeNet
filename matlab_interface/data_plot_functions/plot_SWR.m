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
    t_1 = R.LFP.transient_steps;
    dt = R.dt;
    step_tot = R.step_tot - t_1;

    comments = R.comments;
    
    seg_size = 2*10^4; % 2*10^4 for 2-pop, segmentation size for each plot
    seg_num = ceil(step_tot/seg_size);
    [no, ~] = size(R.LFP.LFP_broad);
    for i = 1:no
        for seg = 1:seg_num
            
            seg_ind = get_seg(step_tot, seg_size, seg);
            
            
            
            h_SWR = figure('NumberTitle','Off','Name',strcat('Raster plot:', R.stamp),'units','normalized','position',[0 0 1 1], ...
                'visible', figure_visibility, 'Color','w', 'PaperPositionMode', 'default');
            
            t = (t_1+seg_ind-1)*dt*1e-3; % second
            ax1 = subplot(8,1,1);
            LFP_tmp = R.LFP.LFP{1}(i,seg_ind);
            plot(t, LFP_tmp);
            ylabel('Unfiltered LFP')
            ax2 = subplot(8,1,2);
            rip_tmp = R.LFP.LFP_ripple(i,seg_ind);
            plot(t,  rip_tmp, 'b');
            ylabel('Rippleband LFP');
            hold on;
            hil_tmp = R.LFP.LFP_ripple_hilbert(i,seg_ind);
            hil_scaled = hil_tmp/max(hil_tmp)*max(rip_tmp);
            plot(t, hil_scaled,'r--'); %ylabel('Ripple Hilbert')
           
            
            ax3 = subplot(8,1,3); % Scaleogram with pseudo-Frquency
            % x_tmp = LFP_tmp; freqrange = [1 R.LFP.hiFreq];
            x_tmp = rip_tmp; freqrange = [R.LFP.lowFreq R.LFP.hiFreq];
            Fs = 1000/dt;
            fc = centfrq('cmor1.5-1');
           
            scalerange = fc./(freqrange*(1/Fs));
           
            scales = scalerange(end):0.5:scalerange(1);
            pseudoFrq = scal2frq(scales,'cmor1.5-1',1/Fs); % pseudo-frequencies
            Coeffs = cwt(x_tmp,scales,'cmor1.5-1');
            imagesc('XData',t,'YData',pseudoFrq,'CData',abs(Coeffs));
            ylim(freqrange)
            ylabel('Hz')
            
          
            ax4 = subplot(8,1,4:7);
            reduced = R.reduced;
            R_LFP.N(1) = sum(R.LFP.LFP_neurons{1}(i,:));
            reduced.spike_hist{1} = reduced.spike_hist{1}(logical((R.LFP.LFP_neurons{1}(i,:))), :);
            R_LFP.reduced = reduced; clear reduced;
            raster_plot(R_LFP, 1, seg, [], 'seg_size', seg_size*(R.dt/R.reduced.dt))
            xlabel('t (sec)');
              
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


end









