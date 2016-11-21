function plot_SWR_multi( R, save_figure, no_given, seg_given )
%   -1 for no figure, 0 for displaying figure, 1 for saving figure
%
%   If the firing history is too long, data will be segmented into several
%   history fractions and plotted separately.
disp('plot_SWR_multi...');


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



seg_size = 0.5*10^4; % 2*10^4 for 2-pop, segmentation size for each plot
seg_num = ceil(step_tot/seg_size);
[no, ~] = size(R.LFP.LFP_broad);
nos = 1:no;

if nargin >= 3
    nos = no_given;
end

if nargin == 4
    seg_num = 1;
end
for seg = 1:seg_num
    
    ii = 0;
    
    figure_width = 8; % cm
    figure_hight = 20.0; %cm
    
    h_SWR = figure('NumberTitle','off','name','SWR-multi' ,'color','w', ...
        'units', 'centimeters', 'position', [3, 2, 3+figure_width, 2+figure_hight], ...
        'PaperSize', [figure_width, figure_hight] ,'PaperPositionMode','auto');
    
    ax = [];
    c_min = inf;
    c_max = -inf;
    for i = nos
        ii = ii + 1;
        ax_tmp = subaxis(no,1,ii,'SV',0.01);
        if ii ~= length(nos)
        axis off;
        end
        ax = [ax ax_tmp];
        hold on;
        
        seg_ind = get_seg(step_tot, seg_size, seg);
        if nargin == 4
            seg_ind = seg_given;
        end
        
        
        
        t = (seg_ind-1)*dt*1e-3; % second
        %
        
        
        %         rip_tmp = R.LFP.LFP_ripple(i,seg_ind);
        %         plot(t,  rip_tmp, 'w');
        
        
        %         rms_tmp = R.LFP.LFP_ripple_rms(i,seg_ind);
        %         rms_scaled = rms_tmp/max(rms_tmp)*max(rip_tmp);
        %         plot(t, rms_scaled,'g'); %ylabel('Ripple Hilbert')
        
        %         hil_tmp = R.LFP.LFP_ripple_hilbert(i,seg_ind);
        %         hil_scaled = hil_tmp/max(hil_tmp)*max(rip_tmp);
        %         plot(t, hil_scaled,'r'); %ylabel('Ripple Hilbert')
        %
        %         plot(t, ones(size(t))*R.LFP.ripple_event.hil_mean_baseline(i,end)/max(hil_tmp)*max(rip_tmp) );
        %         plot(t, ones(size(t))*(R.LFP.ripple_event.hil_mean_baseline(i,end) + R.LFP.ripple_event.no_std*R.LFP.ripple_event.hil_std_baseline(i,end))/max(hil_tmp)*max(rip_tmp) );
        %
        
        %
        scales = R.LFP.wavelet.scales;
        x_tmp = R.LFP.LFP_ripple(i,seg_ind);
        coeffs_tmp = abs(cwt(x_tmp,scales,'cmor1.5-1'))';
        c_min = min(c_min, min(coeffs_tmp(:)));
        c_max = max(c_max, max(coeffs_tmp(:)));
        
        freqrange = [R.LFP.lowFreq R.LFP.hiFreq];
        imagesc('XData',t,'YData',R.LFP.wavelet.pseudoFreq,'CData',transpose(coeffs_tmp));
        ylim(freqrange)
        
        xlim(minmax(t));
        
        
    end
    for iii = 1:ii
        caxis(ax(iii), [c_min c_max]);
    end
    
    
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










