function peak_check = plot_SWR_multi( R, LFP_grid, I,J,seg_ind )
%   -1 for no figure, 0 for displaying figure, 1 for saving figure
%

%   If the firing history is too long, data will be segmented into several
%   history fractions and plotted separately.
disp('plot_SWR_multi...');


dt = R.dt;
fs = 1/(dt*1e-3); % sampling frequency (Hz)




figure_width = 15; % cm
figure_hight = 20.0; %cm

h_SWR = figure('NumberTitle','off','name','SWR-multi' ,'color','w', ...
    'units', 'centimeters', 'position', [3, 2, 3+figure_width, 2+figure_hight], ...
    'PaperSize', [figure_width, figure_hight] ,'PaperPositionMode','auto');

ax = [];
ax_2 = [];
c_min = inf;
c_max = -inf;
nos = length(I);
peak_check= [];
for ii = 1:nos
    
    ax_tmp = subaxis(nos,2,2*(ii-1)+1,'SV',0.01);
    if ii ~= nos
        axis off;
    else
        set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
    end
    
    ax = [ax ax_tmp];
    hold on;
    
    
    
    
    t = (0:(range(seg_ind)))*dt; % msec
    
    scales = R.LFP.wavelet.scales;
    LFP_tmp = LFP_grid(I(ii),J(ii),:);
    LFP_tmp = LFP_tmp(:)';
    coeffs_tmp = abs(cwt(LFP_tmp,scales,'cmor1.5-1'))';
    coeffs_tmp =  coeffs_tmp(seg_ind,:);
    c_min = min(c_min, min(coeffs_tmp(:)));
    c_max = max(c_max, max(coeffs_tmp(:)));
    
    freqrange = [R.LFP.lowFreq R.LFP.hiFreq];
    imagesc('XData',t,'YData',R.LFP.wavelet.pseudoFreq,'CData',transpose(coeffs_tmp));
    ylim(freqrange)
    
    xlim(minmax(t));
    

        

   
    % Butterworth filter
    order = 4; % 4th order
    lowFreq = 100; % ripple band (default values for this function are 150-250 Hz)
    hiFreq = 250;
    Wn = [lowFreq hiFreq]/(fs/2);
    [b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
    gaus_width = 5; %ms
    
    LFP_ripple = filter(b,a,  LFP_tmp);
    
    xlim(minmax(t));
    
    [ Kernel ] = spike_train_kernel_YG( gaus_width, R.dt, 'gaussian_unit' );
    
    hil_tmp = abs(hilbert( LFP_ripple));
    ripple_power = conv(hil_tmp, Kernel,'same');
    
    ax1 = gca;
    Position = get(ax1,'Position');
    ax2tmp = axes('Position',Position);
    ax_2 = [ax_2  ax2tmp];
    set(ax2tmp, 'visible', 'off')
    
    hold on;
    plot(t,  ripple_power(seg_ind), 'Parent',ax2tmp,'Color','w','LineWidth',1)
    plot(t, LFP_ripple(seg_ind), 'Parent',ax2tmp,'Color','w','LineWidth',1)
    peak_check = [peak_check max(ripple_power(seg_ind))];
   
end
for iii = 1:ii
    caxis(ax(iii), [c_min c_max]);
end

linkaxes(ax_2, 'x');

end












