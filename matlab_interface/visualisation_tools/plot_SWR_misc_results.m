% clc; clear; close all;
% R=load('0009-201608161645-31296_in_1471330157928_out_RYG.mat');
clc;close all;

index1 = cell2mat(R.LFP.ripple_event.index1);
index2 = cell2mat(R.LFP.ripple_event.index2);
index3 = cell2mat(R.LFP.ripple_event.index3);

ripple_Hz = R.LFP.ripple_event.Hz(:)
ripple_du = mean(cell2mat(R.LFP.ripple_event.ripple_du_steps)*0.1')
rate_inside = mean(cell2mat(R.LFP.ripple_event.inside_rate))
rate_outside = mean(cell2mat(R.LFP.ripple_event.outside_rate))
rate_inside_std = std(cell2mat(R.LFP.ripple_event.inside_rate))
rate_outside_std = std(cell2mat(R.LFP.ripple_event.outside_rate))

sw_amp = cell2mat( R.LFP.wavelet.peak.sw_amp(~cellfun(@isempty,R.LFP.wavelet.peak.sw_amp)) );
rp_freq = cell2mat( R.LFP.wavelet.peak.rp_freq(~cellfun(@isempty,R.LFP.wavelet.peak.rp_freq)) );
rp_raw_amp = cell2mat( R.LFP.wavelet.peak.rp_raw_amp(~cellfun(@isempty,R.LFP.wavelet.peak.rp_raw_amp)) );
prn = cell2mat( R.LFP.wavelet.peak.prn(~cellfun(@isempty,R.LFP.wavelet.peak.prn)) );
% rp_wl_amp = cell2mat( R.LFP.wavelet.peak.rp_wl_amp(~cellfun(@isempty,R.LFP.wavelet.peak.rp_wl_amp)) );


figure(6)
title('Ripple and sharp waves')
set(gcf, 'color', 'w');
data1 = []
data2 = []
c = 0;
for ii = 1:4
    for jj = 1:4
        c = c +1;
        subplot(4,4,c);
        hold on;
        axis off;
%         plot(R.LFP.wavelet.rp_average{c}/max(R.LFP.wavelet.rp_average{c}))
%         plot(R.LFP.wavelet.sw_average{c}/max(R.LFP.wavelet.sw_average{c}))
         plot(-750:1:750, - R.LFP.wavelet.rp_average{c})
         plot(-750:1:750, - R.LFP.wavelet.sw_average{c})
         data1 = [data1 max(R.LFP.wavelet.rp_average{c})];
         data2 = [data2 min(R.LFP.wavelet.sw_average{c}(750:end))];
         % xlim([-500 500] )
         % ylim(minmax(-R.LFP.wavelet.sw_average{1}))
         yl = get(gca,'ylim');
         plot([0 0], yl, 'r:');
    end
end




figure(5);
subplot(1,2,1);
rp_amp_no_std = [];
for i = 1:16
    [~,~,sw_amp_mean_base] = find(R.LFP.ripple_event.hil_mean_baseline(i,:), 1, 'last');
    [~,~,sw_amp_std_base] = find(R.LFP.ripple_event.hil_std_baseline(i,:), 1, 'last');
    rp_amp_no_std = [rp_amp_no_std (R.LFP.wavelet.peak.rp_raw_amp{i} - sw_amp_mean_base)/sw_amp_std_base];
end
nostd = 1:15;
c = histc(rp_amp_no_std,nostd );
bar(nostd, c);
xlabel('Ripple magnitude (SD of baseline)');

subplot(1,2,2);
data1 = rp_raw_amp(~isnan(prn));
data2 = prn(~isnan(prn));
plot(data1, data2, 'o');
xlabel('ripple peak magnitude')
ylabel('Post-ripple negativtiy (mininum LFP)')
Post_ripple_negativity = corrcoef(data1(:), data2(:));
Post_ripple_negativity = Post_ripple_negativity(1,2);
title(['r = ' num2str(Post_ripple_negativity)])
set(gcf,'color','w')





figure(4); set(gcf,'color','w')
subplot(1,2,1);
X = [sw_amp(~isnan(sw_amp))', rp_raw_amp(~isnan(sw_amp))'];
[N, C] = hist3(X,[10 10]);
imagesc(C{1}, C{2},N);
set(gca, 'YDir','normal')
xlabel('Sharp-wave amplitude (a.u.)')
ylabel('Ripple magnitude (a.u.)')
a = corrcoef(sw_amp(~isnan(sw_amp)), rp_raw_amp(~isnan(sw_amp)));
title(['r = ' num2str(a(1,2))])

subplot(1,2,2);
X = [sw_amp(~isnan(sw_amp))', rp_freq(~isnan(sw_amp))'];
[N, C] = hist3(X,[10 10]);
imagesc(C{1}, C{2}, N);
set(gca, 'YDir','normal')
xlabel('Sharp-wave amplitude (a.u.)')
ylabel('Ripple frequency (Hz)')
a = corrcoef(sw_amp(~isnan(sw_amp)), rp_freq(~isnan(sw_amp)));
title(['r = ' num2str(a(1,2))])
set(gcf,'color','w')




figure(3)
% Prop. of spikes during SPW-Rs
subplot(1,3,1);
x = index1;
bins = logspace(floor(log10(min(x))), ceil(log10(max(x))), 25);   % Define bins
xc = histc(x,bins);    
plot(bins, xc)
set(gca,'xscale','log')
xlabel('Prop. of spikes during SPW-Rs')
% Prop. of SPW-Rs in which neuron fired
subplot(1,3,2);
x = index2;
bins = logspace(-2, 0, 25);    % Define bins
xc = histc(x,bins);    
plot(bins, xc)
set(gca,'xscale','log')
xlabel('Prop. of SPW-Rs in which neuron fired');
% Mean number of spikes per SPW-R
subplot(1,3,3);
x = index3;
bins = logspace(-2, 1, 25);    % Define bins
xc = histc(x,bins);    
plot(bins, xc)
xlim([10^-2 10])
set(gca,'xscale','log')
xlabel('Mean number of spikes per SPW-R');

set(gcf,'color','w')


figure(2)
subplot(1,3,1)
hist(R.neuron_stats.IE_ratio{1},20);
xlabel('I-E ratio')
title(['mean(IE-ratio) = ' num2str(mean(R.neuron_stats.IE_ratio{1}))])

subplot(1,3,2);
CV_ISI = R.Analysis.CV2_ISI{1}.^0.5;
hist(CV_ISI, 25);
xlabel('CV of ISI')
title(['mean(CV) = ' num2str(mean(R.Analysis.CV2_ISI{1}.^0.5))])


subplot(1,3,3);
hist(R.Analysis.rate{1}, 25);
xlabel('Firing rate (Hz) ')
set(gcf,'color','w')
title(['mean(Hz) = ' num2str(mean(R.Analysis.rate{1}))])

