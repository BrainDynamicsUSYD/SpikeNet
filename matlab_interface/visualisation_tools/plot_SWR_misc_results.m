function plot_SWR_misc_results(R)
% clc; clear; close all;
% R=load('0009-201608161645-31296_in_1471330157928_out_RYG.mat');
clc;close all;
text_fontsize = 6;

index1 = cell2mat(R.LFP.ripple_event.index1);
index2 = cell2mat(R.LFP.ripple_event.index2);
index3 = cell2mat(R.LFP.ripple_event.index3);

IE_ratio = R.neuron_stats.IE_ratio{1};
CV2_ISI = R.Analysis.CV2_ISI{1}(~isnan(R.Analysis.CV2_ISI{1}))';
rate_E =  R.Analysis.rate{1}';

ripple_Hz = R.LFP.ripple_event.Hz(:)';
ripple_du = cell2mat(R.LFP.ripple_event.ripple_du_steps)*0.1';
rate_inside = cell2mat(R.LFP.ripple_event.inside_rate)';
rate_outside = cell2mat(R.LFP.ripple_event.outside_rate)';

fprintf('Firing rate of exitatory neurons is %0.3f ± %0.3f \n',mean(rate_E),std(rate_E))
fprintf('CV of ISI is %0.3f ± %0.3f \n',mean(CV2_ISI.^0.5),std(CV2_ISI.^0.5))
fprintf('IE ratio is %0.3f ± %0.3f \n',mean(IE_ratio),std(IE_ratio))
fprintf('Rate inside SWR is %0.3f ± %0.3f (Hz) \n',mean(full(rate_inside)),std(full(rate_inside)))
fprintf('Rate outside SWR is %0.3f ± %0.3f (Hz) \n',mean(full(rate_outside)),std(full(rate_outside)))
fprintf('SWR duration is %0.3f ± %0.3f (ms) \n',mean(ripple_du),std(ripple_du))
fprintf('SWR event frequency is %0.3f ± %0.3f (Hz) \n',mean(ripple_Hz),std(ripple_Hz))



sw_amp = cell2mat( R.LFP.wavelet.peak.sw_amp(~cellfun(@isempty,R.LFP.wavelet.peak.sw_amp)) );
rp_freq = cell2mat( R.LFP.wavelet.peak.rp_freq(~cellfun(@isempty,R.LFP.wavelet.peak.rp_freq)) );
rp_raw_amp = cell2mat( R.LFP.wavelet.peak.rp_raw_amp(~cellfun(@isempty,R.LFP.wavelet.peak.rp_raw_amp)) );
prn = cell2mat( R.LFP.wavelet.peak.prn(~cellfun(@isempty,R.LFP.wavelet.peak.prn)) );


figure('numbertitle','off','name','SWR misc results', 'color','w','position', [88   229   604   856])
subplot(4,3,7);
X = [sw_amp(~isnan(sw_amp))', rp_raw_amp(~isnan(sw_amp))'];
[N, C] = hist3(X,[10 10]);
imagesc(C{1}, C{2},N);
set(gca, 'YDir','normal')
xlabel('Sharp-wave amplitude (a.u.)')
ylabel('Ripple magnitude (a.u.)')
a = corrcoef(sw_amp(~isnan(sw_amp)), rp_raw_amp(~isnan(sw_amp)));
title(['r = ' num2str(a(1,2))])

subplot(4,3,8);
X = [sw_amp(~isnan(sw_amp))', rp_freq(~isnan(sw_amp))'];
[N, C] = hist3(X,[10 10]);
imagesc(C{1}, C{2}, N);
set(gca, 'YDir','normal')
xlabel('Sharp-wave amplitude (a.u.)')
ylabel('Ripple frequency (Hz)')
a = corrcoef(sw_amp(~isnan(sw_amp)), rp_freq(~isnan(sw_amp)));
title(['r = ' num2str(a(1,2))])
set(gcf,'color','w')


subplot(4,3,9);
title('Ripple and sharp waves');
hold on;
plot(mean(cell2mat(R.LFP.wavelet.rp_average')));
plot(mean(cell2mat(R.LFP.wavelet.sw_average')));




subplot(4,3,10);
nostd = 1:20;
rp_amp_no_std = cell2mat( R.LFP.wavelet.peak.rp_amp_no_std(~cellfun(@isempty, R.LFP.wavelet.peak.rp_amp_no_std)) )';
c = histc( rp_amp_no_std,nostd );
bar(nostd, c);
xlabel('Ripple magnitude (SD of baseline)');
xlim(minmax(nostd));

subplot(4,3,11);
data1 = rp_raw_amp(~isnan(prn));
data2 = prn(~isnan(prn));
plot(data1, data2, 'o');
xlabel('ripple peak magnitude')
ylabel('Post-ripple negativtiy (mininum LFP)')
Post_ripple_negativity = corrcoef(data1(:), data2(:));
Post_ripple_negativity = Post_ripple_negativity(1,2);
title(['r = ' num2str(Post_ripple_negativity)])
set(gcf,'color','w')


% Prop. of spikes during SPW-Rs
ax_tmp = subplot(4, 3, 1 );hold on;
%add_figure_label( 'A', anno_shift1, anno_fontsize );
%set(gca,'fontsize',tick_fontsize);
x = index1(index1 ~= 0);
bins = logspace(floor(log10(min(x))), ceil(log10(max(x))), 50);   % Define bins
xc = histc(x,bins);    
plot(bins, xc,'k.','markersize',5);
xlim([10^floor(log10(min(x))), 10^ceil(log10(max(x)))])
set(gca,'xscale','log','xtick', 10.^(-7:2:-3))
xlabel('Prop. of spikes during SPW-Rs','fontsize',text_fontsize)
% logn-fit
[parmhat,parmci] = lognfit(x);
y_f = lognpdf(bins,parmhat(1)+parmhat(2)^2,parmhat(2));% log-bin correction to parameter mu!!
plot(bins,y_f/sum(y_f)*sum(xc),'b-')


% Prop. of SPW-Rs in which neuron fired
ax_tmp = subplot(4, 3, 2 );hold on;
%add_figure_label( 'B', anno_shift1, anno_fontsize );
%set(gca,'fontsize',tick_fontsize);
x = index2(index2>0);
bins = logspace(-2, 0, 50);    % Define bins
xc = histc(x,bins);    
plot(bins, xc,'k.','markersize',5)
set(gca,'xscale','log','xtick', 10.^(-2:1:0))
xlabel('Prop. of SPW-Rs in which neuron fired','fontsize',text_fontsize);
% logn-fit
[parmhat,parmci] = lognfit(x);
y_f = lognpdf(bins,parmhat(1)+parmhat(2)^2,parmhat(2)); % log-bin correction to parameter mu!!
plot(bins,y_f/sum(y_f)*sum(xc),'b-')

% Mean number of spikes per SPW-R
ax_tmp = subplot(4, 3, 3 );hold on;
%add_figure_label( 'C', anno_shift1, anno_fontsize );
%set(gca,'fontsize',tick_fontsize);
x = index3(index3>0);
bins = logspace(-2, 1, 50);    % Define bins
xc = histc(x,bins);    
plot(bins, xc,'k.','markersize',5)
xlim([10^-2 10])
set(gca,'xscale','log','xtick',10.^(-2:1:1))
xlabel('Mean number of spikes per SPW-R','fontsize',text_fontsize);
% logn-fit
[parmhat,parmci] = lognfit(x);
y_f = lognpdf(bins,parmhat(1)+parmhat(2)^2,parmhat(2)); % log-bin correction to parameter mu!!
plot(bins,y_f/sum(y_f)*sum(xc),'b-')
set(gcf,'color','w')



subplot(4,3,4)
hist(R.neuron_stats.IE_ratio{1},20);
xlabel('I-E ratio')
title(['mean(IE-ratio) = ' num2str(mean(R.neuron_stats.IE_ratio{1}))])

subplot(4,3,5);
hist(CV2_ISI.^0.5, 25);
xlabel('CV of ISI')
title(['mean(CV) = ' num2str(mean(CV2_ISI.^0.5))])


subplot(4,3,6);
hist(R.Analysis.rate{1}, 25);
xlabel('Firing rate (Hz) ')
set(gcf,'color','w')
title(['mean(Hz) = ' num2str(mean(R.Analysis.rate{1}))])


x = cumsum(cos(R.grid.jump_dir) .* R.grid.jump_dist);
y = cumsum(sin(R.grid.jump_dir) .* R.grid.jump_dist);
subplot(4,3,12);
plot(x, y)




trajectory_animation(x,y,0.01, 1000);




cc_pop = R.Analysis.CC_pop{1}';
fprintf('Pair-wise spike train correlation is %0.3f ± %0.3f \n',mean(cc_pop),std(cc_pop))

end

