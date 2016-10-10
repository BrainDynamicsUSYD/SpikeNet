clc;clear;close all;

IE_ratio = CollectVectorYG('neuron_stats', 'neuron_stats.IE_ratio{1}');
CV2_ISI = CollectVectorYG('Analysis', 'Analysis.CV2_ISI{1}(~isnan(Analysis.CV2_ISI{1}))' );
rate_E = CollectVectorYG('Analysis', 'Analysis.rate{1}' );
index1 = CollectVectorYG('LFP', 'cell2mat(LFP.ripple_event.index1)');
index2 = CollectVectorYG('LFP', 'cell2mat(LFP.ripple_event.index2)');
index3 = CollectVectorYG('LFP', 'cell2mat(LFP.ripple_event.index3)');
ripple_Hz = CollectVectorYG('LFP', 'LFP.ripple_event.Hz(:)');
ripple_du = CollectVectorYG('LFP', 'cell2mat(LFP.ripple_event.ripple_du_steps)*0.1');
rate_inside = CollectVectorYG('LFP','cell2mat(LFP.ripple_event.inside_rate)');
rate_outside = CollectVectorYG('LFP','cell2mat(LFP.ripple_event.outside_rate)');
sw_amp = CollectVectorYG('LFP', 'cell2mat( LFP.wavelet.peak.sw_amp(~cellfun(@isempty,LFP.wavelet.peak.sw_amp)) )');
rp_freq = CollectVectorYG('LFP', 'cell2mat( LFP.wavelet.peak.rp_freq(~cellfun(@isempty,LFP.wavelet.peak.rp_freq)) )');
rp_amp = CollectVectorYG('LFP', 'cell2mat( LFP.wavelet.peak.rp_raw_amp(~cellfun(@isempty,LFP.wavelet.peak.rp_raw_amp)) )');
prn = CollectVectorYG('LFP', 'cell2mat(LFP.wavelet.peak.prn(~cellfun(@isempty,LFP.wavelet.peak.prn)) )');
cc_pop = CollectVectorYG('Analysis', 'Analysis.CC_pop{1}' );

rp_amp_no_std = CollectVectorYG('LFP', 'cell2mat( LFP.wavelet.peak.rp_amp_no_std(~cellfun(@isempty,LFP.wavelet.peak.rp_amp_no_std)) )');
rp_amp_mean_baseline= CollectVectorYG('LFP', 'LFP.ripple_event.hil_mean_baseline');
rp_amp_std_baseline = CollectVectorYG('LFP', 'LFP.ripple_event.hil_std_baseline');


save('SWR_metaData.mat') %,'-append')

% load('SWR_metaData.mat')
fprintf('CV of ISI is %0.3f ± %0.3f \n',mean(CV2_ISI.^0.5),std(CV2_ISI.^0.5))
fprintf('IE ratio is %0.3f ± %0.3f \n',mean(IE_ratio),std(IE_ratio))
fprintf('Pair-wise spike train correlation is %0.3f ± %0.3f \n',mean(cc_pop),std(cc_pop))
fprintf('Rate inside SWR is %0.3f ± %0.3f (Hz) \n',mean(full(rate_inside)),std(full(rate_inside)))
fprintf('Rate outside SWR is %0.3f ± %0.3f (Hz) \n',mean(full(rate_outside)),std(full(rate_outside)))
fprintf('SWR duration is %0.3f ± %0.3f (ms) \n',mean(ripple_du),std(ripple_du))
fprintf('SWR event frequency is %0.3f ± %0.3f (Hz) \n',mean(ripple_Hz),std(ripple_Hz))

