clear;
%   Detailed explanation goes here
disp('CollectMetaData...');
tic;

[IE_ratio, loop] = CollectVectorYG('neuron_stats', 'mean(neuron_stats.IE_ratio{1})');
CV2_ISI = CollectVectorYG('Analysis', 'mean(Analysis.CV2_ISI{1}(~isnan(Analysis.CV2_ISI{1})))' );
rate_E = CollectVectorYG('Analysis', 'mean(Analysis.rate{1})' );
rate_I = CollectVectorYG('Analysis', 'mean(Analysis.rate{2})' );
non_silence = CollectVectorYG('reduced', 'nnz(reduced.num_spikes{1})/length(reduced.num_spikes{1})');

loop_interesting = loop(rate_E < 20 & CV2_ISI.^0.5 > 0.5 & abs(IE_ratio) > 0.8 & non_silence > 0.8 );
interesting_definition = 'loop_interesting = loop(rate_E < 50 & CV2_ISI.^0.5 > 0.5 & abs(IE_ratio) > 0.8 & non_silence > 0.8 );';

toc;
disp('CollectMetaData Done.');

save('meta_data_tmp.mat');


ripple_Hz = CollectVectorYG('LFP', 'mean(LFP.ripple_event.Hz(:))');
rp_freq = CollectVectorYG('LFP', 'mean(cell2mat( LFP.wavelet.peak.rp_freq(~cellfun(@isempty,LFP.wavelet.peak.rp_freq)) ))');
