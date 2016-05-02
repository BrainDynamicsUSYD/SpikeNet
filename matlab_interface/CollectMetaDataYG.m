function CollectMetaDataYG( varargin )
%
%   Detailed explanation goes here
disp('CollectMetaData...');
tic;

[EI_ratio, loop] = CollectVectorYG('neuron_stats', 'mean(neuron_stats.EI_ratio{1})');
CV2_ISI = CollectVectorYG('Analysis', 'mean(Analysis.CV2_ISI{1}(~isnan(Analysis.CV2_ISI{1})))' );
rate_E = CollectVectorYG('Analysis', 'mean(Analysis.rate{1})' );
rate_I = CollectVectorYG('Analysis', 'mean(Analysis.rate{2})' );

toc;
disp('CollectMetaData Done.');

save('meta_data_tmp.mat');

end


