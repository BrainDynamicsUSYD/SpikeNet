% plot_MetaDataYG( varargin )
%
%   Detailed explanation goes here
load('meta_data_tmp.mat');


figure('NumberTitle','off','color', 'w', 'name', 'Basic Meta Data');
subplot(3,1,1);
hold on;
plot(loop, rate_E,'ro');
plot(loop, rate_I,'bo');

subplot(3,1,2);
plot(loop, IE_ratio,'o');
ylabel('IE ratio');
subplot(3,1,3);
plot(loop, CV2_ISI.^0.5,'o');
ylabel('CV of ISI');


loop_interesting = loop(rate_E < 20 & CV2_ISI.^0.5 > 0.5 & abs(IE_ratio) > 0.8 & abs(IE_ratio) < 2 & non_silence > 0.8 )

