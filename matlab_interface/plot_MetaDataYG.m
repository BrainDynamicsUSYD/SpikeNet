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
plot(loop, EI_ratio,'o');
ylabel('EI ratio');
subplot(3,1,3);
plot(loop, CV2_ISI.^0.5,'o');
ylabel('CV of ISI');


loop_interesting = loop(rate_E < 50 & CV2_ISI.^0.5 > 0.5 & abs(abs(EI_ratio)-1) < 0.2 )

