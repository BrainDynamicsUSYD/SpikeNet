function plot_stPR( R )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

figure('numbertitle','off','name','stPR','color', 'w');

subplot(2,1,1);
plot(R.stPR.lags, R.stPR.stPR_full);
xlabel('ms')
ylabel('Spikes per sec');


subplot(2,1,2);
plot(R.stPR.lags, R.stPR.stPR_full_shuffle);
xlabel('ms')
ylabel('Spikes per sec');

mm = minmax(R.stPR.c_norm(:)');
edges = linspace(mm(1), mm(2), 15);
Y1 = histc(R.stPR.c/nanmedian(R.stPR.c_shuffle), edges);
Y2 = histc(R.stPR.c_shuffle/nanmedian(R.stPR.c_shuffle), edges);

inset();
hold on;
bar(edges, Y1, 'b');
bar(edges, Y2, 'g');
xlim([floor(mm(1)), ceil(mm(2))]);

end

