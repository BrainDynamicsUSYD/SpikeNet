function [bin_count, bin_center] = firing_rate_hist(R, pop)

T = (R.step_tot*R.dt/1000);
firing_rate = full(sum(R.spike_hist{pop},2)/T);

% be careful about neurons that are not firing at all!!!!
fprintf('Number of non-firing neurons: %d \n', nnz(firing_rate == 0))
firing_rate(firing_rate == 0) = []; 

bin_num = 300;


[bin_count, bin_center] = hist(firing_rate, bin_num);

% bin_center = logspace (min(spike_count), max(spike_count), bin_num );
% bin_count = histc(spike_count, bin_center);

bin_count = bin_count/sum(bin_count);
plot(bin_center, bin_count,'.');


%%do some curve fitting!!!
pd = fitdist(firing_rate(:),'lognormal'); % maximum likelihood estimation
%ci = paramci(pd); % confidence interval
% x = linspace(min(spike_count),max(spike_count),1000);
x = bin_center;
y = pdf(pd,x);
% Scale the density by the histogram area for easier display
area = sum(bin_count)*(bin_center(2)-bin_center(1));
hold on;
plot(x,y*area,'r');

% [PARMHAT,PARMCI] = lognfit(X) % another way to fit the distribution

set(gca, 'xscale', 'log');
xlabel('Firing rate (Hz)');
xlim(  10.^[log10(min(bin_center))-1 log10(max(bin_center))+1]  );


end