function corr_hist_plot(R_temp, pop, varargin)
% pop_ind = 0 for network-wide data sample
hold on;



text_fontsize = 12;
for i = 1:(length(varargin)/2)
    eval([varargin{i*2-1}, '=', num2str(varargin{i*2}) ]);
end


% Dump fields
if pop == 0
    % network wide
    CC = R_temp.Analysis.CC_network;
    CC_kernel_type = R_temp.Analysis.CC_network_kernel_type;
    CC_kernel_width = R_temp.Analysis.CC_network_kernel_width;
else
    CC = R_temp.Analysis.CC_pop{pop};
    CC_kernel_type = R_temp.Analysis.CC_pop_kernel_type;
    CC_kernel_width = R_temp.Analysis.CC_pop_kernel_width;
end

% histogram
bin_num = 100;
bin_edge = linspace(min(CC), max(CC), bin_num);
bin_count = histc(CC,bin_edge);
bin_count = bin_count/sum(bin_count);
h1 = plot(bin_edge,bin_count);
% h1 = bar(bin_edge, bin_count, 'histc');
% set(h1,'facecolor','w','EdgeColor','g', 'linestyle','-');% no filling

% denote mean value
CC_mean = mean(CC);
up_shift = 0.1;
%CC_mean_y = interp1(bin_edge,bin_count,CC_mean) + up_shift*max(bin_count);
CC_mean_y = max(bin_count)*(1+up_shift);
plot(CC_mean, CC_mean_y, 'v'); % "v" for triangle (down)

% % axis
% if max(abs(minmax(bin_edge))) < 0.5
%     xlim([-0.5,0.5]);
% elseif max(abs(minmax(bin_edge))) < 1
%     xlim([-1, 1]);
% end
xlim([-1, 1]);

% % remove y axis
% set(gca,'ytick',[], 'YColor','w','box','off', 'TickDir','out');

set(gcf, 'InvertHardCopy', 'off'); % prevent the x-axis line to reappear when printed

% display information
if pop > 0
    fprintf('Spike-train correlation coefficient distribution with mean value %g \n %g pairs sampled from pop %g, %s with T = %g (ms)\n', CC_mean, length(CC), pop, CC_kernel_type, CC_kernel_width);
elseif pop == 0
    fprintf('Spike-train correlation coefficient distribution with mean value %g \n %g pairs sampled from entire network, %s with T = %g (ms)\n', CC_mean, length(CC), CC_kernel_type, CC_kernel_width);
end

end


