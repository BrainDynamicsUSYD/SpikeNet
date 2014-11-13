function R = get_cluster_rate(R)
    % kernel for rate estimation
    CC_kernel_width = 50; % ms, kernel length
    choice = 'gaussian';
    kernel = spike_train_kernel_YG(CC_kernel_width, R.reduced.dt, choice);
    
    % cluster mean rate
    C_label = ceil((1:R.N(1))./round(R.N(1)/R.ExplVar.Mnum)); % cluster membership label
    C_rate = zeros(R.ExplVar.Mnum, R.reduced.step_tot);
    for cc = 1:R.ExplVar.Mnum
        C_begin = find(C_label == cc, 1, 'first');
        C_end = find(C_label == cc, 1, 'last');
        C_rate(cc,:) = SpikeTrainConvolve(sum(R.reduced.spike_hist{1}(C_begin:C_end,:), 1)/(C_end-C_begin+1), kernel);
    end
    % record results
    R.cluster.label = C_label;
    R.cluster.rate = C_rate;
    R.cluster.kernel_width = CC_kernel_width;
    R.cluster.kernel_name = choice;
end