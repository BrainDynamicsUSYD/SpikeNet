function R = get_cluster_rate(R, CC_kernel_width)

    Mnum = R.ExplVar.Mnum;
    N = R.N;
    dt =  R.reduced.dt;
    step_tot = R.reduced.step_tot;
    spike_hist = R.reduced.spike_hist;

    % kernel for rate estimation
    if nargin == 1
        CC_kernel_width = 50; % ms, kernel length
    end
    choice = 'gaussian';
    kernel = spike_train_kernel_YG(CC_kernel_width, dt, choice);
    

    
    % cluster mean rate
    C_label = ceil((1:N(1))./round(N(1)/Mnum)); % cluster membership label
    C_rate = zeros(Mnum, step_tot);
    for cc = 1:Mnum
        C_begin = find(C_label == cc, 1, 'first');
        C_end = find(C_label == cc, 1, 'last');
        C_rate(cc,:) = SpikeTrainConvolve(sum(spike_hist{1}(C_begin:C_end,:), 1)/(C_end-C_begin+1), kernel);
    end
    % record results
    R.cluster.label = C_label;
    R.cluster.rate = C_rate;
    R.cluster.kernel_width = CC_kernel_width;
    R.cluster.kernel_name = choice;
end