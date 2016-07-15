function R = cluster_sorted_rate( R )

    % get cluster rate
    R = get_cluster_rate( R );
    
    % dump constants
    dt = R.reduced.dt;
    step_tot = R.reduced.step_tot;
    rate = R.cluster.rate;
    Mnum = R.ExplVar.Mnum;
    
    % sort the rate and each time point
    [sorted_rate, cluster_sequence] = sort(rate, 'descend');
    
    % cluster rank
    [a_ind, b_ind] = size(cluster_sequence);
    cluster_rank = zeros(a_ind,b_ind);
    for i = 1:b_ind
        cluster_rank(cluster_sequence(:,i),i) = (1:a_ind)';
    end

    
    % symbolic autocorrelation
    sample_step = 10;
    t_sample = 1:sample_step:step_tot; % down sampling the data!
    dt_sample = dt*sample_step;
    lag = 1:1:round(length(t_sample)/3);
    symb_acc = [];
    lag_ms = lag*dt_sample;
    for i = 1:Mnum
        [cc_tmp] = es_acf(cluster_sequence(i,t_sample),lag); % this may take some time
        symb_acc = [symb_acc; cc_tmp];
    end
    
    % correlation matrix for sort rate
    cc_sorted_rate = corrcoef(sorted_rate');
    cc_rate = corrcoef(rate');
    
    
    % threshold the 1st symbolic sequence based on 1st sorted rate
    % this should be done in 3-steps
    
    % step 1: domination requirement
    % normalize the sorted rates at every time step
    % and the largest component must be larger than some threhold
    % this is essentially testing if there is any truly dominating
    % component at all.
    % note that the normalized components are the cos(.) values of the
    % angles between the vector and each axis!
    cos_thre = 0.8;
    sorted_rate_norm = normc(sorted_rate); % normalize each column
    seq_1st_dominant = cluster_sequence(1,:);
    seq_1st_dominant(sorted_rate_norm(1,:) < cos_thre ) = 0;
    
    
    % step 2: absolute firing requirement
    % test if the dominating component has a firing
    % rate high enough
    Hz_thre = 5:1:10; % Hz, threshold
    lt = length(Hz_thre);
    seq_1st = repmat(seq_1st_dominant, lt, 1); 
    for i = 1:lt
        seq_1st(i, sorted_rate(1,:) < Hz_thre(i)  ) = 0; % or NaN?
    end
    
    % step 3: persistence requirement
    % apply persistence requirement on the symbolic sequence
    high_du_min = 200; % ms
    high_du_min_steps =  round(high_du_min/dt);
    for i = 1:lt
        seq_1st(i, :) = persistence_requirement( seq_1st(i, :), high_du_min_steps );
    end
    
    
    % statistical analysis on the thresholded 1st symbolic sequence
    switch_seq = cell(1,lt);
    switch_freq = zeros(1,lt);
    H2 = zeros(1,lt);
    high_du = cell(1,lt);
    low_du = cell(1,lt);
    switch_level = cell(1,lt);
    % order_para = zeros(1,lt);
    scope = zeros(1,lt); % scope is the number of different clusters the system visits
    for i = 1:lt
        [seq_tmp, high_du_tmp, low_du_tmp] = seq_postprocess( seq_1st(i,:), dt );
        
        switch_level{i} = switch_level_detect(seq_tmp);
        switch_seq{i} = seq_tmp;
        high_du{i} = high_du_tmp;
        low_du{i} = low_du_tmp;
        
        % switch frequency (not a very goof order parameter)
        switch_freq(i) = length(seq_tmp)/(sum([high_du_tmp low_du_tmp])*10^-3); % in Hz
        
        % switching dynamics order parameter
        H2(i) = symbolic_block_entropy(seq_tmp, Mnum);
        scope(i) = length( unique(seq_tmp) );
        % order_para(i) = H2(i)*switch_freq(i);
    end
    
    
    % output results
    R.cluster.sorted_rate = sorted_rate;
    R.cluster.sym_seq = cluster_sequence;
    R.cluster.rate_rank = cluster_rank;
    
    R.cluster.symb_acc = symb_acc;
    R.cluster.acc_lag = lag_ms;
    R.cluster.cc_sorted_rate = cc_sorted_rate;
    R.cluster.cc_rate = cc_rate;
    
    R.cluster.high_du_min = high_du_min;
    R.cluster.cos_thre = cos_thre;
    R.cluster.Hz_thre = Hz_thre;
    R.cluster.switch_freq = switch_freq;
    R.cluster.scope = scope;
    R.cluster.H2 = H2;
    % R.cluster.order_para = order_para;
    R.cluster.switch_seq = switch_seq;
    R.cluster.switch_level = switch_level;
    R.cluster.high_du = high_du;
    R.cluster.low_du = low_du;
    R.cluster.sym_seq_1st_theta = seq_1st;
    

    
end



function [switch_level] = switch_level_detect(seq)

level_num = ceil(log2(max(seq))) + 1; % total number of hierarchy levels
seq = seq - 1; % zero based index

% convert to binary strings
seq_bin = cell(size(seq));
for i = 1:length(seq)
    seq_bin{i} = dec2bin(seq(i), level_num);
end

% find level based on binary strings
switch_level = zeros(1, length(seq)-1);
for i = 1:(length(seq)-1)
    cmp = (seq_bin{i} == seq_bin{i+1});
    for j = 2:level_num
        if cmp(j) == 1 && cmp(j-1) == 0
            cmp(j) = 0;
        end
    end
    switch_level(i) = (level_num+1) - find(cmp,1,'last'); 
end

end








