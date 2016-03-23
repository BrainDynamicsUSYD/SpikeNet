function [I, J, K, N, iter_hist] = spatial_embed_in_out_network(hw, P0_init, ...
    degree_CV, tau_d, in_out_corrcoef, K_mu_log, K_sigma_log, cn_scale_wire, cn_scale_weight)
% [I, J, K, N, iter_hist] = spatial_embed_in_out_network(hw, P0_init, ...
%    degree_CV, tau_d, in_out_corrcoef, K_mu_log, K_sigma_log, cn_scale_wire, cn_scale_weight)
% 
% P0: average connection probability
% The final P0 will be different from P0_init
% degree_CV = degree_dist_std/degree_dist_mean
%
% tau_d: decay constant for the exponential distribution (units: neurons)
%
% default: log-normal digree distribution and exponential distance wiring
% rule
%

% hw = 21; % hw = 31; % (31*2+1)^2 = 3969 ~ 4000
% P0_init = 0.1;
% degree_CV = 0.1;
% tau_d = 10; %
% in_out_corrcoef = 0.3;


D = 2; % Pulin: 2D is enough
[Lattice, N] = lattice_nD(D, hw);
fprintf('The number of node is %d. \n' , N);

% % degree distribution
% log-normal
deg_mean = N*P0_init;
deg_std = degree_CV*deg_mean;
[ degree_in_0, degree_out_0 ] = logn_in_out_degree( N, [deg_mean deg_mean], [deg_std,  deg_std], in_out_corrcoef  );

%
for iter = 1:5 % 5 is arbitrary, try something else?
    I = []; % pre node index
    J = []; % post node index
    K = [];
    degree_in_left = degree_in_0;
    for i = randperm(N)
        % distance factor
        post_dist = lattice_nD_find_dist(Lattice, hw, i);
        dist_factor = exp(-post_dist/tau_d); % exponential distribution, will be properly scaled later;
        dist_factor(i) = 0; % no self-connection
        % common neighbour factor
        if iter == 1
            cn_factor = ones(size(dist_factor));
        else
            cn_factor = 1 + (cn(:,i) - cn_minmax(1))/diff(cn_minmax)*(cn_scale_wire - 1);
        end
        % degree_in factor
        degree_in_factor = degree_in_left;
        % joint factor
        %dist_factor = normc(dist_factor);
        %degree_in_factor = normc(degree_in_factor);
        %cn_factor = normc(cn_factor);
        joint_factor = dist_factor.*degree_in_factor.*cn_factor; % will be properly scaled later;
        joint_factor = normc(joint_factor);
        
        % establish connections
        [~, ind] = sort( rand(N,1)./joint_factor, 'ascend' );
        chosen_j = ind(1: degree_out_0(i) );
%         Issue: sometimes what's left will be [NaN NaN 1 1 2 1 NaN] and 5
%         connections are needed. Try to fix this?
%         if sum(isnan(degree_in_left(chosen_j))) ~= 0
%             warning('sth wrong here')
%         end
        
        degree_in_left(chosen_j) = degree_in_left(chosen_j) - 1;
        degree_in_left(degree_in_left == 0) = NaN;
        
        % store results
        chosen_j = chosen_j(:); % column vector
        I = [I; i*ones(size(chosen_j)) ]; %#ok<AGROW>
        J = [J; chosen_j]; %#ok<AGROW>
        K_tmp = post_dist(chosen_j);
        K = [K; K_tmp(:)]; %#ok<AGROW>
    end

    
    %%%%% find number of common pre-synaptic neighbours
    A = sparse(I,J, ones(size(I))); % I:pre, J:post
    cn = A'* A; % A * A' gives common post-synaptic neighbours, note that cn is symmetric
    cn(logical(eye(size(cn)))) = 0; % no self-connection allowed
    % pre-calculate the global common neighbour scaling for the next
    % iteration
    cn_dist = full(triu(cn,1));
    cn_dist = cn_dist(:)';
    edges =  0:max(cn_dist);
    Y = histc(cn_dist, edges);
    cn_minmax = minmax(cn_dist);
    clear cn_dist;
    % store iteration history
    if iter == 1
        iter_hist.edges = {edges};
        iter_hist.Y = {Y};
        iter_hist.cn_minmax = {cn_minmax};
    else
        iter_hist.edges{end+1} = edges;
        iter_hist.Y{end+1} = Y;
        iter_hist.cn_minmax{end+1} = cn_minmax;
    end
%     plot(edges,Y)
%     mean(cn_dist)
%     std(cn_dist)

end

%%%% Generate K

% define and in_degree vs K_neuron_sum scaling 
K_tot = K_mu_log*length(I);
in_degree = full(sum(A, 1));
scale = sqrt(in_degree);
K_neuron_sum_scale = scale*K_tot/sum(scale);
iter_hist.K_neuron_sum_scale = K_neuron_sum_scale;

plot_results = 0;
[ K ] = get_K_given_certain_conditions( I, J, K_mu_log, K_sigma_log, K_neuron_sum_scale, plot_results);

K_minmax = minmax(K);
% now shuffle the K according to common neighbour rule
for i = 1:N
    c_tmp = cn(i, J(I == i));
    c_tmp = c_tmp(:)';
    K_tmp = K(I == i);
    K_tmp = K_tmp(:)';
    K_tmp_shuffle = zeros(size(K_tmp));
    % bilinear factor
    ck_factor = 1 + (c_tmp - cn_minmax(1))./diff(cn_minmax) .* ...
        (K_tmp - K_minmax(1))./diff(K_minmax) * (cn_scale_weight - 1);
    [~, ck_factor_sort_ind] = sort( rand(size(c_tmp))./ck_factor, 'ascend' );
    % shuffle
    K_tmp_shuffle(ck_factor_sort_ind) = sort( K_tmp, 'descend' );
    K(I == i) =  K_tmp_shuffle;
end




% % output and check results
% in = full(sum(A,1));
% out = full(sum(A,2));
% c_out = corrcoef(out, degree_out);
% c_out = c_out(1,2);
% c_in = corrcoef(in, degree_in);
% c_in = c_in(1,2);
% out_in = corrcoef(in,out);
% out_in = out_in(1,2);
% in_out_corr_comp = [in_out_corrcoef  out_in];
% 
% % Check the distance-dependent rule (exponential)
% dist_avail = lattice_nD_find_dist(Lattice, hw, 1);
% edges = linspace(0,max(dist_avail),50);
% N_avail = histc(dist_avail,edges );
% N_K = histc(K,edges );
% N_exp_check = N_K(:)./N_avail(:);
% % plot(edges, N_exp_check,'o');

% end
end



