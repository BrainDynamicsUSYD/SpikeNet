% construct spatially embedded network with given in and out-degree and distance
% distributions

%function [ll_matrix, N, Lattice] = spatial_embed_in_out_network(hw, P0_init, degree_CV, tau_d)
% P0: average connection probability
% The final P0 will be different from P0_init
% degree_CV = degree_dist_std/degree_dist_mean
%
% tau_d: decay constant for the exponential distribution (units: neurons)
%
% default: log-normal digree distribution and exponential distance wiring
% rule


hw = 21; % 31
P0_init = 0.1;
degree_CV = 0.1;
tau_d = 10; %
in_out_corrcoef = 0.3;


D = 2; % Pulin: 2D is enough
% hw = 31; % (31*2+1)^2 = 3969 ~ 4000

[Lattice, N] = lattice_nD(D, hw);
% fprintf('The number of node is %d. \n' , N);

% link-length matrix (sparse)
I = []; % pre node index
J = []; % post node index
K = []; % link length

% % degree distribution
% log-normal
deg_size = N;
deg_mean = N*P0_init;
deg_std = degree_CV*deg_mean;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_tol = deg_mean/2;
diff = deg_mean;
while diff > diff_tol
    [ degree_in_out ] = ceil(my_logn_rand( [deg_mean deg_mean], [deg_std deg_std], ...
        [1 in_out_corrcoef; in_out_corrcoef, 1], N, 'mu_sigma_log' ));
    diff = abs(sum(degree_in_out(1:end-1,1)) - sum(degree_in_out(1:end-1,2)));
end
% Regenerate the last pair to best match sum(in) == sum(out)
% If there is still any mismatch, manually correct it.
% A lousy technique, cannot believe it works fine
[ degree_in_out_re ] = ceil(my_logn_rand( [deg_mean deg_mean], [deg_std deg_std], ...
    [1 in_out_corrcoef; in_out_corrcoef, 1], N, 'mu_sigma_log' ));
diff = sum(degree_in_out(1:end-1,1)) - sum(degree_in_out(1:end-1,2));
corr = degree_in_out_re(:,2) - degree_in_out_re(:,1);
[resi, corr_best_ind] = min( abs(diff - corr) );
last_pair = degree_in_out_re(corr_best_ind,:);
corr_best = last_pair(2) - last_pair(1);
if resi > 0 && diff > corr_best
    last_pair(2) =  last_pair(2) + (diff - corr_best);
    degree_in_out(end,:) =  last_pair;
elseif resi > 0 && diff < corr_best
    last_pair(1) =  last_pair(1) + (corr_best - diff);
    degree_in_out(end,:) =  last_pair;
else
    degree_in_out(end,:) =  last_pair;
    if sum(degree_in_out(:,1) - degree_in_out(:,2)) ~= 0
        degree_in_out(end,:) =  flipud(last_pair);
    end
end


degree_in = ceil(degree_in_out(:,1));
degree_out = ceil(degree_in_out(:,2));


%%%%%%%%%%%%%%%%%%%
% The following solution ensures that degree_in and degree_out are strictly
% used.


degree_in_left = degree_in;
progress = 0;
for i = randperm(N)
    post_dist = lattice_nD_find_dist(Lattice, hw, i);
    dist_factor = exp(-post_dist/tau_d); % exponential distribution, will be properly scaled later;
    dist_factor(i) = 0; % no self-connection
    
    joint_factor = dist_factor.*degree_in_left; % will be properly scaled later;
    
    % establish connections
    [~, ind] = sort( rand(N,1)./joint_factor, 'ascend' );
    chosen_j = ind(1: degree_out(i) );
    
    degree_in_left(chosen_j) = degree_in_left(chosen_j) - 1;
    degree_in_left(degree_in_left == 0) = NaN;
    
    % store results
    chosen_j = chosen_j(:); % column vector
    I = [I; i*ones(size(chosen_j)) ]; %#ok<AGROW>
    J = [J; chosen_j]; %#ok<AGROW>
    K_tmp = post_dist(chosen_j);
    K = [K; K_tmp(:)]; %#ok<AGROW>
    
    progress = progress + 1;
    progress/N %#ok<NOPTS>
end
clc;

% %%%%%%%%%%%%%%%%
% degree_in_prob = degree_in; % will be properly scaled later;
% bad_prob = 0;
%
% for i = 1:N % randomize it
%
%     post_dist = lattice_nD_find_dist(Lattice, hw, i);
%     dist_prob = exp(-post_dist/tau_d); % exponential distribution, will be properly scaled later;
%     dist_prob(i) = 0; % no self-connection
%
%     joint_prob = dist_prob.*degree_in_prob; % will be properly scaled later;
%
%     % establish connections
%     joint_prob_scaled =  joint_prob/sum(joint_prob) * degree_out(i);
%     chosen_j = find( rand(N, 1) < joint_prob_scaled );
%
%     % check bad prob
%     bad_prob = bad_prob + sum(joint_prob_scaled  > 1);
%
%
%     % store results
%     chosen_j = chosen_j(:); % column vector
%     I = [I; i*ones(size(chosen_j)) ]; %#ok<AGROW>
%     J = [J; chosen_j]; %#ok<AGROW>
%     K_tmp = post_dist(chosen_j);
%     K = [K; K_tmp(:)]; %#ok<AGROW>
%
%     i/N
% end
% %%%%%%%%%%%%%%%%
% bad_prob

clc;


% output and check results
ll_matrix = sparse(J,I,K);
A = ll_matrix > 0;

in = full(sum(A,2));
out = full(sum(A,1));
c_out = corrcoef(out, degree_out);
c_out = c_out(1,2)
c_in = corrcoef(in, degree_in);
c_in = c_in(1,2)

out_in = corrcoef(in,out);
out_in = out_in(1,2);

in_out_corr_comp = [in_out_corrcoef in_out_corr out_in]
mean_comp = [deg_mean mean(in) mean(out)]
std_comp = [deg_std std(in) std(out)]

% Check the distance-dependent rule (exponential)
edges = linspace(0,max(dist_avail),50);
dist_avail = lattice_nD_find_dist(Lattice, hw, 1);
N_avail = histc(dist_avail,edges );
N_K = histc(K,edges );
N_exp_check = N_K./N_avail;
plot(edges, N_exp_check,'o');




