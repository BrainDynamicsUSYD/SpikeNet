% construct spatially embedded network with given degree and distance
% distributions

% ref
% Structural properties of spatially embedded networks, Kosmidis, Kosmas et
% al, 2008, Europhysics Letters


function [ll_matrix, N, Lattice] = spatial_embed_network(hw, P0_init, degree_CV, tau_d)
% P0: average connection probability 
% The final P0 will be different from P0_init
% degree_CV = degree_dist_std/degree_dist_mean
%
% tau_d: decay constant for the exponential distribution (units: neurons)
%
% default: log-normal digree distribution and exponential distance wiring
% rule

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
deg_range = 1:1:N;
deg_size = N;
deg_mean = N*P0_init;
deg_std = degree_CV*deg_mean;
degree = pmf_yg('logn', [deg_mean deg_std], deg_range, deg_size);
% hist(degree, 1000);set(gca,'xscale','log') %, 'yscale','log')
% mean(degree)
% std(degree)

for i = 1:N
    post_dist = lattice_nD_find_dist(Lattice, hw, i);
    post_prob_0 = exp(-post_dist/tau_d); % exponential distribution
    post_prob = post_prob_0/sum(post_prob_0)*degree(i);
    post_prob(i) = 0; % no self-connection
    
    % establish connections
    chosen_j = find( rand(N, 1) < post_prob );
    % store results
    chosen_j = chosen_j(:); % column vector
    I = [I; i*ones(size(chosen_j)) ];
    J = [J; chosen_j];
    K_tmp = post_dist(chosen_j);
    K = [K; K_tmp(:)];
end

% % out-degree log-normal:
% ll_matrix = sparse(I,J,K);

% in-degree log-normal:
ll_matrix = sparse(J,I,K);


end






