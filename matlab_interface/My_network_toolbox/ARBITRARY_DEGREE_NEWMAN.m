function [ A ] = ARBITRARY_DEGREE_NEWMAN( varargin )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

addpath(genpath('/import/yossarian1/yifan/Project1/My_network_toolbox'));


% Generate random directed graphs with given joint in-and-out
% degree distribution p(j,k)
% If p(j,k) is unknown, assume p(j) and p(k) is uncorrelated, i.e.,
% p(i,j) = p(j)*p(k), which usually is not true.
% Ref. Newmann et al., 2000, Random graphs with arbitrary degree
% distribution and their applications

% Default parameters
order_index_wrt_in_degree = 0;
plot = 0;
N = 2000;
directed = 1; % directed graph or undirected
pj = []; % in-degree
pk = []; % out-degree
pjk = [];
max_degree = 300;
min_degree = 50;
gamma = 2; % default distirbution is power-law

% Read parameter/distribution setting
for i = 1:length(varargin)/2
    temp = varargin{2*i};
    eval([varargin{2*i-1}, '= temp;']);
end

% Get joint degree distribution
% In and out degree distribution example
% Note that the in/out degree distribution vector pj/pk should be
% distribution_function(1:1:max_degee)!!!!

fprintf('The in/out degree distribution vector pj/pk should be distribution_function(1:1:max_degee)!\n');
if isempty(pj) && isempty(pk) && isempty(pjk)
    pj = [zeros(1,min_degree-1),(min_degree:max_degree).^(-gamma)];
    pj = pj/sum(pj);
    pk = pj;
end
% If joint distribution not specifid, assume non-correlation and
% generate it from simply multiplying pj and pk
if isempty(pjk)
    [pj_mesh,pk_mesh] = meshgrid(pj,pk);
    pjk = pj_mesh.*pk_mesh;
end
N, directed
if plot == 1
    figure('NumberTitle', 'Off', 'Name', 'Degree distribution')
    imagesc(pjk);
end

% Generate N in-out degree pairs i(j,k) for each point randomly
num_iteration = 30
while_count = 0;
pair_count = 0;
pairs = zeros(N,2);
while pair_count < N
    randmat = rand(size(pjk))/(N/num_iteration);
    [ind_i, ind_j] = find(randmat <= pjk);
    pairs(pair_count+1:pair_count+length(ind_i),:) = [ind_i ind_j];
    pair_count = pair_count+length(ind_i);
    while_count = while_count+1;
end
while_count

% randomly delete extra pairs
if pair_count > N
    delete = randperm(pair_count,pair_count-N);
    pairs(delete,:) = [];
end

% Regenerate the last pair until sum(in) == sum(out)
% A lousy technique, cannot believe it works fine
disp('Equalising in and out degree...');
diff = 1;
diff_tol = max(20,min_degree); % If under this tolerance, manually correct it.
while diff ~= 0
    randmat = rand(size(pjk))/(N/num_iteration);
    [ind_i, ind_j] = find(randmat<=pjk);
    for ii = 1:round(length(ind_i)/2) % re-use randmat
        reproduced = randperm(length(ind_i),1);
        rejected = randperm(N,1);
        pairs(rejected,:) = [ind_i(reproduced) ind_j(reproduced)];
        in_sum = sum(pairs(:,1));
        out_sum = sum(pairs(:,2));
        diff = in_sum-out_sum;
        % Forced correction, may slightly deviate from the desired
        % distribution but increase the speed!
        if  abs(diff) <= diff_tol
            if diff > 0
                pairs(rejected,2) = pairs(rejected,2)+diff;
            else
                pairs(rejected,1) = pairs(rejected,1)-diff;
            end
            disp('Forced correction done.');
            diff = 0;
            break;
        end
    end
end
disp('Equalising done.');
% Randomly connect nodes according to pairs
disp('Generating adjacency matrix...');
no_self = 1;
[A] = generate_graph_given_in_out_degree(pairs(:,1), pairs(:,2), no_self);
disp('Adjacency matrix generated.');
% note that A(i,j) = 1 means connection exists from pre-synaptic neuron i to post-synaptic neuron j
if order_index_wrt_in_degree == 1
    In_degree = sum(A,1);
    % Get new neuron index according to their average firing rate (a bit tricky!!)
    [~,in_degree_ind] = sort(In_degree);
    A = A(:,in_degree_ind);
end

%%%%%%%%%%%%%%%%%%%%%% Check the result %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show the basic graph measurement
min_in_degree_generated = min(pairs(:,1))
max_in_degree_generated = max(pairs(:,1))
min_out_degree_generated = min(pairs(:,2))
max_out_degree_generated = max(pairs(:,2))
mean_degree_generated = mean(pairs(:,1))
density_generated = mean_degree_generated /(N-1)
% Degree
[In_degree,Out_degree] = degrees_dir(A); % (1-by-N) vectors
[Joint_degree,~,~,~] = jdegree(A);
% Joint_degree is a matrix in which the value of each element (u,v)
% corresponds to the number of nodes that have u-1 outgoing connections
% and v-1 incoming connections.



end


