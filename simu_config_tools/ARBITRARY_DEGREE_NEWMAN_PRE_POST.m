function [ A ] = ARBITRARY_DEGREE_NEWMAN_PRE_POST( varargin )
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
plot = 0;
N_pre = 500;
N_post = 2000;
order_index_wrt_in_degree = 0;

pj = []; % in-degree
pk = []; % out-degree

max_degree_in = 300;
min_degree_in = 50;
in_degree_vector = min_degree_in:1:max_degree_in; % can be non-integer

gamma_in = 2; % default distirbution is power-law

% Read parameter/distribution setting
for i = 1:length(varargin)/2
    temp = varargin{2*i};
    eval([varargin{2*i-1}, '= temp;']);
end


if isempty(pj) && isempty(pk)
    % in_degree
    
    pj = in_degree_vector.^(-gamma_in);
    pj = pj/sum(pj);
    
    % out_degree
    % scale the expected values of the in/out degree distribution
    % according to N_post/N_pre
    % (Be careful here!!!!!!!!!!)
    pk = pj;
    out_degree_vector = in_degree_vector*(N_post/N_pre);
    
end


% Generate N_post in-degree randomly
num_iteration = 30
while_count = 0;
pair_count = 0;
pairs_in = zeros(N_post,1);
while pair_count < N_post
    randmat = rand(size(pj))/(N_post/num_iteration);
    in_degree_temp = cont2int_rand(in_degree_vector(randmat <= pj));
    pairs_in(pair_count+1:pair_count+length(in_degree_temp)) = in_degree_temp;
    pair_count = pair_count+length(in_degree_temp);
    while_count = while_count+1;
end
while_count
% randomly delete extra pairs
if pair_count > N_post
    delete = randperm(pair_count,pair_count-N_post);
    pairs_in(delete,:) = [];
end


% Generate N_pre out-degree randomly
while_count = 0;
pair_count = 0;
pairs_out = zeros(N_pre,1);
while pair_count < N_pre
    randmat = rand(size(pk))/(N_pre/num_iteration);
    out_degree_temp = cont2int_rand(out_degree_vector(randmat <= pk));
    pairs_out(pair_count+1:pair_count+length(out_degree_temp)) = out_degree_temp;
    pair_count = pair_count+length(out_degree_temp);
    while_count = while_count+1;
end
while_count
% randomly delete extra pairs
if pair_count > N_pre
    delete = randperm(pair_count,pair_count-N_pre);
    pairs_out(delete,:) = [];
end

% Regenerate randomly selected pair until sum(in) == sum(out)
% A lousy technique, cannot believe it works fine
disp('Equalising in and out degree...');
diff = 1;
diff_tol = max([min(in_degree_vector), min(out_degree_vector)]); % If under this tolerance, manually correct it.
while diff ~= 0
    
    if rand > 0.5
        % N_post, pj, pairs_in
        randmat = rand(size(pj))/(N_post/num_iteration);
        in_degree_temp = cont2int_rand(in_degree_vector(randmat <= pj));
        for ii = 1:round(length(in_degree_temp)/2) % re-use randmat
            reproduced = randperm(length(in_degree_temp),1);
            rejected = randperm(N_post,1);
            pairs_in(rejected,:) = in_degree_temp(reproduced);
            in_sum = sum(pairs_in);
            out_sum = sum(pairs_out);
            diff = in_sum-out_sum;
            % Forced correction, may slightly deviate from the desired
            % distribution but increase the speed!
            if  abs(diff) <= diff_tol && diff < 0
                pairs_in(rejected) = pairs_in(rejected)-diff;
                disp('Forced correction done.');
                diff = 0;
                break;
            end
        end
        
    else
        % N_pre, pk, pairs_out
        randmat = rand(size(pk))/(N_pre/num_iteration);
        out_degree_temp = cont2int_rand(out_degree_vector(randmat <= pk));
        for ii = 1:round(length(out_degree_temp)/2) % re-use randmat
            reproduced = randperm(length(out_degree_temp),1);
            rejected = randperm(N_pre,1);
            pairs_out(rejected,:) = out_degree_temp(reproduced);
            in_sum = sum(pairs_in);
            out_sum = sum(pairs_out);
            diff = in_sum-out_sum;
            % Forced correction, may slightly deviate from the desired
            % distribution but increase the speed!
            if  abs(diff) <= diff_tol && diff > 0
                pairs_out(rejected) = pairs_out(rejected)+diff;
                disp('Forced correction done.');
                diff = 0;
                break;
            end
        end
    end
    
end


disp('Equalising done.');
% Randomly connect nodes according to pairs
disp('Generating adjacency matrix...');
stub_in = degree2stub( pairs_in );
stub_out = degree2stub( pairs_out );
A = full(sparse(stub_out,stub_in,1,N_pre,N_post)); % generate adjacency matrix
disp('Connection matrix generated.');
% note that A(i,j) = 1 means connection exists from pre-synaptic neuron i to post-synaptic neuron j

if order_index_wrt_in_degree == 1
    In_degree = sum(A,1);
    % Get new neuron index according to their average firing rate (a bit tricky!!)
    [~,in_degree_ind] = sort(In_degree);
    A = A(:,in_degree_ind);
end

%%%%%%%%%%%%%%%%%%%%%% Check the result %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show the basic graph measurement
min_in_degree_generated = min(pairs_in)
max_in_degree_generated = max(pairs_in)
min_out_degree_generated = min(pairs_out)
max_out_degree_generated = max(pairs_out)
mean_degree_generated = mean(pairs_in)
density_generated = mean_degree_generated /(N_post-1)
end


% Helper function
function [ stub ] = degree2stub( degree )
% This function turn the in/out-degree vector into a stub vector
% The stub vector will be randomly permuted.
% Example:
%       if degree = [3 2 5 1]', for node i = 1, 2, 3, 4
%       the stub vector is
%          stub = [1 1 1 2 2 3 3 3 3 3 4]'; % before permutation
% Note that degree must be column vector!
% stub is also column vector.

Ind = 1:length(degree); % index vector for the nodes
stub = cell2mat(arrayfun(@(x, y) repmat(x, [1 y]), Ind, degree', 'UniformOutput', false));
stub = stub';

% random permutation
edge_num = length(stub);
ind_perm = randperm(edge_num);
stub = stub(ind_perm);

end

function V_out = cont2int_rand(V)
% if V(1) = 2.3, then V_out(1) = 2 with p = 0.7 and V_out(1) = 3 with p = 0.3
%

V_f = floor(V);
p = (V-V_f);
V_r = double(rand(size(V)) < p);
V_out = V_f+V_r;

end