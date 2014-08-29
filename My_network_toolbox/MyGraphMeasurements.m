function [ GraphMeasure ] = MyGraphMeasurements( A, varargin )
%Measurements for directed and weighted graph (network), especially chosen for
%neuron networks.
%   Detailed explanation goes here

% Note: many of the measurements are not directly applicable to directed
% and weighted network!!
% Even some measurements for directed and weighted network are not defined
% for neuron networks! (transportation, communication, social...)
% Neuron networks are more like transportation networks..

disp('MyGraphMeasurements...');

% default measurement choices
get_clustering_coef = 1;
get_strength = 1;
get_path_length = 0;
get_associativity = 0;
get_modularity = 1;
get_rich_club_coef = 1;
get_degree_and_density = 1;

% Read parameter setting (measurement choices)
for i = 1:length(varargin)/2
    temp = varargin{2*i};
    eval([varargin{2*i-1}, '= temp;']);
end

% % Sampling???
% sample_size = 1000; % if node number is larger than it, sample the connection matrix!
% node_number = length(A(1,:));
% if  node_number > sample_size
%     sample_ind = sort(randperm(node_number, sample_size));
%     A_sample = A(sample_ind,sample_ind);
% else
%     A_sample = A;
% end



%--------------------------------------------------------------------------------------------%
%--------------------------------------------------------------------------------------------%
% (Both weighted and directed!)


% Clustering coefficients
if get_clustering_coef == 1;
    tic;
    disp('\t Calculating clustering coefficients...')
    [ GraphMeasure.clustering_coef ] = clustering_coef_wd(A); % (N-by-1) vector, Fagiolo (2007) as well
    % figure('NumberTitle','Off','Name','Clustering coefficient')
    % Calculate C(K) ??
    toc;
end

% Strength (weighted degree)
if get_strength == 1;
    tic;
    disp('\t Calculating node strength...')
    [GraphMeasure.in_strength, GraphMeasure.out_strength] = strengths_dir(A);
    % However, note that strength (total connection weight) and degree (total connection number) should be both examined!
    GraphMeasure.strength_total = sum(GraphMeasure.in_strength); % also = sum(instrength), a measure of global connection strength
    toc;
end

% Average path length
if get_path_length == 1;
    tic;
    disp('\t Calculating path length...')
    L = 1./A; % L: mapping from weight to distance (usually weight inversion)
    [GraphMeasure.distance_weighted_path, GraphMeasure.distance_number_of_edges] = distance_wei(L);
    % D,      distance (shortest weighted path) matrix
    % B,      number of edges in shortest weighted path matrix
    [GraphMeasure.char_path_length,GraphMeasure.efficiency,GraphMeasure.ecc,GraphMeasure.radius,GraphMeasure.diameter] = charpath(GraphMeasure.distance_weighted_path);
    toc;
end


% Associativity
%   The assortativity coefficient is a correlation coefficient between the
%   strengths (weighted degrees) of all nodes on two opposite ends of a link.
%   A positive assortativity coefficient indicates that nodes tend to link to
%   other nodes with the same or similar strength.
if get_associativity == 1;
    tic;
    disp('\t Calculating associativity...');
    GraphMeasure.associativity_out_in = assortativity_wei(A,1);%   1, directed graph: out-strength/in-strength correlation
    GraphMeasure.associativity_in_out = assortativity_wei(A,2);%   2, directed graph: in-strength/out-strength correlation
    GraphMeasure.associativity_out_out = assortativity_wei(A,3);%   3, directed graph: out-strength/out-strength correlation
    GraphMeasure.associativity_in_in = assortativity_wei(A,4);%   4, directed graph: in-strength/in-strength correlation
    toc;
end

% Modularity
if get_modularity == 1;
    tic;
    disp('\t Calculating modularity...');
    gamma = 1;
    [~, GraphMeasure.modularity_louvain] = modularity_louvain_dir(A,gamma); % seem to deal with direction more explicitly
    [~, GraphMeasure.modularity_newman] = modularity_newman_dir(A,gamma);
    toc;
end











%--------------------------------------------------------------------------------------------%
%--------------------------------------------------------------------------------------------%
% Weighted + Total degree (no distinction between in/out degree)


% Rich-club coefficient
if get_rich_club_coef == 1;
    tic;
    disp('\t Calculating rich-club coefficient...');
    [ GraphMeasure.rich_club_coef ] = rich_club_wd(A); % vector
    % The rich club coefficient at level k is the fraction of edges that
    % connect nodes of degree k (indegree + outdegree) or higher out of the maximum number of edges
    % that such nodes might share.
    toc;
end











%--------------------------------------------------------------------------------------------%
%--------------------------------------------------------------------------------------------%
% Un-weighted + directed
% Not considering weights! Only topological meaningful!
if get_degree_and_density == 1;
    tic;
    disp('\t Calculating node degree and density...');
    % Degree (Connection weights are ignored)
    [GraphMeasure.in_degree,GraphMeasure.out_degree] = degrees_dir(A); % (1-by-N) vectors
    [GraphMeasure.joint_degree,~,~,~] = jdegree(A);
    % Joint_degree is a matrix in which the value of each element (u,v)
    % corresponds to the number of nodes that have u-1 outgoing connections
    % and v-1 incoming connections.
    [GraphMeasure.density,~,~] = density_dir(A); % Density is the fraction of present connections to possible connections (
    toc;
end


disp('MyGraphMeasurements done.');

end

