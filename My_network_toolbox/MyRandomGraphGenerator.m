function [ A ] = MyRandomGraphGenerator( Method, varargin )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

addpath(genpath('/import/yossarian1/yifan/Project1/My_network_toolbox'));

switch upper(Method)
    %������������������������������������������%
    case 'E_R'
        % Erdos-Renyi random graph (directed or undirected)
        % All possible pairs of nodes are connected with probability 'p'
        % Features:
        %    a) degree distribution is a Poisson distribution, i.e., every node has
        %    about the same degree
        %    b) low clustering coefficient
        %    c) no hubs
        %    d) see http://en.wikipedia.org/wiki/Giant_component
        
        % Default parameters
        N = 1000;
        p = 0.3;
        directed = 1; % directed graph or undirected
        
        
        % Read parameter setting
        for i = 1:length(varargin)/2
            temp = varargin{2*i};
            eval([varargin{2*i-1}, '= temp;']);
        end
        
        fprintf('E-R Graph: N=%d, p=%g, directed=%d\n', N, p, directed);
        
        % Generate graph
        A = zeros(N,N);
        A(rand(N) <= p) = 1; % directed graph (not symmetric)
        A(logical(eye(N))) = 0; % set diagonal entries to zero, no self-connection
        if directed == 0 % symmetrise it
            A = triu(A,1); % extract upper triangular part
            A = A + A';
        end
        
    %������������������������������������������%
    case 'E_R_PRE_POST'
        % Erdos-Renyi model extended to connection between two population
        % (from pre to post-synaptic population)
        
        % Default parameters
        N_pre = 1000;
        N_post = 1000;
        p = 0.3;

        % Read parameter setting
        for i = 1:length(varargin)/2
            temp = varargin{2*i};
            eval([varargin{2*i-1}, '= temp;']);
        end
        fprintf('E-R Graph pre-post: N_pre=%d, N_post=%d, p=%g\n', N_pre, N_post, p);
        
        % Generate graph
        A = zeros(N_pre,N_post);
        A(rand(N_pre, N_post) <= p) = 1;
        %disp('Note that A is connection from pre to post-synaptic population!\n');
        %disp('A is thus not square!');
    %������������������������������������������%    
    case 'MODULAR_E_R'
        % Build a random modular graph (with equal module size!)
        % INPUTs: number of nodes, number of modules, total link density,
        %         and proportion of links within modules compared to links across
        % OUTPUTs: adjacency matrix, modules to which the nodes are assigned
        
        % N - number of nodes
        % c - number of clusters/modules
        % p - overall probability of attachment
        % r - proportion of links within modules
        
        % Default parameters
        N = 1000;
        m = 4;
        p = 0.1;
        r = 0.8;
        directed = 1; % directed graph or undirected
        
        % Read parameter setting
        for i = 1:length(varargin)/2
            temp = varargin{2*i};
            eval([varargin{2*i-1}, '= temp;']);
        end
        fprintf('Modular E-R Graph: N=%d, m=%d, p=%g, r=%g, directed=%d\n', N, m, p, r, directed);
        
        % Generate graph
        z = p*(N-1);   % Erdos-Renyi average degree
        module_size = N/m; % average size of each module
        p_within = r*z/(module_size-1); % probability of attachment within module = allowed connections within module/total possible connections within node
        p_between = (1-r)*z/(N-module_size); % probability of attachment across modules

        modules = ceil((1:N)*m/N); % node module assignment vector
        [Module_i, Module_j] = meshgrid(modules,modules);
        randmat = rand(N);
        A = zeros(N); % initialize adjacency matrix
        A((randmat <= p_within) & (Module_i == Module_j)) = 1;
        A((randmat <= p_between) & (Module_i ~= Module_j)) = 1;
        
        if directed == 0 % symmetrise it
            A = triu(A,1); % extract upper triangular part
            A = A + A';
        end
        
        
    %������������������������������������������%
    case 'SQUARE_EXPONENTIAL'
        % Default parameters
        a = 50; % (a-by-a) neurons
        radius = 20; % maximum length of connection
        density = 0.01; % connection density
        lambda = 10; % distance constant
        Num_iteration = 50; % number of iteration for network generation
        
        % Read non-default inputs
        for i = 1:length(varargin)/2
            temp = varargin{2*i};
            eval([varargin{2*i-1}, '= temp;']);
        end
        N = a*a; % total number of neurons
        fprintf('SQUARE_EXPONENTIAL graph: N=%d, radius=%d, denisty=%g,lambda=%d\n', N, radius, density, lambda);
        
        [ C_matrix , D_vector ] = CircularConnection( a, radius, 'periodic' ); % Set up spatially extented (2D square) neuron sheet
        [ K_matrix ] = EDR_Connection( density, lambda, Num_iteration, C_matrix, D_vector );% EDR_Connection
        
        % Change data format
        % Re-format the data into vectors C_i, C_j, K_ij
        [C_i,~,C_j] = find(C_matrix.*(K_matrix>0));
        [~,~,K_ij] = find(K_matrix);
        
        % Re-format the data into A
        A_sparse = sparse(C_i,C_j,K_ij,N,N); % sparse adjacency matrix
        clear C_matrix K_matrix C_i C_j K_ij
        A = full(A_sparse);
        clear A_sparse

end
        
end