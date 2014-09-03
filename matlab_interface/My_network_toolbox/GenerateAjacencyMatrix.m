function [ A, C_i, C_j, K_ij ] = GenerateAjacencyMatrix( choice, varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
  
addpath(genpath('C:\Users\yigu8115\Desktop\My network toolbox'));

switch lower(choice)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'homo_gaussian'
        % Default Gaussian coupling and scale parameters
        a = 40; % (a-by-a) neurons
        sigma = 5;
        J = 1;
        
        % Read non-default inputs
        for i = 1:length(varargin)/2
            temp = varargin{2*i};
            eval([varargin{2*i-1}, '= temp;']);
        end
        sigma
        J
        
        % Set up spatially extented (2D square) neuron sheet
        N = a*a % total number of neurons
        radius = 3*sigma; % maximum length of connection
        [ C_matrix , D_vector ] = CircularConnection( a, radius, 'periodic' );
        K_vector = J*exp(-0.5*D_vector.^2/sigma^2);
        K_matrix = ones(N,1)*K_vector;
        
        % Change data format
        % Re-format the data into vectors C_i, C_j, K_ij
        [C_i,~,C_j] = find(C_matrix.*(K_matrix>0));
        [~,~,K_ij] = find(K_matrix);
        
        % Re-format the data into A
        A_sparse = sparse(C_i,C_j,K_ij,N,N); % sparse adjacency matrix
        clear C_matrix K_matrix
        A = full(A_sparse);
        clear A_sparse
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'square_exponential'
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
        radius
        density
        lambda
        
        N = a*a % total number of neurons  
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
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case  'gaussian_embedding_exponential'
        %Default parameters
        num_control = 0.1; % parameter controls total number of nodes
        sigma = 0.3; % X and Y coordinate range is [0,1]
        length_decay = 2; % length constant for EDR
        num_iteration = 50; % control the graph generation speed
        density = 0.01;
        
        % Read non-default inputs
        for i = 1:length(varargin)/2
            temp = varargin{2*i};
            eval([varargin{2*i-1}, '= temp;']);
        end
        num_control
        sigma
        length_decay
        num_iteration
        density
        
        
        %-------------------------------------------------------%
        x1 = -1:.01:1; x2 = -1:.01:1;
        [X,Y] = meshgrid(x1,x2);
        X = X(:);
        Y = Y(:);
        F = mvnpdf([X Y],[0 0], [sigma 0;0 sigma]);
        randMat = rand(length(X),1)/num_control;
        cut = find(randMat > F);
        X(cut) = [];
        Y(cut) = [];
        % Inter-Node Distance Distribution
        % figure('Name','2D Node distribution','NumberTitle','off');plot(X,Y,'.');
        N = length(X) % total number of nodes
        DistMat = zeros(N,N);
        for i = 1:N
            for j = i:N
                DistMat(i,j) = sqrt((X(i)-X(j))^2+(Y(i)-Y(j))^2);
                DistMat(j,i) = DistMat(i,j);
            end
        end
        figure('Name','Inter-node distance distribution','NumberTitle','off');hist( DistMat(:),50);
        % Apply Exponential Distance Rule (EDR)
        %-------------------------------------------------------%
        pMat = exp(-DistMat./length_decay);
        pMat = pMat/(sum(sum(pMat))); % normalisation
        A = zeros(N,N); % adjacency matrix
        count = 0;
        while_count = 0;
        tot_syn = N^2*density;
        while count < tot_syn
            randMat = rand(N,N)*num_iteration/tot_syn;
            hit = randMat<pMat;
            A(hit) = A(hit)+1;
            count = sum(sum(A));
            while_count = while_count+1;
        end
        while_count
        A = A/sum(sum(A));
        [C_i,C_j,K_ij] = find(A);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    case 'macaque71'
        load('macaque71.mat');
        A = CIJ;
        N = size(A,1)
        [C_i,C_j,K_ij] = find(A);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    case 'coactivation_matrix'
        load('Coactivation_matrix.mat');
        A = Coactivation_matrix;
        N = size(A,1)
        [C_i,C_j,K_ij] = find(A);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'cat_all'
        load('cat.mat');
        A = CIJall;
        N = size(A,1)
        [C_i,C_j,K_ij] = find(A);
end

