% Generate an instance of a random graph
% with given properties accirding to
% graph type:
% ER (Erdos-Renyi) param is link
% probability
% LATTICE (nearest neighbour) param
% is dimension
% SO (spread-out) param is distance L of
% furthest neighbours (dimension is 2)
% AB (Albert-Barabasi) param is node
% degree
% WS (Watts-Strogatz) param is
% probability of longrange edge (base
% graph is 2d)
function [E]=generate_random_graph(graph_type, N, param)
graph_ER = 0; graph_LATTICE=1; graph_BA=2; graph_SO=3; graph_WS=4; graph_DEGSEQ=5; graph_REGULAR=6;
TOL=0.0000000000000000001;


switch graph_type
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Erdos-Renyi random graph
    case graph_ER
        sprintf('Generating Erdos-Renyi random graph\n')
        p = param; % parameter for ER graphs is link probability
        q = sqrt(p); % Compensate for randomizing twice

        E=rand(N) < q;
        E=E.*E';
        for i=1:N
            E(i,i)=0;
        end
        E = sparse(E); % turn into sparse

        return;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Z^d grid
    case graph_LATTICE
        sprintf('Generating nearest neighbour lattice graph\n')
        d=floor(param);
        if (~isreal(d) || (d ~= param) || (d>50))
            sprintf('Error! Dimension must be positive integer smaller than 50\n')
            E=[];
            return;
        end
        if (  ((round(N^(1/d)))  - (N^(1/d)))^2 > TOL  )
            sprintf('Error! For a lattice Z^d N must be an integer d-power\n')
            E=[];
            return;
        end

        % Set edge length
        L = round(N^(1/d));

        % Start filling the lattice. Use toric boundary conditions, so each
        % vertex has 2*d neighbors.
        E_num = N*d; % number of edges

        ii = repmat([0:N-1], 1, d);
        J = zeros(1,N*d);

        for i=1:d
            target = [0:N-1] + L^(i-1); % the target without modulu
            bad_ind = find(mod(floor(target/L^(i-1)) , L)==0);
            target(bad_ind) = target(bad_ind) - L^i;
            J((i-1)*N+1:i*N) = target;
        end

        % Transfer I,J vectors to a sparse matrix
        E = sparse(ii+1,J+1,ones(1,E_num),N,N);
        if(L > 2) % The hyper-cube is a seperate case
            E = E+E'; % Make symmetric and undirected
        end
        return;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Barabasi-Alberts preferential-attachment network
    case  graph_BA
        sprintf('Generating Barabasi-Albert random graph\n')
        m = param; % parameter is out degree of new nodes

        gamma = 1; % This should be changed and given as an input argument!!!

        if ( (~isreal(m)) || (m ~= param) || (m < 1) || (floor(m) ~= m) )
            sprintf('Error! Degree must be positive integer\n')
            E=[];
            return;
        end
        if (N < m+1)
            sprintf('Error! For m-out-degree Albert-Barabasi model you need at least m+1 nodes\n')
            E=[];
            return;
        end

        rowmat = zeros(1,N*m); % This is the smaller index
        colmat = zeros(1,N*m); % This is the larger index
        deg_vec = zeros(1,N); % The degree vector

        deg_vec(1:m) = m-1; % We start with the clique C_m

        cur_edge=1;     % The current edge
        for i=1:m-1
            rowmat(cur_edge:cur_edge+m-i-1)=i;
            colmat(cur_edge:cur_edge+m-i-1)=i+1:m;
            cur_edge = cur_edge + m-i;
        end



        % Enumerate the vertices we add to the network
        for i=m+1:N

            % Pick at random who to connect to
            p_cum_vec = cumsum( deg_vec(1:i-1).^gamma ./ sum(deg_vec(1:i-1).^gamma) );

            m_chosen_vertices = zeros(1,m); % Vector of vertices to link to our vector

            r = rand(1,m);
            for edge=1:m
                m_chosen_vertices(edge) = min(find(r(edge) < p_cum_vec));
            end
            m_chosen_vertices = unique(m_chosen_vertices);

            rowmat(cur_edge:cur_edge+length(m_chosen_vertices)-1) = i;
            colmat(cur_edge:cur_edge+length(m_chosen_vertices)-1) = m_chosen_vertices;
            cur_edge = cur_edge + length(m_chosen_vertices);

            deg_vec(m_chosen_vertices) = deg_vec(m_chosen_vertices)+1;
            deg_vec(i) = length(m_chosen_vertices);

        end

        rowmat = rowmat(1:cur_edge-1); colmat = colmat(1:cur_edge-1);


        % Note: The average degree is a little less than m because of repeating edges
        E  = sparse(rowmat,colmat,ones(1,cur_edge-1),N,N);

        E = E+E'; % make symmetric
        return;

        % % % % % % Danny's buggy code:
        % % % % % %%%%%%%%%%%%%%%%%%%%
        % % % % %     rowmat = zeros(1,N*m);
        % % % % %     colmat = zeros(1,N*m);
        % % % % %
        % % % % %     ii = 1;
        % % % % %     % start with graph K_(m+1)
        % % % % %     for k=1:m+1
        % % % % %         for j=1:m+1
        % % % % %             if (k ~= j)
        % % % % %                 rowmat(ii) = k;
        % % % % %                 colmat(ii) = j;
        % % % % %                 ii = ii+1;
        % % % % %             end
        % % % % %         end
        % % % % %     end
        % % % % %     degseq = [zeros(1,m+1)+m,zeros(1,N-(m+1))];
        % % % % %     degtot = m*(m+1);
        % % % % %
        % % % % %     % attach rest of the N-(m+1) nodes
        % % % % %     for srcnode=m+2:N
        % % % % %         % attach new node
        % % % % %         for edge=1:m
        % % % % %             % choose edge according to degree
        % % % % %             % sequence
        % % % % %             rowmat(ii) = srcnode;
        % % % % %             dstnode = srcnode; % hack to start the next while
        % % % % %             while (dstnode == srcnode)
        % % % % %                 r = randint(1,1,degtot);
        % % % % %                 dstnode = 1;
        % % % % %                 while (r>degseq(dstnode))
        % % % % %                     r = r - degseq(dstnode);
        % % % % %                     dstnode = dstnode+1;
        % % % % %                 end
        % % % % %             end
        % % % % %             colmat(ii) = dstnode;
        % % % % %             ii=ii+1;
        % % % % %             degseq(srcnode) = degseq(srcnode)+1;
        % % % % %             degseq(dstnode) = degseq(dstnode)+1;
        % % % % %             degtot = degtot+1;
        % % % % %         end
        % % % % %     end
        % % % % %
        % % % % %     % XXX the generated degrees might be
        % % % % %     % smaller than the specified param
        % % % % %     % because of randomly repeating edges.
        % % % % %     % this should be fixed.
        % % % % %     E  = sparse(rowmat,colmat,ones(1,N*m),N,N,'unique');
        % % % % %     return;
        % % % % %     % End of Danny's buggy code
        % % % % %     %%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate Watts-Strogatz graph, which is a square grid with some random
        % edges ???
    case graph_WS
        sprintf('Generating Watts-Strogatz random graph\n')
        p = param;
        if ((round(N^(1/2))) ~= (N^(1/2)))
            sprintf('Error! For a lattice Z^2 N must be a square\n')
            E=[];
            return;
        end

        % Set edge length
        R = N^(1/2);
        if ((p<0) || (p>=1.0))
            sprintf('Error! Probability mush be in [0,1)\n')
            E=[];
            return;
        end

        % fill plane. and add longrange edges
        max_ii = 4*N + 10*p*N; % an upper bound for number of edges

        rowmat = ones(1,max_ii);
        colmat = ones(1,max_ii);
        ii = 1;
        for m=1:R
            for n=1:R
                if (ii>max_ii-6)
                    sprintf('Error! more than expected edges\n')
                    E=[];
                    return;
                end
                % add the 4 lattice edges
                rowmat(ii) = (m-1)*R+n-1+1;
                colmat(ii) = mod((m-1+R-1),R)*R+mod(n+0+R-1,R)+1;
                ii = ii+1;
                rowmat(ii) = (m-1)*R+n-1+1;
                colmat(ii) = mod((m+1+R-1),R)*R+mod(n+0+R-1,R)+1;
                ii = ii+1;
                rowmat(ii) = (m-1)*R+n-1+1;
                colmat(ii) = mod((m+0+R-1),R)*R+mod(n-1+R-1,R)+1;
                ii = ii+1;
                rowmat(ii) = (m-1)*R+n-1+1;
                colmat(ii) = mod((m+0+R-1),R)*R+mod(n+1+R-1,R)+1;
                ii = ii+1;
                % add longrange edge with probablity p
                if (rand()<p)
                    rowmat(ii) = (m-1)*R+n-1+1;
                    colmat(ii) = mod((m+randint(1,1,R)+R-1),R)*R+mod(n+randint(1,1,R)+R-1,R)+1;
                    while (rowmat(ii) == colmat(ii))
                        colmat(ii) = ...
                            mod((m+randint(1,1,R)+R-1),R)*R+mod(n+randint(1,1,R)+R-1,R)+1;
                    end
                    ii = ii+1;
                end
            end
        end

        E  = sparse(rowmat,colmat,[ones(1,ii),zeros(1,max_ii-ii)],N,N,'unique');
        return;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate a spread-out lattice graph - what the hell is this?
    case graph_type == graph_SO
        sprintf('Generating spread-out lattice graph\n')
        L = param
        if ((round(N^(1/2))) ~= (N^(1/2)))
            sprintf('Error! For a lattice Z^2 N must be a square\n')
            E=[];
            return;
        end

        % Set edge length
        R = N^(1/2)
        if (2*L+1 >= R)
            sprintf('Error! Neighbour distance must be smaller than half grid edge\n')
            E=[];
            return;


        end

        % fill plance. use toric boundary
        % conditions, so each vertex has
        % (2L+1)^2-1 neighbours.
        E_num = (2*L+1)^2-1 % number of edges

        N*E_num
        rowmat = zeros(1,N*E_num);
        colmat = zeros(1,N*E_num);
        ii = 1;
        for m=1:R
            for n=1:R
                for x=-L:1:L
                    for y=-L:1:L
                        if ((x~=0) || (y~=0))
                            rowmat(ii) = (m-1)*R+n-1+1;
                            colmat(ii) = mod((m+x+R-1),R)*R+mod(n+y+R-1,R)+1;
                            ii = ii+1;
                        end
                    end
                end
            end
        end

        E  = sparse(rowmat,colmat,ones(1,N*E_num),N,N);
        return;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate a graph with a given degree sequence
    case graph_DEGSEQ
        % Start by generating one specific graph with the given degree
        % sequence
        deg_seq = sort(param, 'descend'); % Sort the degrees
        % Check if this is a graphic sequence, i.e. a valid degree sequence
        % for a graph (Erdos and Galai, 1960)
        cumsum_degs = cumsum(deg_seq);
        for r=1:N-1
            if(cumsum_degs(r) > r*(r-1) + sum(min(r, deg_seq(r+1:end))))
                sprintf('Error! Not a valid degree sequence!!!')
                E=[];
                return;
            end
        end
        E = speye(N)-speye(N); % Start with all zero sparse matrix
        Indexes = [1:N]+1; Index_i=1;
        for i=1:N
            [max_deg max_ind] = max(deg_seq);
            [sorted_degs sorted_ind] = sort(deg_seq, 'descend');
            if(max_deg > 0)
                targets = sorted_ind(2:max_deg+1);
                E(max_ind,targets)=1; E(targets,max_ind)=1;
                deg_seq(targets)=deg_seq(targets)-1; deg_seq(max_ind)=0;
            end
        end
        % Now flip edges randomly to get a random graph of the ensamble
        flip_iters = 1; % Choose how many flips to perform. This is multiplied by edges number!!!
        [E, actual_flips]=flip_random_edges(E, flip_iters);
        return;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate d-regular graph. Simply call (using recrsion) the routine
        % generating graph with a given degree sequence
    case graph_REGULAR
        deg_seq = zeros(1,N)+param;
        E=generate_random_graph(graph_DEGSEQ, N, deg_seq);
        return;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end % switch












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function: Flip random couples of edges of E
function [E_out actual_flips] = flip_random_edges(E, flip_iters)

% Get the edge indexes, edge&vertices number
[I J S] = find(triu(E)); N = size(E,1); M = length(I);
r_vec = floor(rand(2,flip_iters*M)*M)+1;

actual_flips=0;

% Loop to flip edges
for i=1:flip_iters*M
    % Save current configuration
    if(mod(i,N) == 1)
        saved_E = E; saved_J = J;
    end

    % Check that chosen random indexes can flip
    % First check we have four distinct vertices and then if there are no
    % 'interferring' edges
    I1 = I(r_vec(1,i));
    I2 = I(r_vec(2,i));
    J1 = J(r_vec(1,i));
    J2 = J(r_vec(2,i));
    flip_flag = (I(r_vec(1,i))-I(r_vec(2,i)))*(J(r_vec(1,i))-J(r_vec(2,i)))*(I(r_vec(1,i))-J(r_vec(2,i)))*(J(r_vec(1,i))-I(r_vec(2,i)));
    flip_flag = flip_flag * (1-E(I(r_vec(1,i)),J(r_vec(2,i)))) *  (1-E(J(r_vec(1,i)),I(r_vec(2,i))));

    % Now do the flip
    if(flip_flag)
        actual_flips = actual_flips+1;

        % Flip in the graph matrix
        E(I(r_vec(1,i)),J(r_vec(1,i))) = 0;
        E(I(r_vec(2,i)),J(r_vec(2,i))) = 0;
        E(I(r_vec(1,i)),J(r_vec(2,i))) = 1;
        E(I(r_vec(2,i)),J(r_vec(1,i))) = 1;
        % The symmetric part
        E(J(r_vec(1,i)),I(r_vec(1,i))) = 0;
        E(J(r_vec(2,i)),I(r_vec(2,i))) = 0;
        E(J(r_vec(1,i)),I(r_vec(2,i))) = 1;
        E(J(r_vec(2,i)),I(r_vec(1,i))) = 1;
        % Flip in the IJ vectors, here only one side is used
        J(r_vec(1,i)) = J2;
        J(r_vec(2,i)) = J1;

        % Check if graph remained connected. If not, return to previous
        % saved graph
        if(mod(i,N) == 0)
            if(~is_connected(E))
                E=saved_E; J=saved_J; actual_flips = actual_flips-N;
            end
        end



    end
end

E_out = E;


