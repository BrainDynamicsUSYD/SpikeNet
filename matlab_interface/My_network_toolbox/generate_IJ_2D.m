function [ I, J, dist_IJ, iter_hist ] = generate_IJ_2D( degree_in_0, degree_out_0, tau_d, cn_scale_wire, iter_num )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


hw = (sqrt(length(degree_in_0)) - 1)/2;
D = 2; % Pulin: 2D is enough
[Lattice, N] = lattice_nD(D, hw);

show_wait_bar = 1;

if nargin == 4
    iter_num = 5;% 5 is arbitrary, try something else?
end

for iter = 1:iter_num
    fprintf('Generate I and J: Iteration = %d \n', iter)
    if show_wait_bar == 1
        wb_h = waitbar(0,'Please wait...');
        wb_p = 0;
    end
    I = []; % pre node index
    J = []; % post node index
    dist_IJ = [];
    degree_in_left = degree_in_0;
    for i = randperm(N)
        if show_wait_bar == 1
            wb_p = wb_p + 1;
            waitbar(wb_p / N)
        end
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
        dist_IJ_tmp = post_dist(chosen_j);
        dist_IJ  = [dist_IJ ; dist_IJ_tmp(:)]; %#ok<AGROW>
    end
    if show_wait_bar == 1
        close(wb_h)
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
    
end


end

