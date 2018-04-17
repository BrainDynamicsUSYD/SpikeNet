function [ I, J, dist_IJ, iter_hist, Lattice ] = generate_IJ_2D( degree_in_0, degree_out_0, tau_d, cn_scale_wire, iter_num, dist_cutoff, record_cc )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


hw = (sqrt(length(degree_in_0)) - 1)/2;
D = 2; % Pulin: 2D is enough
[Lattice, N] = lattice_nD(D, hw);

show_wait_bar = 0;

if nargin < 6
    dist_cutoff = Inf; 
end
if nargin < 7
    record_cc = 0;% 5 is arbitrary, try something else?
end

if cn_scale_wire < 1
    warning('cn_scale_wire cannot be less than 1!');
end

for iter = 1:iter_num
    fprintf('Generate I and J: Iteration = %d \n', iter)
    if show_wait_bar == 1
        wb_h = waitbar(0,'Please wait...');
        
    end
    wb_p = 0;
    I = []; % pre node index
    J = []; % post node index
    dist_IJ = [];
    degree_in_left = degree_in_0;
    for i = randperm(N)
        wb_p = wb_p + 1;
        if show_wait_bar == 1
            
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
        degree_in_factor(degree_in_factor <= 4) = degree_in_factor(degree_in_factor <= 4).^6/(4^6); % try to solve the issue: sometimes what's left will be [NaN NaN 1 1 2 1 NaN] and 5
        % joint factor
        %dist_factor = normc(dist_factor);
        %degree_in_factor = normc(degree_in_factor);
        %cn_factor = normc(cn_factor);
        joint_factor = dist_factor.*degree_in_factor.*cn_factor; % will be properly scaled later;
        joint_factor = normc(joint_factor);
        
        % distance cutoff
        if nargin == 6
            joint_factor(post_dist > dist_cutoff) = NaN;
        end
    
        % establish connections
        [~, ind] = sort( rand(N,1)./joint_factor, 'ascend' );
        chosen_j = ind(1: degree_out_0(i) );
        
        % The following way is wrong!!
%         [~, ind] = sort( rand(N,1).*joint_factor, 'descend' );
%         chosen_j = ind(1: degree_out_0(i) );
        
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
    cn_std = std(cn_dist);
    if record_cc == 1
        if hw > 31 % only calculate it for a circular region with radius less than 31 for computational reasons
            sample_ind = (Lattice(:,1).^2 + Lattice(:,2).^2) <= 31^2;
            ccoef = full(clustering_coef_wd(A(sample_ind,sample_ind)));
        else
            ccoef = full(clustering_coef_wd(A));
        end
    else
        ccoef = 0;
    end
    clear cn_dist;
    % store iteration history
    if iter == 1
        iter_hist.edges = {edges};
        iter_hist.Y = {Y};
        iter_hist.cn_minmax = {cn_minmax};
        iter_hist.cn_std = {cn_std};
        iter_hist.ccoef = {ccoef};
    else
        iter_hist.edges{end+1} = edges;
        iter_hist.Y{end+1} = Y;
        iter_hist.cn_minmax{end+1} = cn_minmax;
        iter_hist.cn_std{end+1} = cn_std;
        iter_hist.ccoef{end+1} = ccoef;
    end
    
end


end

