function [A] = generate_graph_given_in_out_degree(degree_in, degree_out, no_self)

if sum(degree_in) ~= sum(degree_out)
    error('The sums of in-degree and out-degree do not match!')
end

if nargin == 2
    no_self = 1; % no self-connection
end

N = length(degree_in);
degree_in = degree_in(:); % column vector
degree_out = degree_out(:);

A = NaN;

t_max = 10;
trials = 0;

while trials <= t_max
    trials = trials + 1;
    failed = 0;
    
    I = []; % pre node index
    J = []; % post node index
    
    degree_in_left = degree_in;
    for i = randperm(N)
        
        % degree_in factor
        in_factor = degree_in_left;
        in_factor(in_factor <= 4) = in_factor(in_factor <= 4).^10/(4^10); % try to solve the issue: sometimes what's left will be [NaN NaN 1 1 2 1 NaN] and 5
        
        % take care of self connection
        self_factor = ones(N,1);
        if no_self == 1
            self_factor(i) = NaN;
        end
        
        % establish connections
        [~, ind] = sort( rand(N,1)./(in_factor.*self_factor), 'ascend' );
        chosen_j = ind(1: degree_out(i) );
        %Issue: sometimes what's left will be [NaN NaN 1 1 2 1 NaN] and 5
        %connections are needed. Try to fix this?
        if sum(isnan(degree_in_left(chosen_j))) ~= 0
            failed = 1;
            break;
            
        end
        
        degree_in_left(chosen_j) = degree_in_left(chosen_j) - 1;
        degree_in_left(degree_in_left == 0) = NaN;
        
        % store results
        chosen_j = chosen_j(:); % column vector
        I = [I; i*ones(size(chosen_j)) ]; %#ok<AGROW>
        J = [J; chosen_j]; %#ok<AGROW>
    end
    
    if failed == 0
        A = sparse(I,J,1,N,N);
        break
    end

end

  
  
  