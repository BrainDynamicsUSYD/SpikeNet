function [pairs] = rand_unique_pairs(N, pair_num, ordered, self_pair)
% This function randomly generates unique integer number pairs
% either ordered or unordered.
% 
% Input arguments:
%          N: defines the range of the numbers (from 1 to N)
%   pair_num: the number of pairs 
%    ordered: 1 for ordered pairs and 0 for unordered. The default is 0.
%  self_pair: 1 for allowing self-pair (e.g, 3-3 is a self-pair). The
%             default is 0. 
%
% Output argument:
%      pairs: a 2-by-pair_num matrix with each column containing a unique pair
%
% Yifan Gu, 6-Feb-2017
% yigu8115@gmail.com


% Input check
if nargin == 2
    ordered = 0;
end
if nargin == 3
    self_pair = 0;
end


% Generate pairs
if ordered == 0 % (i,j) is the same as (j,i)
    max_pair_num = N*(N-1)/2;
    if pair_num > max_pair_num
        pair_num = max_pair_num;
        warning('pair_num is larger than allowed!');
    end
    NN = ones(N);onnection
    if self_pair == 0
        NN(logical(eye(N))) = 0; % set diagonal entries to zero, no self-pair
    end
    NN = triu(NN); % extract upper triangular part for non-ordered pair
    [I,J] = find(NN);
    ind = randperm(length(I),pair_num);
    pairs = [I(ind)'; J(ind)'];

elseif ordered == 1
    max_pair_num = N*(N-1);
    if pair_num > max_pair_num
        pair_num = max_pair_num;
        warning('pair_num is larger than allowed!');
    end
    NN = ones(N);
    if self_pair == 0
        NN(logical(eye(N))) = 0; % set diagonal entries to zero, no self-connection
    end
    [I,J] = find(NN);
    ind = randperm(length(I),pair_num);
    pairs = [I(ind)'; J(ind)'];
end
    

end
