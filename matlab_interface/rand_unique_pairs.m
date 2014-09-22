function [pairs] = rand_unique_pairs(N, pair_num, ordered)
% given N, generate random pair (i,j)
% where 1 <= i, j <= N and i!=j
% output pairs is (2 x pair_num)
% and each pair is unique

if nargin == 2
    ordered = 0;
end

if ordered == 0 % (i,j) is the same as (j,i)
    max_pair_num = N*(N-1)/2;
    if pair_num > max_pair_num
        pair_num = max_pair_num;
        warning('pair_num is larger than allowed!');
    end
    NN = ones(N);
    NN(logical(eye(N))) = 0; % set diagonal entries to zero, no self-connection
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
    NN(logical(eye(N))) = 0; % set diagonal entries to zero, no self-connection
    [I,J] = find(NN);
    ind = randperm(length(I),pair_num);
    pairs = [I(ind)'; J(ind)'];
end
    

end
