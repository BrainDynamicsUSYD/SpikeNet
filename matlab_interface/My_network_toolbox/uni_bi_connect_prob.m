function [p, bi_p, uni_p, z_bi, z_uni] =  uni_bi_connect_prob(A)

A = double(A > 0);
[n,~] = size(A);
N = n*(n-1); % total number of  pairs

p = nnz(A(:))/(n*(n-1));

if sum(A(logical(eye(n)))) > 0
    warning('There is self connection!')
end


A = A + A';


bi_n = sum(A(:) == 2);
bi_p = bi_n/(n*(n-1));
uni_n = sum(A(:) == 1);
uni_p = uni_n/(n*(n-1));

% Given connection probability p and total number of pairs N, the expected
% number of unconnected pairs should be N*(1-p)^2. The expected number of
% unidirectionally connected pairs should be 2*N*p*(1-p), and the expected
% number of bidirectionally connected pairsshould be N*p^2.

z_bi = bi_n / (N*p^2);
z_uni = uni_n / (2*N*p*(1-p));

% testing
% [p, bi_p, uni_p, z_bi, z_uni] =  uni_bi_connect_prob(MyRandomGraphGenerator('E_R'))
end