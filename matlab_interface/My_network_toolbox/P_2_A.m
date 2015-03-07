function A = P_2_A(P, Msize, non_self)
% Given inter-modular connection probability matrix P
% and module size vector Msize, generate full random connection matrix
% based on E-R random graph.
%
% By defualt, no_self = 1 and all of the self-connection entries are
% removed



% Generate node module index
Mnum = length(Msize);
Msize = Msize(:)'; % row vector
modules = cell2mat(arrayfun(@(x, y) repmat(x, [1 y]), 1:Mnum, Msize, 'UniformOutput', false));

% Generate connection matrix according to P matrix
[I,J,~] = find(rand(sum(Msize)) <= P(modules,modules));
A = logical(sparse(I,J,ones(size(I)),sum(Msize),sum(Msize)));

% No self-connection
if nargin < 3
    non_self = 1;
end

if non_self == 1
    A(logical(eye(size(A)))) = false; 
end

end
