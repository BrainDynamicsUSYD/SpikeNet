function [ stub ] = degree2stub( degree )
% This function turn the in/out-degree vector into a stub vector
% The stub vector will be randomly permuted.
% Example:
%       if degree = [3 2 5 1]', for node i = 1, 2, 3, 4
%       the stub vector is
%          stub = [1 1 1 2 2 3 3 3 3 3 4]'; % before permutation
% Note that degree must be column vector!
% stub is also column vector.

Ind = 1:length(degree); % index vector for the nodes
stub = cell2mat(arrayfun(@(x, y) repmat(x, [1 y]), Ind, degree', 'UniformOutput', false));
stub = stub';

% random permutation
edge_num = length(stub);
ind_perm = randperm(edge_num);
stub = stub(ind_perm);

end