function [ A ] = MyHierarchyGraph(varargin)
% 1) full binary tree: module number must be 2^n
% 2) equal module size
%
% % Example
% % overall connection probability
% P0 = 0.1;
% % hierarchy model parameter
% r = 0.7;
% % number of nodes in each module in first hierarchy
% N = [500; % N1
%      500; % N2
%      500; % N3
%      500];% N4
% % Hierarchy Group Label
% %     H1  H2  H3
% HGL = [1   1   1 ; %N1
%        2   1   1 ; %N2
%        3   2   1 ; %N3
%        4   2   1 ];%N4
% %
% % H3:             1
% %          ______________
% % H2:      1             2
% %      _________     _________
% % H1:  1       2     3       4
% %     N1       N2    N3      N4


% Default parameters
N = 4000; % total number of nodes
Mnum = 8; % number of modules in first hierarchy
P0 = 0.1;
r = 0.7;

% Read parameter setting (N, P0)
for i = 1:length(varargin)/2
    temp = varargin{2*i};
    eval([varargin{2*i-1}, '= temp;']);
end

% Generate vector of first-level module size
Msize = Mnum_2_Msize(Mnum, N);

% Generate inter modular connection probability matrix
[P, CL] = inter_module_Pmatrix(Msize, P0, r);
% % Plot P matrix
% imagesc(P);colorbar;

% Generate full connection matrix from P
A = P_2_A(P,Msize);

% Display
fprintf('Hierarchical Graph: N=%d, Mnum=%d, P0=%g, r=%g\n', N, Mnum, P0, r);

end




