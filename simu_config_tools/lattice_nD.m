function [Lattice, N] = lattice_nD(D, hw)


% % parameters
% D = 2; % dimension of embedding space (regular square lattice)
% hw = 10; % half-width of the embedding space, width is 2*hw+1


% put nodes on a regular square lattice
if D == 2
    [X1,X2] = ndgrid(-hw:hw); % coordinates
    Lattice = [X1(:) X2(:)]; % every row is the coordinate of one node
elseif D == 3
    [X1,X2,X3] = ndgrid(-hw:hw); % coordinates
    Lattice = [X1(:) X2(:) X3(:)]; % every row is the coordinate of one node
end

% show number of nodes
N = length(Lattice(:,1));

end