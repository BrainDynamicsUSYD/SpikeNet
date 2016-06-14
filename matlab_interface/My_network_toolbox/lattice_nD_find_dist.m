function [dist] = lattice_nD_find_dist(Lattice, hw, i)
    % dist(j) is the Euclidean distance between node i and node j in the Lattice
    % Lattice
    %
    % hw is the half-width of the lattice
    %
    % the boundary condition is periodic
    %
    % p0 is the coordinates (x,y) of any point (not necessarily on the lattice)
    [N, D] = size(Lattice);
    p0 = Lattice(i,:);
    dist = zeros(N,1);
    for d = 1:D
        Xd = mod(Lattice(:,d)-p0(d)+hw, 2*hw+1) - hw; % periodic boundary condition
        dist = dist + Xd.^2;
    end
    dist = dist.^0.5; % Euclidean distance
end