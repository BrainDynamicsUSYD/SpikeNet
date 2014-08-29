function [ Ind_matrix , Distance_vector] = CircularConnection( n, radius, BoundaryCondition )
% This function returns index connection matrix and according distance
% matrix;
% Input: >>"n" defines the size of the (n x n) grid
%        >>"radius" defines the circular range within which neurons are
%           connected.
%        >>"BoundaryCondition": free/periodic

% First build the subscript connection matrix and according distance matrix
index = (1:n^2)';
[x, y] = ind2sub(n,index);
range = -radius:radius;
M = length(range);
range2d = zeros(2,M^2);
for b=1:M
    for a=1:M
        range2d(:,(b-1)*M+a) = [range(a); range(b)];
    end
end
% Crop any self-connecting values that exceed the circular radius given by
% 'extent'
dist = sqrt(range2d(1,:).^2+range2d(2,:).^2);
crop = find(dist==0 | dist>radius);
dist(:,crop) = [];
range2d(:,crop) = [];
% Find subscripts of nearby neurons by adding the coordinates (x,y) of
% each target neuron to the 'stereotypical' distances given by range2d.
% Matrix multiplication by L and M to calculate for all target neurons
% without using a for loop
L = ones(1,size(range2d,2));
M = ones(length(x),1);
X_matrix = x*L+M*range2d(1,:);
Y_matrix = y*L+M*range2d(2,:);
% Distance_matrix = M*dist;
Distance_vector = dist;

% Then convert the subscript connection matrix to index matrix.
switch lower(BoundaryCondition);
    case 'free'
        X_matrix( X_matrix <= 0 ) = NaN;
        X_matrix( X_matrix >  n ) = NaN;
        Y_matrix( Y_matrix <= 0 ) = NaN;
        Y_matrix( Y_matrix >  n ) = NaN;
        Ind_matrix  = sub2ind([n n],X_matrix,Y_matrix);
        % Define an additional point (n^2+1) as the rubbish bin of all the out-of-boundary points
        Ind_matrix (isnan(Ind_matrix)) = n^2+1;
    case 'periodic'
        % Error check
        % Connection overlaping itself due to large radius is not allowed.
        if radius >= 0.5*n
            error('Radius of local circular connection is too large!')
        end
        X_matrix( X_matrix <= 0 ) = X_matrix( X_matrix <= 0 )+n;
        X_matrix( X_matrix >  n ) = X_matrix( X_matrix >  n )-n;
        Y_matrix( Y_matrix <= 0 ) = Y_matrix( Y_matrix <= 0 )+n;
        Y_matrix( Y_matrix >  n ) = Y_matrix( Y_matrix >  n )-n;
        Ind_matrix  = sub2ind([n n],X_matrix,Y_matrix);
end

end



