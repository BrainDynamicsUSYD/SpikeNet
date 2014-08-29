% Test graph-layout algorithms for neuron networks
% Yifan GU, Complex System, School of Physics, USYD
% 06-Jan-2014
% yigu8115@uni.sydney.edu.au

function TestGraphLayout_Yifan()
clc;clear all; close all;

% Generate Adjacency matrix
% Creat Homogeneous Gaussian Connection for 2D square neuron sheet
a = 25; % n-by-n gird of neurons
%a = 20,30,40...gives seemingly different layout (but highly likely the same in high-dimensional space)
Space_time.n = a; 
[ Connect ] = Generate_circular_coupling ( Space_time, 'periodic', 'gaussian' ); % 2D Gaussian mask: peak strength = 1, sigma = 2 neurons

% Re-format the data into adjacency matrix (A)
C = Connect.C_E;
K = Connect.K_E;
[C_i,~,C_j] = find(C);
[~,~,K_ij] = find(K);
A_sparse = sparse(C_i,C_j,K_ij,a^2,a^2); % sparse adjacency matrix
A = full(A_sparse);
clear a Space_time Connect C K C_i C_j K_ij

% Calcualte Lap, Lapn, Laps
%number of nodes
n = size(A,1)
Iden = eye(n);  % Identity matrix
%degree matrix
Deg = diag(sum(A));   %  sum each Col,  then place on diagonals of a matrix 
%Laplacian matrix
Lap = Deg - A;
% % simple degree-normalized graph Laplacian - wrt row averageas - asymm. 
%   % Lapn = inv(Deg)*A + Iden;  % nb. this is working from the "unsigned" def'n of L = D+A
% Lapn = Deg\Lap; %inv(Deg)*Lap;  % use L=D-A (as above) - the is L-rw 

% symmetrise this - via deg-i, pre&post multiply :
Degh = Iden;  Deginv = eye(n,n);% start w identity matrix
for i=1:n
    Degh(i,i) = 1/sqrt(Deg(i,i));
    Deginv(i,i) = 1/(Deg(i,i));
end
Laps = (Degh*Lap)*Degh;  % this is Lap-sym

% check Laps
figure('Name','1','NumberTitle','off');
spy(Laps);
title('Spy(), SymLap','FontWeight','bold','FontSize',14 )
%rank(Laps);

%% EVD with Laps
% EigVecs {ie the resultant coords, x} also degree normalised
[vx,~] = eig(Laps);  % , 'nobalance');  %  EVD may be poor 

%  prepare the plot coords
xevd = diag(Degh).*vx(:,2);  % normalise each Node coord by its Degree (nb. ".*")  ?? not sqrt?
yevd = diag(Degh).*vx(:,3);  % mu-1 = 0; then next 3 smallest
zevd = diag(Degh).*vx(:,4);

% Plot EVD with Laps in 3D
figwidth = 768;
figheight = 768;
figposition = [100, 100, figwidth, figheight];
figure('Name','2','NumberTitle','off','position',figposition, 'units','pixels');
hold on
view(47,38);

% Plot the nodes
for ii = 1:n
    plot3(xevd(ii),yevd(ii),zevd(ii),'.','Color',[0 0 0.7],'MarkerSize',15)
end
xlabel('EigVec-2','FontWeight','bold','FontSize',14)
ylabel('EigVec-3','FontWeight','bold','FontSize',14)
zlabel('EigVec-4','FontWeight','bold','FontSize',14)
title('SymLap - EVD - in 3D','FontWeight','bold','FontSize',14 )

% Plot the connections (too many!? not working?)
% [ii,jj] = find(A);
% line([x(ii)'; x(jj)'],[y(ii)'; y(jj)'], [z(ii)'; z(jj)'],'Color', [0.65 0.65 0.65],'LineWidth',0.25)

% show/animate the transformation from 2D square sheet to 3D EVD layout
for ii = 2:n
    line([xevd(ii-1); xevd(ii)],[yevd(ii-1); yevd(ii)], [zevd(ii-1); zevd(ii)],'Color', [0.65 0 0],'LineWidth',0.5)
    % pause(0.01);
end

%% SVD with Laps
[ux,~,~] = svd(Laps);

xsvd = diag(Degh).*ux(:,n-1);  %  take smallest EigVals (here the later ones) drop 1 x zero
ysvd = diag(Degh).*ux(:,n-2);  %  & normalise by  sqrt(deg-i) - not deg-i
zsvd = diag(Degh).*ux(:,n-3);

% plot in 3D
% Plot EVD with Laps in 3D
figwidth = 768;
figheight = 768;
figposition = [100, 100, figwidth, figheight];
figure('Name','3','NumberTitle','off','position',figposition, 'units','pixels');
hold on
view(59,32);

% Plot the nodes
for ii = 1:n
    plot3(xsvd(ii),ysvd(ii),zsvd(ii),'.','Color',[0 0 0.7],'MarkerSize',15)
end
xlabel('EigVec-2','FontWeight','bold','FontSize',14)
ylabel('EigVec-3','FontWeight','bold','FontSize',14)
zlabel('EigVec-4','FontWeight','bold','FontSize',14)
title('SymLap - SVD - in 3D','FontWeight','bold','FontSize',14 )

% Plot the connections (too many!? not working?)
% [ii,jj] = find(A);
% line([x(ii)'; x(jj)'],[y(ii)'; y(jj)'], [z(ii)'; z(jj)'],'Color', [0.65 0.65 0.65],'LineWidth',0.25)

% show/animate the transformation from 2D square sheet to 3D EVD layout
for ii = 2:n
    line([xsvd(ii-1); xsvd(ii)],[ysvd(ii-1); ysvd(ii)], [zsvd(ii-1); zsvd(ii)],'Color', [0.65 0 0],'LineWidth',0.5)
    % pause(0.01);
end

%% SDE ( Ali Civril'06)
%  uses Paths
Paths=zeros(n); Pred=zeros(n);
h_wb= waitbar(0,'Please wait for breath(A,i)...');
for i=1:n
    [Paths(i,:), Pred(:,i)] = breadth(A,i);
    Paths(i,i)=0;  % need to eliminate souce (spurious entry??)
    waitbar(i/n);
end
close(h_wb);
% and need CleanPaths to elim. Inf's
onevec = ones(1,n);
L = Paths.^2; % array of (d-ij)^2, graph dist's.
gamma = eye(n)- onevec'*onevec/n;  % form projection matrix, to get centered coords (mean 0) 
M=-gamma*L*gamma/2;  %  normalised dist matrix (of graph path lengths)
[ux,w,~] = svd(M);
% for g-unweighted:
sl1=sqrt(w(1,1)); sl2=sqrt(w(2,2)); sl3=sqrt(w(3,3)); % sqrt(Lamda-1) etc
xsde = sl1*ux(:,1);  %  take largest 2 EigVals
ysde = sl2*ux(:,2);  
zsde = sl3*ux(:,3);
clear onevec;      % saves space

% Plot - SDE in 3D
figwidth = 768;
figheight = 768;
figposition = [100, 100, figwidth, figheight];
figure('Name','4','NumberTitle','off','position',figposition, 'units','pixels');
hold on
view(59,32);

% Plot the nodes
for ii = 1:n
    plot3(xsde(ii),ysde(ii),zsde(ii),'.','Color',[0 0 0.7],'MarkerSize',15)
end
xlabel('EigVec-2','FontWeight','bold','FontSize',14)
ylabel('EigVec-3','FontWeight','bold','FontSize',14)
zlabel('EigVec-4','FontWeight','bold','FontSize',14)
title('SymLap - SDE - in 3D','FontWeight','bold','FontSize',14 )

% Plot the connections (too many!? not working?)
% [ii,jj] = find(A);
% line([x(ii)'; x(jj)'],[y(ii)'; y(jj)'], [z(ii)'; z(jj)'],'Color', [0.65 0.65 0.65],'LineWidth',0.25)

% show/animate the transformation from 2D square sheet to 3D EVD layout
for ii = 2:n
    line([xsde(ii-1); xsde(ii)],[ysde(ii-1); ysde(ii)], [zsde(ii-1); zsde(ii)],'Color', [0.65 0 0],'LineWidth',0.5)
    % pause(0.01);
end

end



















%% Supplementary functions
function [ Connect ] = Generate_circular_coupling ( Space_time, BoundaryCondition, ConnectType, varargin )
%This function generate circular coupling for 2-D neuron circuits, which
%means that the coupling is homogeneous and isotropical.
% 

n = Space_time.n;

switch lower(ConnectType)
    case 'gaussian'
        % (Default) Parameters for Gaussian connections
        Para.J = 1; % strength
        Para.sigma = 2; % s.t.d.
        % Read assigned parameters if given
        for i = 1:length(varargin)/2
            temp = varargin{2*i};
            eval(['Para.', varargin{2*i-1}, '= temp;']);
        end
        % Range of connections included (beyond this point values are
        % negligible and it is computationally faster to ignore them)
        Para.extent = 3*Para.sigma;
        % Determine index connection matrix and according distance matrix
        [Connect_matrix, Distance_matrix] = CircularConnection( n, Para.extent, BoundaryCondition );
        % Calculate according connection strength matrix given by Mexican hat function
        % The Gaussian form is also used in (Wenhao Zhang & Si Wu, 2012)
        % Strength_matrix = Para.J/(Para.sigma*sqrt(2*pi))*exp(-0.5*Distance_matrix.^2/Para.sigma^2);
        Strength_matrix = Para.J*exp(-0.5*Distance_matrix.^2/Para.sigma^2);
        
    case 'mexhat'
        % Parameters for Mexican hat connections
        Para.WE = 0.24;
        Para.sigmaE = 20;
        Para.WI = 0.17;
        Para.sigmaI = 43;
        % Read assigned parameters if given
        for i = 1:length(varargin)/2
            temp = varargin{2*i};
            eval(['Para.', varargin{2*i-1}, '= temp;']);
        end
        % Range of connections included (beyond this point values are
        % negligible and it is computationally faster to ignore them)
        Para.extent = 15;
        % Determine index connection matrix and according distance matrix
        [ Connect_matrix, Distance_matrix ] = CircularConnection( n, Para.extent, BoundaryCondition );
        % Calculate according connection strength matrix given by Mexican hat function
        Strength_matrix = Para.WE*exp(-Distance_matrix.^2/Para.sigmaE) ...
            - Para.WI*exp(-Distance_matrix.^2/Para.sigmaI);

end % switch

% Seperate excitatory and inhibitory parts as required by the neuron model
[ Connect ] = SeperateEI( Connect_matrix, Strength_matrix );

% Record parameters
Connect.Para = Para;

end

function [ Ind_matrix , Distance_matrix] = CircularConnection( n, radius, BoundaryCondition )
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
Distance_matrix = M*dist;

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

function [ Connect ] = SeperateEI( Connect_matrix, Strength_matrix )
% This function seperates excitatory and inhibitory components of the
% connections matrix and according strength matrix

n = length(Connect_matrix(:,1));
K_typical = Strength_matrix(1,:);
pos = find(K_typical>0);
neg = find(K_typical<0); 

% Plot the strength distribution

% Create a separate matrices for excitatory connections, with each row
% corresponding to the index of the pre/post-synaptic neurons 
lengthE = length(pos);
C_E = Connect_matrix(Strength_matrix>0);% Find excitatory indices
C_E = reshape(C_E,n,lengthE);% Reshape into matrix
K_E = Strength_matrix(Strength_matrix>0);% Find excitatory strengths
K_E = reshape(K_E,n,lengthE);% Reshape into matrix

% Create separate matrices for inhibitory connections, with each row
% corresponding to the index of the pre/post-synaptic neurons 
lengthI = length(neg);
C_I = Connect_matrix(Strength_matrix<0);% Find inhibitory indices
C_I = reshape(C_I,n,lengthI);% Reshape into matrix
K_I = -Strength_matrix(Strength_matrix<0);% Find inhibitory strengths Note that the strength is positive value
K_I = reshape(K_I,n,lengthI);% Reshape into matrix

% % Sort C and K matrix
% for i = 1:n^2
%     [sort_result, sort_ind] = sort( C_E(i,:),'ascend');
%     C_E(i,:) = sort_result;
%     K_E(i,:) = K_E(i,sort_ind);
%     [sort_result, sort_ind] = sort( C_I(i,:),'ascend');
%     C_I(i,:) = sort_result;
%     K_I(i,:) = K_I(i,sort_ind);
% end

% Record outputs into structure for code simplisity
Connect.C_E = C_E;
Connect.K_E = K_E;
Connect.C_I = C_I;
Connect.K_I = K_I;

end

function [distance,branch] = breadth(CIJ,source)
%BREADTH        Auxiliary function for breadthdist.m
%
%   [distance,branch] = breadth(CIJ,source);
%
%   Implementation of breadth-first search.
%
%   Input:      CIJ,        binary (directed/undirected) connection matrix
%               source,     source vertex
%
%   Outputs:    distance,   distance between 'source' and i'th vertex
%                           (0 for source vertex)
%               branch,     vertex that precedes i in the breadth-first search tree
%                           (-1 for source vertex)
%        
%   Notes: Breadth-first search tree does not contain all paths (or all 
%   shortest paths), but allows the determination of at least one path with
%   minimum distance. The entire graph is explored, starting from source 
%   vertex 'source'.
%
%
%   Olaf Sporns, Indiana University, 2002/2007/2008

N = size(CIJ,1);

% colors: white, gray, black
white = 0; 
gray = 1; 
black = 2;

% initialize colors
color = zeros(1,N);
% initialize distances
distance = inf*ones(1,N);
% initialize branches
branch = zeros(1,N);

% start on vertex 'source'
color(source) = gray;
distance(source) = 0;
branch(source) = -1;
Q = source;

% keep going until the entire graph is explored
while ~isempty(Q)
   u = Q(1);
   ns = find(CIJ(u,:));
   for v=ns
% this allows the 'source' distance to itself to be recorded
      if (distance(v)==0)
         distance(v) = distance(u)+1;
      end;
      if (color(v)==white)
         color(v) = gray;
         distance(v) = distance(u)+1;
         branch(v) = u;
         Q = [Q v];
      end;
   end;
   Q = Q(2:length(Q));
   color(u) = black;
end
end
