% Test graph layout using sub-optimal circle packing
function TestCirclePackingLayout()

clc;close all;clear;


% Input data
% Descending-sorted vector of radii (plural of radius)
% Radii can be directed mapped from strength of graph nodes
N = 10000; % number of nodes
Radii = rand(1,N);

% Output circle centre coordinates X, Y
tic;
[X,Y] = CirclePackingLayout(Radii);
toc;

% Plot Final result
PlotCircles(X,Y,Radii);

end

%-------------------------------------------------------------------------%
% Plot circles
function PlotCircles(X,Y,Radii)
figure('NumberTitle','Off','Name','Circles','position',[100,100,900,900]);
axis equal;
hold on;

N = length(X);
theta = linspace(0,2*pi,100);
dx = cos(theta);
dy = sin(theta);
for i = 1:N
    plot(X(i)+dx*Radii(i),Y(i)+dy*Radii(i));
end

end
























