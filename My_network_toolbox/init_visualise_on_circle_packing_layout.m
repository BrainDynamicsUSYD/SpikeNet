function [ p, cmap, h ] = init_visualise_on_circle_packing_layout( X, Y, Radii )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

cmap = colormap(jet); % jet
close % close the figure opened by colormap()
% prepare figure
figwidth = 768;
figheight = 768;
figposition = [100, 100, figwidth, figheight];
h = figure('Name','Simulatioin Visualisation on  Circle-Packing Layout','NumberTitle','off','position',figposition, 'units','pixels');
axis equal; hold on;
set(gcf,'Renderer','OpenGL');  % to exploit grahpics cards (using built-in OpenGL in graphic card, i.e., firmware), making things magnitudes faster
% set(gca,'Color',[0 0 0]); % black background for better contrast
% white background seems better since the main idea is "the more active, the more visible"

% prepare scatter plot
max_size = 200;
s = pi*Radii.^2/max(pi*Radii.^2)*max_size;
p = scatter(X,Y,s,'filled');


end

