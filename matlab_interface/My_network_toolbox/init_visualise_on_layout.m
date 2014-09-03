function [ p, cmap, h ] = init_visualise_on_layout( X,Y,Z )
%initialisation for visualise_simulation()
%   Detailed explanation goes here

cmap = colormap(jet);
close % close the figure opened by colormap()
% prepare figure
figwidth = 768;
figheight = 768;
figposition = [100, 100, figwidth, figheight];
h = figure('Name','Simulatioin Visualisation on  Layout','NumberTitle','off','position',figposition, 'units','pixels');
axis equal; hold on;
set(gcf,'Renderer','OpenGL');  % to exploit grahpics cards (using built-in OpenGL in graphic card, i.e., firmware), making things magnitudes faster
view(32,-18);
% set(gca,'Color',[0 0 0]); % black background for better contrast
% white background seems better since the main idea is "the more active, the more visible"

% prepare patches
r = 0.02*(max(X)-min(X))*ones(length(X),1);
[ p ] = patchmarker(X,Y,Z,r, 'FaceAlpha', 0.1);

end

