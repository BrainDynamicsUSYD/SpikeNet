
clc;clear all;close all;

addpath(genpath('C:\Users\yigu8115\Desktop\My network toolbox'));

% generate network
[ A, ~,~ , ~ ] = GenerateAjacencyMatrix( 'cat_all');
% solve layout problem
[~, ~, ~] = SVDLayout(A, 'ShowPlots', 1, 'ShowConnection', 1);


% % generate network
% [ A, C_i, C_j, K_ij ] = GenerateAjacencyMatrix( 'homo_gaussian', 'a', 40, 'sigma', 5 );
% % solve layout problem
% [X, Y, Z] = SVDLayout(A, 'ShowPlots', 0);
% % simulatoin
% V_range = [0,1];
% for i = 1:100
%     V = rand(length(X),1);
%     if i == 1
%         % initialisation
%         [ p, cmap, h ] = init_visualise_simulation( X,Y,Z );
%     else
%         % visualisation
%         visualise_simulation( p, V, V_range, cmap )
%         pause(0.01);
%     end
% end



        
        