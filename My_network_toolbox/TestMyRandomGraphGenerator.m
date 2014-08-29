% Test

clc;clear all;close all;
addpath(genpath('C:\Users\yigu8115\Desktop\My network toolbox'));
[ A ] = MyRandomGraphGenerator('arbitrary_degree_newman');

A = A./repmat(sum(A, 1), [length(A(:, 1)), 1]); % make the summation of in-degree weight one

% [ A ] = GenerateAjacencyMatrix( 'cat_all' );
% [~, ~, ~] = SVDLayout(A, 'ShowPlots', 1, 'ShowConnection', 1, 'Method', 'Laps');
% [~, ~, ~] = SVDLayout(A, 'ShowPlots', 1, 'ShowNode', 1, 'ShowModularity', 1, 'ShowDegree', 1, 'ShowConnection', 1, 'Method', 'Laps');
% MyGraphMeasurements( A )
