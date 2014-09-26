function [ CC, RateA, RateB ] = CorrCoefYG( SpikeTrainA, SpikeTrainB, kernel, time_series )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
% SpikeTrainA/B is two sparse logical row vectors, coloumn index is time step
% [0 0 1 1 1 0 1 0 0 0;
%  1 0 0 0 0 0 1 0 0 0]
%
% Kernel is row vector, column index is time step, e.g.
% [0.1  0.2  0.4  0.2  0.1]
%  -2dt  -dt  0   +dt  +2dt
%
% TestCorrCoefYG.m:
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% SpikeLength = 10^6;
% SpikeNum = 10^2;
% kernelLength = 10^2;
% 
% SpikeTrainA = sparse(1, SpikeLength);
% SpikeTrainA(randperm(SpikeLength,SpikeNum)) = true;
% 
% SpikeTrainB = sparse(1, SpikeLength);
% SpikeTrainB(randperm(SpikeLength,SpikeNum)) = true;
% 
% kernel = ones(1,kernelLength)/kernelLength;
% 
% tic;
% [ CC ] = CorrCoefYG( SpikeTrainA, SpikeTrainB, kernel )
% toc;
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

kernel = kernel(:)'; % row vector

if length(kernel) == 1 && kernel(1) == 0 % no kernel required
    RateA = SpikeTrainA;
    RateB = SpikeTrainB;
else
    RateA = SpikeTrainConvolve(SpikeTrainA, kernel);
    RateB = SpikeTrainConvolve(SpikeTrainB, kernel);
end



if nargin == 3 || (nargin == 4 && time_series == 0) %
    CC = COV(RateA,RateB)/sqrt(COV(RateA,RateA)*COV(RateB,RateB));
    % CC_matlab = corrcoef([RateA(:) RateB(:)])
elseif nargin == 4 && time_series == 1 % return the time series
    CC =  (RateA-mean(RateA)).*(RateB-mean(RateB))/sqrt(COV(RateA,RateA)*COV(RateB,RateB));
end
    

end

% Covariance
function out = COV(V1,V2)
out = mean(V1.*V2)-mean(V1)*mean(V2);
% note that cov(x,x) = std(x)^2 = variance(x)
end
