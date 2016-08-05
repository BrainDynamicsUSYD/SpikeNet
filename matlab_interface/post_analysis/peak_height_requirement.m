function [ seq ] = peak_height_requirement( seq, data, peak_min )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

 cutoff = 1;
 [~, one_du, ~, one_start, ~] = seq_postprocess(seq, 1, cutoff);
 
 for i = 1:length(one_du)
     if max(data(one_start(i):one_start(i)+one_du(i)-1)) < peak_min
         seq(one_start(i):one_start(i)+one_du(i)-1) = 0;
     end
 end
end

