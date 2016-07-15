function [ cc, lag ] = semmat_autocc( sym_seq, semmat, lag )
% auto-correlation coefficients for symbolic sequence based on  symbolic
% multiplication defined by given semantic relatedness matrix of the
% alphabet 


if nargin == 2
    lag = 1:100;
end

cc = zeros(size(lag));
sym_seq = sym_seq(:)'; % row vector



for i = 1:length(lag)
    
    l = lag(i);
    
    % symbolic multiplication defined by semantic relatedness matrix of the alphabet
    cc_tmp = [];
    for j = 1:(length(sym_seq) - l)
        cc_tmp = [cc_tmp semmat( sym_seq(j), sym_seq(l+j)  )] ;  % table look-up
    end
    
    cc(i) = mean(cc_tmp);
end


end
