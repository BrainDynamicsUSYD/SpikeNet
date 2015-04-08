function [ cc, lag ] = es_acf( sym_seq, lag )
% equal-symbol auto-correlation function
% reference: Voss, R.F., 1992, Evolution of Long-Rang Fractal Correlation
% and 1/f Noise in DNA Base Sequences
% sym_seq: symbolic sequence


if nargin == 1
    lag = 1:100;
end

cc = zeros(size(lag));
sym_seq = sym_seq(:)'; % row vector
for i = 1:length(lag)
    l = lag(i);
    % equal-symbol multiplication
    cc(i) = mean(sym_seq(l+1:end) == sym_seq(1:end-l)); % the bigger the lag, the less the data!
end


end

