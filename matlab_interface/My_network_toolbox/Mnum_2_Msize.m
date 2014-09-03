function Msize = Mnum_2_Msize(Mnum, N)
% Generate vector of first-level module size
Msize = zeros(1,Mnum);
Msize(1:end-1) = round(N/Mnum);
Msize(end) = N-sum(Msize);
Msize = Msize(:); % make sure it's column vector
end