function [centers, counts] = log_bin_histc(x, bin_num)


if nnz(x <= 0)
    x = x(x > 0);
    disp('Non-positive values ignored!');
end

bmin = log( min(x) );
bmax = log( max(x) );
bin_edge = linspace(bmin, bmax, bin_num);
bin_step = diff(bin_edge(1:2));
bin_edge = [bin_edge,  bin_edge(end)+bin_step];
centers = bin_edge(2:end) - bin_step/2; % mid-points in log scale
centers = exp(centers);

counts = histc( log(x), bin_edge );
counts = counts(1:end-1); % last value should be zero

%     N = histc(X,EDGES), for vector X, counts the number of values in X
%     that fall between the elements in the EDGES vector (which must contain
%     monotonically non-decreasing values).  N is a LENGTH(EDGES) vector
%     containing these counts.
%     
%     N(k) will count the value X(i) if EDGES(k) <= X(i) < EDGES(k+1).  The
%     last bin will count any values of X that match EDGES(end).  Values
%     outside the values in EDGES are not counted.  Use -inf and inf in
%     EDGES to include all non-NaN values.


end