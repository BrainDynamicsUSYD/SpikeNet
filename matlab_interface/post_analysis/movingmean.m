function y = movingmean(x, k)
   k = ones(1, k) / k;
   y = conv(x, k, 'same');
end