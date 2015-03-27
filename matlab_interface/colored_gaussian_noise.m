function [out,err] = colored_gaussian_noise(acf, Nsamples)
[a,err] = ac2poly(acf);
% err is the expected error in the filter - can be quite high
% generate gaussian noise
noise = randn(1,Nsamples);
% apply via internal Matlab filter
filter_output = filter(1,a,noise);
% normalisation is total power, from frequency response
norm_factor = mean(abs(freqz(1,a,Nsamples)).^2);
% max(imag(filter_output)./real(filter_output)) %possible small complex component
% Take square root of power normalisation to obtain unit variance
out = real(filter_output) / sqrt(norm_factor); 

% variance and mean
out = out - mean(out);
out = out/std(out);
end
