
function [ Kernel ] = spike_train_kernel_YG( width, dt, choice )
%UNTITLED2 Summary of this function goes here
%   width and dt in msec
%   convolution estimates instantaneous firing rates (Hz)

switch lower(choice)
    % Square kernel
    case 'square'
    square_width = width;
    SquareKernelLength = round(square_width/dt);
    UnitSquareKernel = ones(1,SquareKernelLength)/SquareKernelLength;
    Square_kernel = UnitSquareKernel * (SquareKernelLength/(square_width/1000)); % so that convolution estimates instantaneous firing rates (Hz)
    Kernel = Square_kernel;
    
    % Gaussian filter
    case 'gaussian'
    sigma_gaussian = width; % ms,  width here is std of Gaussian function!
    num_sigma_gaussian = 4; % 3?
    sigma_gaussian_steps = round(sigma_gaussian/dt);
    t_k = -num_sigma_gaussian*sigma_gaussian_steps:num_sigma_gaussian*sigma_gaussian_steps;
    unitkernel = 1/(sqrt(2*pi)*sigma_gaussian_steps)*exp(-t_k.^2/(2*sigma_gaussian_steps^2));
    Gaussian_kernel = unitkernel * (length(t_k)/(2*num_sigma_gaussian*sigma_gaussian/1000)); % so that convolution estimates instantaneous firing rates (Hz)
    Kernel = Gaussian_kernel;
end

end

