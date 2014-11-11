function [ Y ] = SpikeTrainConvolve(SpikeTrain, kernel)
% SpikeTrain is (sparse) row vectors, coloumn index is time step, e.g
%
% [0 0 1 1 1 0 1 0 0 0; for single neuron spike train
% [0 0 1 2 1 9 0 0 9 0] for neuron-cluster spike train
%
% Kernel is row vector, column index is time step, e.g.
% [0.1  0.2  0.4  0.2  0.1]
%  -2dt  -dt  0   +dt  +2dt

if mod(length(kernel),2) == 0
    kernel = [kernel 0];
end
kernel_length = length(kernel);
kernel_half_length = (kernel_length-1)/2;
kernel_t = -kernel_half_length:kernel_half_length;

[num, q] = size(SpikeTrain);
Y = zeros(num,q);



% if nnz(SpikeTrain) > 0.1*num*q % ???
%     
%     SpikeTrain = double(full(SpikeTrain));
%     for n = 1:num
%         Y(n,:) = conv( SpikeTrain(n,:), kernel_t, 'same');
%     end
%     
% else
    RatePadded = zeros(1,q+2*kernel_half_length);% padding
    for n = 1:num
        [~,timing,num_spikes] = find(SpikeTrain(n,:));
        if ~isempty(timing)
            for i = 1:kernel_length
                RatePadded(1, timing+kernel_half_length+kernel_t(i)) = RatePadded(1, timing+kernel_half_length+kernel_t(i))+num_spikes*kernel(i);
            end
        end
        Y(n,:) = RatePadded((kernel_half_length+1):(end-kernel_half_length));
    end
% end


end

