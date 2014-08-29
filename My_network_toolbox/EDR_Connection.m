function [ K_matrix ] = EDR_Connection( density, lambda, Num_iteration, C_matrix, D_vector )
%Apply Exponential Distance Rule
%   Detailed explanation goes here
% density: connection density
% lambda: distance constant
% Num_iteration controls the balance between 
% a) approximation level of the method here compared with the published one
% b) computational expense
% if Num_iteration == M_tot, the same as published result

N = length(C_matrix(:,1)); % total number of neurons
P_vector = zeros(1,length(D_vector)); % probability vector
for i = 1:length(D_vector)
    P_vector(i) = exp(-D_vector(i)/lambda); % exponential rule
end
P_vector = P_vector/sum(P_vector); % Normalisation wrt distance
P_matrix = ones(N,1)*P_vector/N; % probability uniformly broken down to each pair


M_tot = N^2*density; % total number of connection
M_count = 0;
threshold = P_matrix*(M_tot/Num_iteration); % if use P as threshold, only one pair will be connected at each iteration
max_threshold = max(max(threshold));
if max_threshold >= 0.5
    disp('Num_iteration is too small!');
end
% if Num_iteration/M_tot >= 0.1
%     disp('Num_iteration is too big!');
% end

% Generate the network
[x,y] = size(C_matrix);
K_matrix = zeros(x,y);
while_count = 1;
while M_count < M_tot
    RandomMatrix = rand(x,y);
    K_matrix = K_matrix+(RandomMatrix < threshold);
    M_count = sum(sum(K_matrix));
    while_count = while_count+1;
end
while_count
Num_iteration
disp('final while_count should be almost the same as Num_iteration');
max_K = max(max(K_matrix));
K_matrix = K_matrix/max_K; % normalise K, entry range[0,1]

end

