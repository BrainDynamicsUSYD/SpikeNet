function [ I,J ] = Lattice2Lattice( L_pre, L_post, hw, tau_c, p0, dist_cutoff )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


N_pre = length(L_pre(:,1));
N_post = length(L_post(:,1));
I = [];
J = [];
for i = 1:N_pre

    dist = lattice_nD_find_dist(L_post, hw, L_pre(i,1), L_pre(i,2) );
    dist_factor = exp(-dist/tau_c);
    % distance cutoff
    if nargin == 6
         dist_factor(dist > dist_cutoff) = NaN;
    end
    
    [~, ind] = sort( rand(N_post,1)./dist_factor, 'ascend' );
    chosen_j = ind(1:poissrnd(N_post*p0));
    I = [I; i*ones(size(chosen_j))]; %#ok<AGROW>
    J = [J; chosen_j]; %#ok<AGROW>
end

end

