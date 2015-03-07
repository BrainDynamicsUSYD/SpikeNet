function [ syn ] = do_lesion_adhoc( syn, Msize, P, CL, lesion_level, lesion_ratio )
% syn_cell{1} = do_lesion(syn_cell{1}, Msize, P, CL, lesion_level, lesion_ratio);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Mnum = length(Msize);
Msize = Msize(:)';
M_a = 1 + [0, cumsum(Msize(1:end-1)) ];
M_b = cumsum(Msize);

post_a = repmat(M_a, Mnum, 1);
post_b = repmat(M_b, Mnum, 1);
pre_a = repmat(M_a', 1, Mnum);
pre_b = repmat(M_b', 1, Mnum);


% check P and CL
P0 = unique( P(CL == lesion_level) );
if length(P0) ~= 1
    error('something wrong with P(CL == lesion_level)');
end
[pres, posts, ~] = find(CL == lesion_level);


for ii = 1:length(pres)
    i = pres(ii); % for this syntax to work, must use row vector
    j = posts(ii); % for this syntax to work, must use row vector
    
    IJ_ind = find(syn.I >= pre_a(i,j) & syn.I <= pre_b(i,j) & ...
        syn.J >= post_a(i,j) & syn.J <= post_b(i,j));
    
    if lesion_ratio <= 1 && lesion_ratio > 0 % randomly delete synapses
        
        IJ_delete = IJ_ind( rand( size(IJ_ind) ) > lesion_ratio );
        syn.I( IJ_delete ) = [];
        syn.J( IJ_delete ) = [];
        syn.K( IJ_delete ) = [];
        syn.D( IJ_delete ) = [];
        
    elseif lesion_ratio > 1 % randomly add synapses
        
        % check if K0 is identical
        K0 = unique( syn.K(IJ_ind) );
        if length(K0) ~= 1
            error('Synapse strength is not identical');
        end
        %
        P1 = P0*lesion_ratio;
        P = (P1 - P0) / (1 - P0); % from (1-P0)*N non-synapses, we choose (P1 - P0)*N new ones
        syn_add = rand( pre_b(i,j)-pre_a(i,j)+1,  post_b(i,j)-post_a(i,j)+1 ) < P; % some of them are existing ones
        I0 = syn.I( IJ_ind );
        J0 = syn.J( IJ_ind );
        syn_existing = full( sparse(I0-pre_a(i,j)+1, J0-post_a(i,j)+1, true) );
        [I_add, J_add, ~] = find(syn_add & ~syn_existing);
        %
        I_add = I_add + pre_a(i,j) - 1;
        J_add = J_add + post_a(i,j) - 1;
        % remove self-connection
        self_c = (I_add == J_add);
        I_add(self_c) = [];
        J_add(self_c) = [];
        % add to syn
        syn.I = [syn.I;I_add];
        syn.J = [syn.J;J_add];
        syn.K = [syn.K;ones(size(I_add))*K0];
        syn.D = [syn.D;rand(size(I_add))];
        
    end
end


warning('Uniform distribution within [0 1] for delay are assumed!');

end

