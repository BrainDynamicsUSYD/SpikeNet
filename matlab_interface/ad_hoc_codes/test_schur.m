function test_schur()
% test Schur decomposition

Kmat = [2.4*0.6  1.4;
        4.5  5.7*0.8 ]*10^-3; % miuSiemens
Pmat = [0.2 0.5;
        0.5 0.5];
rr = 0.6;
Mnum = 8;
N = [4000;1000];

% Note that in my notation, (i,j) means synapse goes FROM i to j.
%11
P0 = Pmat(1,1);
Msize = Mnum_2_Msize(Mnum, N(1));
[P, CL] = inter_module_Pmatrix(Msize, P0, rr);
P(CL==4) = P(CL==4)*0.5; % highest level connectionP
A11 = P_2_A(P,Msize);
[I11, J11, ~] = find(A11);
K11 = ones(size(I11))*Kmat(1,1);
% 12
[I12, J12, ~] = find(MyRandomGraphGenerator('E_R_pre_post','N_pre',N(1),'N_post',N(2),'p',Pmat(1,2)));
J12 = J12 + N(1);
K12 = ones(size(I12))*Kmat(1,2);
% 21
[I21, J21, ~] = find(MyRandomGraphGenerator('E_R_pre_post','N_pre',N(2),'N_post',N(1),'p',Pmat(2,1)));
I21 = I21 + N(1);
K21 = ones(size(I21))*Kmat(2,1)*-1;
% 22
[I22, J22, ~] = find(MyRandomGraphGenerator('E_R_pre_post','N_pre',N(2),'N_post',N(2),'p',Pmat(2,2)));
I22 = I22 + N(1);
J22 = J22 + N(1);
K22 = ones(size(I22))*Kmat(2,2)*-1;

% Now switch back to the convection that synapse goes from j to i
J = [I11; I12; I21; I22];
I = [J11; J12; J21; J22];
K = [K11; K12; K21; K22];
A = full(sparse(I,J,K));



disp('this may take a minute...');
[U, T0]=schur(A); % 5000x5000 may take about 1min

[S, T] = sep_schur(T0);

%  [U,T] = schur(X) produces a quasitriangular Schur matrix T and
%   a unitary matrix U so that X = U*T*U' and U'*U = EYE(SIZE(U)).

end



function [S, T] = sep_schur(T0)
% T0 is given by [U,T0]=schur(A);
%
% output S is block diagonal with 
% a) 2x2 blocks containing complex conjuate pairs of eigenvalus as its own eigenvalues 
% b) 1x1 blocks the real eigenvalues
%
% T is strictly upper-triangular.
%
% S + T = T0

[N,~] = size(T0);

S = zeros(N,N);

for i = 1:N
    S(i,i) = T0(i,i);
    if i < N && T0(i+1,i) ~= 0
        S(i+1,i) = T0(i+1,i);
        S(i,i+1) = T0(i,i+1);
    end
end

T = T0 - S;

end



