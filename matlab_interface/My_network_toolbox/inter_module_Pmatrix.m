function [P, CL] = inter_module_Pmatrix(Msize, P0, r)

Mnum = length(Msize); % number of modules in first hierarchy

HGL = Mnum_2_HGL(Mnum);



% % Lesion: use NaN in HGL
% HGL(:,end) = NaN;
% disp('Lesion in the hierarchical network!');

% r: change of inter-modular connnection PROBABILITY as hierarchy goes up. H1: P1, H2:P1*r, H3: P1*r^2,...
% Build the connection level matrix
CL = HGL_2_CL(HGL);

% Build the un-normalised inter-modular connection probability matrix P
P = zeros(Mnum,Mnum);
for i = 1:Mnum
    for j = 1:Mnum
        if CL(i,j) > 0 % 0 means no connection
            P(i,j) = r.^(CL(i,j)-1); % note the "-1" here
        end
    end
end


% Normalise P using overall probability P0
for j = 1:Mnum % note that we use i_pre, j_post
    P(:,j) = P(:,j)/mtimes(P(:,j)',Msize)*sum(P0*Msize);
end



% % Lesion: use NaN in HGL
% HGL(:,end) = NaN;
% CL = HGL_2_CL(HGL);
% Lesion = CL>0;
% P = P.*Lesion; % directly remove connection
% disp('Lesion in the hierarchical network!');


% check P matrix
if sum(sum(P>=1)) > 0
    disp('Error: entry in P cannot be larger than one!');
    P=[];
end

end

function HGL = Mnum_2_HGL(Mnum)
% Complete binary tree
% % Hierarchy Group Label
% %     H1  H2  H3
% HGL = [1   1   1 ; %N1
%        2   1   1 ; %N2
%        3   2   1 ; %N3
%        4   2   1 ];%N4
Hnum = log2(Mnum)+1; % number of hierarchies; so Mum must be 2^n
HGL = zeros(Mnum,Hnum);
for i = 1:Mnum
    for j = 1:Hnum
        HGL(i,j) = floor((i-1)/2^(j-1))+1;
    end
end
end

function CL = HGL_2_CL(HGL)

% Build the connection level matrix
Mnum = length(HGL(:,1));
CL = zeros(Mnum,Mnum);
for i_pre = 1:Mnum
    for j_post = 1:Mnum
        % find the lowest hierarchy where pre- and post-module share group label
        ind = find(HGL(i_pre,:) == HGL(j_post,:)); % note that (NaN != NaN) is true
        if ~isempty(ind)
            lowestH = ind(1);
            CL(i_pre,j_post) = lowestH;
        end
    end
end

end






