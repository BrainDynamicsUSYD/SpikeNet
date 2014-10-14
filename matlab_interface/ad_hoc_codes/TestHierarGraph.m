% Hierarchical graph test
clc;clear;close all;

N = 4000;
Mnum = 8;
P0 = 0.1;
r = 0.6;


load('/import/yossarian1/yifan/Project1/data/1-20140709-1144-13158A11.mat');
A = A11;

%[ A ] =  MyHierarchyGraph('N', N, 'Mnum', Mnum, 'P0', P0, 'r', r);
figure(1);
spy(A);

C_label = ceil((1:N)./round(N/Mnum)); % cluster membership label
for i = 1:Mnum
    cluster = find(C_label == i);
    cc(i) = mean(clustering_coef_wd((A(cluster, cluster))));
end
cc

% % Clustering coefficients
% figure(2);
% color_vec = {'b','g','r','c','m','y','k','b--','g--','r--','c--','m--','y--','k--'};
% for i = 1:Mnum
%     cluster = find(C_label == i);
%     subplot(2,1,1);hold on;
%     cc = clustering_coef_wd((A(cluster, cluster))); % (N-by-1) vector, Fagiolo (2007) as well
%     plot(1:length(cluster),sort(cc),'Color',color_vec{mod(i-1,14)+1});
%     subplot(2,1,2);hold on;
%     [in_degree,~] = degrees_dir(A(cluster, cluster));
%     plot(1:length(cluster),sort(in_degree),'Color',color_vec{mod(i-1,14)+1});
% end