clc;close all;

% Input check and default values

seg_size = 4*10^4; % 2*10^4 for 2-pop, segmentation size for each plot
seg = 1;

% Dump fields
dt = R_temp.dt;
step_tot = R_temp.step_tot;

% Segmetation
seg_num = ceil(step_tot/seg_size);
if seg < seg_num
    seg_ind = ((seg-1)*seg_size+1):(seg*seg_size);
else
    seg_ind = ((seg-1)*seg_size+1):(step_tot);
end

T = seg_ind*dt;

rate = R_temp.C_rate(:,seg_ind);

% T = length(rate(1,:));
% 
% figure(1);
% 
% cutoff = 1.2;
% D = 3;
% data = [];
% for i = 1:100:T
%     d = pdist(rate(:,i));
%     
%     ds = squareform(d);
%     l = linkage(d);
%     
%     %I = inconsistent(l);
%     
%     T = cluster(l,'cutoff', cutoff, 'Depth',D)'
% 
%     %if sum(T==1) == 8
%         p = data(:,i)'
%         I = inconsistent(l,D)
%         % x=input('next?');
%          clc
%     %end
%     
%     dendrogram(l);
%     set(gca, 'xlim', [0,12]);
%     pause(0.5);
%     
% end
% 
% plot(1:length(data),data);

myylim = [-5,5];

figure(1)
for i = 1:4
    subplot(4,1,i);
    i
    [ CC ] = CorrCoefYG( rate(i*2-1,:), rate(i*2,:), 0, 1);
    mean(CC)
    plot(T,CC);
    ylim(myylim);
end


figure(2)
for i = 1:2
    subplot(2,1,i);
    i
    [ CC ] = CorrCoefYG( mean(rate(i*4-3:i*4-2,:)), mean(rate(i*4-1:i*4,:)), 0, 1);
    mean(CC)
    plot(T,CC);
    ylim(myylim);
end


figure(3);
[ CC ] = CorrCoefYG( mean(rate(1:4,:)), mean(rate(5:6,:)), 0, 1);
mean(CC)
plot(T,CC);
ylim(myylim);


