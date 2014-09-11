function [S]=SwitchCount(C_label,loop_num,r_num, Xpara)
%r_num: realization number

plot_fig = 1;

loop_max = max(loop_num);
S{1} = zeros(1,loop_max);
S{2} = zeros(1,loop_max);
S{3} = zeros(1,loop_max);
S{4} = zeros(1,loop_max);
Up_tot = zeros(1,loop_max);
for i = 1:loop_max
    %%%% Total
    Up_tot(i) = length(C_label(loop_num == i));
    %%%% L4
    % (1,2,3,4) <-> (5,6,7,8)
    C_temp = C_label(loop_num == i);
    S{4}(i) = sum(FindSwitching(C_temp, 4));
    
    %%%% L3
    % (1,2) <-> (3,4)
    C_temp = C_label(loop_num == i);
    C_temp(C_temp>4) = NaN;
    S{3}(i) = S{2}(i) + sum(FindSwitching(C_temp, 2));
    
    % (5,6) <-> (7,8)
    C_temp = C_label(loop_num == i);
    C_temp(C_temp<=4) = NaN;
    S{3}(i) = S{2}(i) + sum(FindSwitching(C_temp, 6));
    
    %%%% L2
    % 1 <-> 2
    C_temp = C_label(loop_num == i);
    C_temp(C_temp>2) = NaN;
    S{2}(i) = S{1}(i) + sum(FindSwitching(C_temp, 1));
    
    % 3 <-> 4
    C_temp = C_label(loop_num == i);
    C_temp(C_temp<=2) = NaN;
    C_temp(C_temp>4) = NaN;
    S{2}(i) = S{1}(i) + sum(FindSwitching(C_temp, 3));
    
    % 5 <-> 6
    C_temp = C_label(loop_num == i);
    C_temp(C_temp<=4) = NaN;
    C_temp(C_temp>6) = NaN;
    S{2}(i) = S{1}(i) + sum(FindSwitching(C_temp, 5));
    
    % 7 <-> 8
    C_temp = C_label(loop_num == i);
    C_temp(C_temp<=6) = NaN;
    S{2}(i) = S{1}(i) + sum(FindSwitching(C_temp, 7));
    
    %%%% L1
    % i <-> i
    C_temp = C_label(loop_num == i);
    S{1}(i) = sum([false, C_temp(1:end-1) == C_temp(2:end)]);
        
    
    
end

for i = 1:4
    S{i} = transpose(vec2mat(S{i},r_num));
end

Up_tot = transpose(vec2mat(Up_tot,r_num));


% plot results
if plot_fig == 1
%     Labels = {'r--','b:','k-.'};
%     figure(1);
%     subplot(4,1,1:3);hold on;
%     
%     step_para = Xpara(2)-Xpara(1);
%     xlim([Xpara(1)-step_para,Xpara(end)+step_para]);
%     % xlabel('Relative L4-connection probability');
%     ylabel('Switching frequency');
%     for i = 1:3
%         errorbar(Xpara,mean(S{i}),std(S{i}),Labels{i});
%     end
%     legend('L2','L3','L4');
%  
%     subplot(4,1,4);
%     errorbar(Xpara,mean(Up_tot),std(Up_tot));
%     xlim([Xpara(1)-step_para,Xpara(end)+step_para]);
%     xlabel('Relative connection probability');
%     ylabel('Up-state count');
%     



% Showing all the data points, not good enough
    Labels = {'r.','b.','k.'};
    figure(1);hold on;
    sp = Xpara(2)-Xpara(1);
    spp = sp/3;
    for j = 1:r_num
        for i = 1:length(Xpara)
            for k = 1:3
                plot(Xpara(i)+(k-1)*spp, S{k}(j,i), Labels{k},'MarkerSize',30);
            end
        end
    end
    
end



    function S = FindSwitching(V,threshold)
        S1 =  [false, V(1:end-1) <= threshold  & V(2:end) > threshold];
        S2 =  [false, V(1:end-1) > threshold  & V(2:end) <= threshold];
        S = S1 | S2;
    end





end