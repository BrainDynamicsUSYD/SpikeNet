% built the table of [v, CV] for quick look-up

warning off;
i =0;

% for miu_V = 0.5:0.5:30
%     i = i + 1;
%     j = 0;
%     for sigma_V = 0.5:0.5:20
%         j = j+1;
%         out = get_eq3_4(miu_V, sigma_V, 'both');
%         v = out(1);
%         CV = out(2);
% 
%         miu_V_mat(i,j) = miu_V;
%         sigma_V_mat(i,j) = sigma_V;
%         v_mat(i,j) = v;
%         CV_mat(i,j) = CV;
%         
%         fprintf('%d,%d;',i,j)
%     end
%     
%     save('v_CV_tables.mat', 'miu_V_mat', 'sigma_V_mat', 'v_mat', 'CV_mat')
% end


% for miu_V = 0.5:0.5:30
%     i = i + 1;
%     j = 0;
%     for sigma_V = 0.3:0.05:1
%         j = j+1;
%         out = get_eq3_4(miu_V, sigma_V, 'both');
%         v = out(1);
%         CV = out(2);
% 
%         miu_V_mat(i,j) = miu_V;
%         sigma_V_mat(i,j) = sigma_V;
%         v_mat(i,j) = v;
%         CV_mat(i,j) = CV;
%         
%         fprintf('%d,%d;',i,j)
%     end
%     
%     save('v_CV_tables_aug.mat', 'miu_V_mat', 'sigma_V_mat', 'v_mat', 'CV_mat')
% end

for miu_V = 17:0.2:22
    i = i + 1;
    j = 0;
    for sigma_V = 0.5:0.025:0.8
        j = j+1;
        out = get_eq3_4(miu_V, sigma_V, 'both');
        v = out(1);
        CV = out(2);

        miu_V_mat(i,j) = miu_V;
        sigma_V_mat(i,j) = sigma_V;
        v_mat(i,j) = v;
        CV_mat(i,j) = CV;
        
        fprintf('%d,%d;',i,j)
    end
    
    save('v_CV_tables_fig4.mat', 'miu_V_mat', 'sigma_V_mat', 'v_mat', 'CV_mat')
end


% %%%%%%%%%%%% table should be loaded as this
% table = load('v_CV_tables.mat');
% % manually fix some (possibly) numerical errors 
% table.CV_mat(1:13,1) = 1;
% table.CV_mat(end-2:end,1) = 0;
% table.v_mat(end-2:end,1) = 0.5;
% % something is still wrong here!!!
% 
% % transpose for interp2
% table.miu_V_mat = table.miu_V_mat'; 
% table.sigma_V_mat = table.sigma_V_mat'; 
% table.CV_mat = table.CV_mat'; 
% table.v_mat = table.v_mat';

% % show plots
% if show_figure == 1
% figure(1);
% set(gcf,'color','w');
% subplot(1,2,2);
% surf(table.miu_V_mat,table.sigma_V_mat,table.CV_mat );
% xlabel('\mu_V');
% ylabel('\sigma_V');
% zlabel('CV');
% 
% subplot(1,2,1);
% surf(table.miu_V_mat,table.sigma_V_mat,table.v_mat*1000 );
% xlabel('\mu_V');
% ylabel('\sigma_V');
% zlabel('Firing Rate (Hz)');
% end
%%%%%%%%%%%% 

