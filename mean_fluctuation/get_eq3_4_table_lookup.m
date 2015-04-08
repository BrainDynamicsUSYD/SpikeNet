function [ out ] = get_eq3_4_table_lookup(miu_V, sigma_V, table,  arg)
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


v = interp2(table.miu_V_mat, table.sigma_V_mat, table.v_mat, miu_V, sigma_V);
CV = interp2(table.miu_V_mat, table.sigma_V_mat, table.CV_mat, miu_V, sigma_V);


if strcmp(arg,'v')
    out = v;
elseif strcmp(arg,'CV')
    out = CV;
elseif strcmp(arg,'both')
    out = [v, CV];
elseif strcmp(arg,'v*CV^2')
    out = v*CV^2;
    
end

end