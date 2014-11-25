% ad hoc
% plot example dynamics stats


% clc;clear;close all;
%load 010-201409041850-21690_1409844753508_RYG.mat;
%load 008-201409041850-21044_1409846182212_RYG.mat; % similar results

text_fontsize = 12;
anno_fontsize = 16;


% Dump fields
dt = R_temp.dt;
step_tot = R_temp.step_tot;
neuron_sample = R_temp.neuron_sample;

h_EI = figure('NumberTitle','Off','Name','E-I balance plot','units','pixels',...
    'position',[ 645   308   677   602], ...
            'Color','w', 'PaperPositionMode', 'auto');



%% plot E-I current correlation
subaxis(2,2,1,'PB',0.1,'PL',0.05);hold on;
XCC = [];
for pop = 1:2
    n = length(neuron_sample.neuron_ind{pop});
    
    for i = 1:n
        x = neuron_sample.I_AMPA{pop}(i,:);
        y = -(neuron_sample.I_GABA{pop}(i,:));
        x = x(1:10:end); % down-sampling
        y = y(1:10:end);
        [xcc,lags] = crosscorr(x,y,round(100/dt)); % do NOT use xcorr!!
        %plot(lags*dt,xcc);
        
        if i == 1 && pop == 1
            XCC = xcc;
        else
            XCC = [XCC; xcc];
        end
        
    end
   
    %set(gca,'xscale','log');
end

xcc_mean = mean(XCC);
plot(lags*dt,xcc_mean);
xlabel('lag (ms)','fontsize',text_fontsize);
ylabel('cross-correlation','fontsize',text_fontsize);

anno_shift1 = [-0.09,0.05];
add_figure_label( 'A', anno_shift1, anno_fontsize );

%% plot histogram of corrcoef
subaxis(2,2,2,'PB',0.1,'PL',0.05);
corr_hist_plot(R_temp, 1);
xlabel('correlation','fontsize',text_fontsize);
ylabel('distribution','fontsize',text_fontsize)
% set(gca,'ytick',[]);
anno_shift2 = [-0.09,0.05];
add_figure_label( 'B', anno_shift2, anno_fontsize );




%%
subaxis(2,2,3,'PB',0.1,'PL',0.05);
avalanche = avalanche_detect(R_temp, 1);
ylabel('distribution','fontsize',text_fontsize)
xlabel('avalanche size','fontsize',text_fontsize);
anno_shift2 = [-0.09,0.05];
add_figure_label( 'C', anno_shift2, anno_fontsize );


% %%
% subaxis(2,2,4,'PB',0.1,'PL',0.05);
% pop_ISI_dist(R_temp, 1);
% 
% xlabel('correlation','fontsize',text_fontsize);
% ylabel('ISI (ms)','fontsize',text_fontsize)
% set(gca,'ytick',[]);
% anno_shift2 = [-0.07,0.05];
% add_figure_label( 'D', anno_shift2, anno_fontsize );



%% Results not good!!!
subaxis(2,2,4,'PB',0.1,'PL',0.05);
firing_rate_hist(R_temp, 1);

xlabel('firing rate (Hz)','fontsize',text_fontsize);
ylabel('distribution','fontsize',text_fontsize)
%set(gca,'ytick',[]);
anno_shift2 = [-0.09,0.05];
add_figure_label( 'D', anno_shift2, anno_fontsize );


% %%
% subaxis(2,3,5,'PB',0.1,'PL',0.05);
% spike_count_hist(R_temp, pop)
% 
% xlabel('spike_count','fontsize',text_fontsize);
% ylabel('distribution','fontsize',text_fontsize)
% %set(gca,'ytick',[]);
% anno_shift2 = [-0.09,0.05];
% add_figure_label( 'D', anno_shift2, anno_fontsize );


%%
set(gcf, 'InvertHardCopy', 'off'); % prevent the x-axis line to reappear when printed
fig_name = 'example_stats';



% print(h_EI, '-dpdf', fig_name );





%         
% %% plot Vand AMPA+GABA example
% seg = 1;
% pop_ind = 1;
% sample_ind = 6;
% 
% subaxis(2,2,1);
% neuron_V_plot(R_temp, pop_ind, sample_ind, seg);
%  
% subaxis(2,2,3);
% current_type = {'I_AMPA','I_GABA'};
% marker = {'r','g','b','y','k','c','m'}; % for 6 different currents
% plot_abs = 1; % plot abs() of currents
% plot_smooth = 1; % smooth curves 
% kernel_width = 10; % ms
% seg_size = 4*10^4; % 2*10^4 for 2-pop, segmentation size for each plot
% 
% % Segmetation
% seg_ind = get_seg(step_tot, seg_size, seg);
% T = seg_ind*dt;
% % plot
% hold on;
% for type = 1:length(current_type)
%         I = neuron_sample.(current_type{type}){pop_ind}(sample_ind,seg_ind);
%         if plot_abs == 1
%             I = abs(I);
%         end
%         if plot_smooth == 1
%              I = SpikeTrainConvolve(I, spike_train_kernel_YG( kernel_width, dt, 'Gaussian' ));
%         end
%         cut = 1:length(T);
%         cut(1:(kernel_width/dt)) = [];
%         cut(end-(kernel_width/dt):end) = [];
%         plot(T(cut),I(cut),marker{type});
%         xlim([min(T),max(T)]);
% end
% % legend('I_E','I_I');
% % set(lh,'Interpreter','none'); % plain text instead of Tex
% set(gca,'box','off', 'TickDir','out')
% scalebar([200 0], {'200 msec ','a'})
