function VI_SamplePlot( Result_cell )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here




Result_num = length(Result_cell);


figure_visibility = 'on';

for r_num = 1:Result_num
    if  any( strcmp(fieldnames(Result_cell{r_num}), 'VI_sample') );
        R_temp = Result_cell{r_num};
        % Dump fields
        dt = R_temp.dt;
        step_tot = R_temp.step_tot;
        spike_hist = R_temp.spike_hist;
        VI_sample = R_temp.VI_sample;
        PopPara = R_temp.PopPara;
        comments = R_temp.comments;
        
        T = (1:step_tot)*dt;
        
        % Gaussian filter
        sigma_gaussian = 50; % ms, which is width???
        kernel = spike_train_kernel_YG(sigma_gaussian, dt, 'gaussian');
        
        
        % Plot VI
        sample_size = length(VI_sample.ind(:,1));
        for i = 1:sample_size
            
            ind = VI_sample.ind(i,2);
            pop = VI_sample.ind(i,1);
            figure_name = sprintf('Activity of sample cell No.%g of Pop %g',ind, pop);
            h_VI = figure('NumberTitle','off','Name', figure_name, 'units','normalized','position',[0 0 1 1], 'visible', figure_visibility, 'Color','w');
            
            axes_matrix = zeros(3, 1);
            
            
            % V plot
            axes_matrix(1) = subplot(6, 8, 1:7 );hold on;
            V_th = PopPara.V_th(pop);
            plot(T, ShowSpike(VI_sample.V(i,:),V_th));
            plot(T, V_th*ones(size(T)), 'r--');
            ylim([  min(-70, min(VI_sample.V(i,:))) V_th+10  ]); % [reset spike_peak]
            set(gca,'xtick',[],'box','off')
            
            % V distribution bar chart
            subplot(6,8,8);hold on;
            bin_num = 100;
            tau_ref = PopPara.tau_ref(pop);
            ref_steps = round(tau_ref/dt);
            V_remove = RemoveSpikeRef(VI_sample.V(i,:), V_th, ref_steps); % remove spike and refractory data
            bin_edge = linspace(min(V_remove), max(V_remove), bin_num+1);
            bin_count = histc(V_remove,bin_edge);
            barh(bin_edge, bin_count, 'histc');
            ylim(ylim(axes_matrix(1)));
            set(gca,'xtick',[],'box','off')
            
%             % multitaper estimation of power spectrum
%             pmtm(V_remove);set(gca,'xscale','log');
            
            
            % rate plot
            axes_matrix(2) = subplot(6, 8, 9:15);hold on;grid on;
            rate = SpikeTrainConvolve(spike_hist{pop}(ind, :), kernel);
            plot(T, rate);
            set(gca,'xtick',[],'box','off')
            
            
            % current plot
            axes_matrix(3) = subplot(6, 8, [17:23, 25:31, 33:39]); hold on;grid on;
            marker = {'r','g','b','y','k','c','m'}; % for 6 different currents
            fname_cell = fieldnames(VI_sample);
            for f = 1:length(fname_cell)-2 % skip first two fields
                fn = fname_cell{f+2};
                if nnz(VI_sample.(fn)(i,:)) > 0
                    plot(T,VI_sample.(fn)(i,:),marker{f});
                else
                    plot(0,0,'w');
                end
            end
            legend('leak','AMPA','GABA','NMDA','GJ','ext');
            set(gca,'box','off')
            
            
            
            
            % Link axes to synchronise them when zooming
            linkaxes(axes_matrix(1:3),'x');
            
            % Keep tick lables while remove tick marks
            for f = 1:3
                set(axes_matrix(f), 'Ticklength', [0 0]);
            end
            
            
            % Comments (Explanatory variables)
            subplot(6,1,6, 'visible','off');
            text(0.5, 0.5, comments, ...
                'VerticalAlignment', 'top', ...
                'HorizontalAlignment', 'center',...
                'FontSize',10,'FontWeight','normal', 'interpreter', 'none'); % ...'interpreter', 'none'... to show underscore
       
            % next figure
            next = input('\t Next figure?');
            delete(h_VI);
        end
        
        
       
    end
end

end

function V = ShowSpike(V, theta)
peak = 0;
V([(V(1:end-1) < theta) & (V(2:end) >= theta), false]) = peak;
end


function V = RemoveSpikeRef(V, theta, ref_steps)
% logical vector of spiking event
spike = [(V(1:end-1) < theta) & (V(2:end) >= theta), false]; 
% logical vector of being refractory
ref = spike; 
for i = 1:ref_steps
    ref = ref | [false(1,i) spike(1:end-i)];
end
% delete refractoriness from V
V(ref) = [];
end


