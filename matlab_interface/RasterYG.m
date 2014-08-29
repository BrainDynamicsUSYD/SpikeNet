function RasterYG( Result_cell, save_figure )
%UNTITLED Summary of this function goes here
%   If the firing history is too long, data will be segmented into several
%   history fractions and plotted separately.
disp('RasterYG...');
tic;

sample_size = 500; % sample neurons for raster plot
seg_size = 4*10^4; % 2*10^4 for 2-pop, segmentation size for each plot

if nargin == 0
    Result_cell = CollectRYG();
end

if nargin <= 1
    save_figure = 1; % default
end

if save_figure == 1
    figure_visibility = 'off'; % 'on', 'off'
else
    figure_visibility = 'on';
end

Result_num = length(Result_cell);

for r_num = 1:Result_num
    R_temp = Result_cell{r_num};
    
    % Dump fields
    reduced_dt = R_temp.reduced_dt;
    reduced_step_tot = R_temp.reduced_step_tot;
    N = R_temp.N;
    comments = R_temp.comments;
    rate_sorted = R_temp.Analysis.rate_sorted;
    
    % Num_pop = R_temp.Num_pop; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Num_pop = 1;
    
    % Segmetation
    seg_num = ceil(reduced_step_tot/seg_size);
    for seg = 1:seg_num
        if seg < seg_num
            seg_ind = ((seg-1)*seg_size+1):(seg*seg_size);
        else
            seg_ind = ((seg-1)*seg_size+1):(reduced_step_tot);
        end
        % Dump fields
        
        reduced_num_spikes = cell(1,Num_pop);
        reduced_num_ref = cell(1,Num_pop);
        reduced_spike_hist = cell(1,Num_pop);
        for pop_ind = 1:Num_pop
            reduced_num_spikes{pop_ind} = R_temp.reduced_num_spikes{pop_ind}(seg_ind);
            reduced_num_ref{pop_ind} = R_temp.reduced_num_ref{pop_ind}(seg_ind);
            reduced_spike_hist{pop_ind} = R_temp.reduced_spike_hist{pop_ind}(:,seg_ind);
        end
        reduced_T = seg_ind*reduced_dt;
        
        %-------------------------------------------------------------------------%
        % Plot
        h_raster = figure('NumberTitle','Off','Name','Raster plot','units','normalized','position',[0 0 1 1], 'visible', figure_visibility, 'Color','w');
        
        axes_matrix = zeros(3,Num_pop);
        for pop_ind = 1:Num_pop
            %%%%%%%%%% Plot raster plot
            axes_matrix(1,pop_ind) = subplot(7,Num_pop,(0:3)*Num_pop+pop_ind);hold on;
            if nnz(reduced_num_spikes{pop_ind}) > 0
                % down-sampling
                if N(pop_ind) >= sample_size
                    ind_sample = sort(randperm(N(pop_ind),sample_size));
                else
                    ind_sample = 1:1:N(pop_ind);
                end
                %             % additional color coding for clustered structure
                %             if pop_ind == 1 && any( strcmp(fieldnames(R_temp.ExplVar), 'Mnum') ) % if clustered/hierarchical structure
                %                 Mnum = R_temp.ExplVar.Mnum;
                %                 color_vec = {'b','g','r','c','m','y','k'}; % following matlab default color coding, the colors are not good for raster plot!!!!
                %                 C_label = ceil((1:N(1))./(N(1)/Mnum)); % cluster membership label
                %                 sample_C_label = C_label(ind_sample);
                %                 shift_Y = 0;
                %                 for cc = 1:Mnum
                %                     ind_sample_cc = ind_sample(sample_C_label == cc);
                %                     [Y0,X,~] = find(reduced_spike_hist{pop_ind}(ind_sample_cc,:));
                %                     Y = Y0 + shift_Y;
                %                     shift_Y = shift_Y + length(ind_sample_cc);
                %                     line([X(:)'; X(:)']*reduced_dt,[Y(:)'-1;Y(:)'],'Color',color_vec{mod(cc-1,7)+1});
                %                 end
                %             else % no color coding
                %                 [Y,X,~] = find(reduced_spike_hist{pop_ind}(ind_sample,:));
                %                 line([X(:)'; X(:)']*reduced_dt,[Y(:)'-1;Y(:)'],'Color','k');
                %             end
                [Y,X,~] = find(reduced_spike_hist{pop_ind}(ind_sample,:));
                line(([X(:)'; X(:)']+seg_ind(1)-1)*reduced_dt,[Y(:)'-1;Y(:)'],'Color','k');
                ylim([0,length(ind_sample)]);
                if pop_ind == 1
                    if rate_sorted == 1
                        ylabel('Rate-sorted sample neuron index');
                    else
                        ylabel('Sample neuron index');
                    end
                end
                set(gca, 'xtick', []);
                xlim([reduced_T(1), reduced_T(1)+(length(seg_ind)-1)*reduced_dt]); % make sure all the plots have the same axis scale
            end
            
            %%%%%%%%%%  Plot number of spikes
            axes_matrix(2,pop_ind) = subplot(7,Num_pop,4*Num_pop+pop_ind);hold on;
            line([reduced_T; reduced_T], [zeros(1, length(reduced_T)); reduced_num_spikes{pop_ind}/N(pop_ind)*100], 'Color', [255 30 30]/255);
            if pop_ind == 1
                ylabel('% Firing');
            end
            set(gca, 'xtick', []);
            
            %%%%%%%%%%  Plot number of refractory neurons
            axes_matrix(3,pop_ind) = subplot(7,Num_pop,5*Num_pop+pop_ind);hold on;
            line([reduced_T; reduced_T], [zeros(1, length(reduced_T)); reduced_num_ref{pop_ind}/N(pop_ind)*100], 'Color', [35 163 200]/255);
            if pop_ind == 1
                ylabel('% Refractory');
            end
            xlabel(['t (ms), pop No.', num2str(pop_ind)]);
            
        end
        
        % Link axes to synchronise them when zooming
        for pop_ind = 1:Num_pop
            linkaxes(axes_matrix(1:3,pop_ind),'x');
        end
        
        % Keep tick lables while remove tick marks
        for pop_ind = 1:Num_pop
            for i = 1:3
                set(axes_matrix(i,pop_ind), 'Ticklength', [0 0]);
            end
        end
        
        % Set initial y-axes to be same, since only one property can be linked
        for i = 2:3
            ylimData = [inf -inf]; % initialise ylim
            % find common ylim
            for pop_ind = 1:Num_pop
                ylimData_new = ylim(axes_matrix(i,pop_ind));
                ylimData(1) = min(ylimData(1),ylimData_new(1));
                ylimData(2) = max(ylimData(2),ylimData_new(2));
            end
            % set common ylim
            for pop_ind = 1:Num_pop
                ylim(axes_matrix(i,pop_ind),ylimData);
            end
        end
        
        
        % Write comments
        subplot(7,Num_pop,6*Num_pop+(1:Num_pop), 'visible','off')
        text(0.5, 0.5, comments, ...
            'VerticalAlignment', 'top', ...
            'HorizontalAlignment', 'center',...
            'FontSize',10,'FontWeight','normal', 'interpreter', 'none'); % ...'interpreter', 'none'... to show underscore
        
        
        % save figure
        if save_figure == 1
            fprintf('\t Saving figure...');
            
            % Bad solutions
            % saveas(h_raster,strcat('data/', R_temp.name{1}, '_raster'), 'fig'); % use openfig('filename.fig','new','visible') to show
            % imwrite(zbuffer_cdata(h_raster),strcat('data/', R_temp.name{1}, '_raster.jpg')); %, 'Quality', 100);
            % saveas(h_raster,strcat('data/', R_temp.name{1}, '_raster'), 'jpg');
            
            % Good solution
            % set(h_raster,'PaperPosition',[0.25 0.25 6 4]); % default PaperUnits is inches, defualt PaperOrientation is portrait
            print(h_raster, '-dpsc2', strcat( R_temp.stamp, '_raster_',num2str(seg)));
            % Generating PostScript file only demands a memory proportional to
            % the information content of the plot in instead of to the
            % unnecessarily large amount of pixels when MATLAB generating .jpg
            % file.
            % Also, rendering PostScipt file is fast enough for a local
            % machine and it's free of concern about resolution.
            % In short, vector image is much better in this case!
            
            delete(h_raster);
            fprintf('done.\n');
        else
            next = input('\t Next figure?');
            delete(h_raster);
        end
        
    end
end

toc;

end

