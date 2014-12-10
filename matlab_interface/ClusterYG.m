function Result_cell = ClusterYG( Result_cell, save_figure )
% assuming equal cluster size

disp('cluster_rate...');
tic;

% cluster_V_stat = 0;

% input check
if nargin <= 1 && nargout == 0
    save_figure = 1; % default
end

if nargin <= 1 && nargout > 0
    save_figure = -1; % default
end

if save_figure == 1
    figure_visibility = 'off'; % 'on', 'off'
else
    figure_visibility = 'on';
end


Result_num = length(Result_cell);
for r_num = 1:Result_num
    Result_cell{r_num} = cluster_sorted_rate( Result_cell{r_num} );
end


end




% % threshold for up-down analysis
% theta = 10; % Hz
% 
% 
% Result_num = length(Result_cell);
% 
% for r_num = 1:Result_num
%     
%     %%%%%%%%%%% cluster mean rate
%     Result_cell{r_num} = get_cluster_rate( Result_cell{r_num} );
%     
%     
%     %%%%%%%%%%%% cluster up and down state analysis
%     Result_cell{r_num}.up_down_analysis = cluster_up_down_state_analysis(Result_cell{r_num}, theta);
%     R = Result_cell{r_num};
%     
%     %%%%%%%%%%%% plot rate
%     if save_figure ~= -1 && ~isempty(Result_cell{r_num}.up_down_analysis.up_C_label)
%         
%         h_rate = figure('NumberTitle','Off','Name','cluster mean rate','units','normalized','position',[0 0 1 1], ...
%             'visible', figure_visibility, 'Color','w','PaperPositionMode', 'default');
%         subplot(5,1,1:4);hold on;
%         
%         edges = 1:R.ExplVar.Mnum;
%         count = histc(Result_cell{r_num}.up_down_analysis.up_C_label,edges);
%         bar(edges,count,'histc');
% %         T = (1:R.reduced.step_tot)*R.reduced.dt/1000;
% %         plot(T,C_rate); ylabel('Hz');xlabel('sec');
% %         legend(num2str(transpose(1:R.ExplVar.Mnum)));
% %         
% %         % show up-down dectection results
% %         if Result_cell{r_num}.up_down_analysis.up_overlapping == 0
% %             upA = Result_cell{r_num}.up_down_analysis.upA;
% %             upB = Result_cell{r_num}.up_down_analysis.upB;
% %             up_C_label = Result_cell{r_num}.up_down_analysis.up_C_label;
% %             for a = 1:length(upA)
% %                 plot(T(upA(a)),C_rate(up_C_label(a), upA(a)),'ko',T(upB(a)),C_rate(up_C_label(a), upB(a)),'ko');
% %             end
% %         end
%         
%         % Write comments
%         subplot(5,1,5, 'visible','off')
%         text(0.5, 0.5, R.comments, ...
%             'VerticalAlignment', 'top', ...
%             'HorizontalAlignment', 'center',...
%             'FontSize',10,'FontWeight','normal', 'interpreter', 'none'); % ...'interpreter', 'none'... to show underscore
%         
%         
%         % save figure
%         if save_figure == 1
%             fprintf('\t Saving figure...');
%             print(h_rate, '-dpsc2', strcat( R.stamp, '_cluster_rate'));
%             delete(h_rate);
%             fprintf('Saving done.\n');
%         else
%             next = input('\t Next figure?');
%             delete(h_rate);
%         end
%     end
%     
%     
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if cluster_V_stat == 1 && ~isempty(R.pop_V_sample)
%         
%         % cluster mean&std of membrane potential distribution
%         C_label = R.C_label;
%         % C_label = ceil((1:R.N(1))./round(R.N(1)/R.ExplVar.Mnum)); % cluster membership label
%         C_mean = zeros(R.ExplVar.Mnum, length(R.pop_V_sample{1}(1,:))); % 1st moment
%         C_std = zeros(R.ExplVar.Mnum, length(R.pop_V_sample{1}(1,:))); % 2nd moment
%         C_skewness = zeros(R.ExplVar.Mnum, length(R.pop_V_sample{1}(1,:))); % 3rd moment
%         C_kurtosis = zeros(R.ExplVar.Mnum, length(R.pop_V_sample{1}(1,:))); % 4th moment
%         
%         for cc = 1:R.ExplVar.Mnum
%             C_begin = find(C_label == cc, 1, 'first');
%             C_end = find(C_label == cc, 1, 'last');
%             C_mean(cc,:) = mean(R.pop_V_sample{1}(C_begin:C_end,:));
%             C_std(cc,:) = std(R.pop_V_sample{1}(C_begin:C_end,:));
%             C_skewness(cc,:) = skewness(R.pop_V_sample{1}(C_begin:C_end,:));
%             C_kurtosis(cc,:) = kurtosis(R.pop_V_sample{1}(C_begin:C_end,:));
%         end
%         
%         % plot mean&std
%         %T = (1:R.reduced.step_tot)*R.reduced.dt/1000;
%         h_pop_V = figure('NumberTitle','Off','Name','cluster potential distribution: 1st to 4th moment','units','normalized','position',[0 0 1 1], 'visible', figure_visibility, 'Color','w');
%         
%         %     linprop={'r','g','b','r','c','m','k','y'};
%         %     for i = 1:R.ExplVar.Mnum
%         %         shadedErrorBar([],C_mean(i,:),C_std(i,:), linprop{i}, 0.2);%plot(T,C_mean);
%         %         nextline = input('next line?');
%         %     end % this plot is not very clear
%         subplot(5,1,1);hold on;
%         plot(1:length(C_mean(1,:)),C_mean);ylabel('mean mV');set(gca, 'xtick', []);
%         subplot(5,1,2);hold on;
%         plot(1:length(C_mean(1,:)),C_std);ylabel('std mV');set(gca, 'xtick', []);
%         subplot(5,1,3);hold on;
%         plot(1:length(C_mean(1,:)),C_skewness);ylabel('skewness');set(gca, 'xtick', []);
%         subplot(5,1,4);hold on;
%         plot(1:length(C_mean(1,:)),C_kurtosis);ylabel('kurtosis ("peakedness") ');set(gca, 'xtick', []); %xlabel('sec');
%         
%         %legend(num2str(transpose(1:R.ExplVar.Mnum)));
%         
%         % Write comments
%         subplot(5,1,5, 'visible','off')
%         text(0.5, 0.5, R.comments, ...
%             'VerticalAlignment', 'top', ...
%             'HorizontalAlignment', 'center',...
%             'FontSize',10,'FontWeight','normal', 'interpreter', 'none'); % ...'interpreter', 'none'... to show underscore
%         
%         
%         % save figure
%         if save_figure == 1
%             fprintf('\t Saving figure...');
%             print(h_pop_V, '-dpdf', strcat( R.stamp, '_cluster_potential'));
%             delete(h_pop_V);
%             fprintf('Saving done.\n');
%         else
%             next = input('\t Next figure?');
%             delete(h_pop_V);
%         end
%         
%         
%     end
% end
% 
% toc;
% 
% end
% 
% 
% 
% 
% 
% 
% 
% 
% function result = cluster_up_down_state_analysis(R, theta)
% % A: cluster up-state beginning time vector index
% % B: cluster up-state endding time vector index
% 
%     dt = R.reduced.dt;
% 
%     % thresholding
%     result.theta = theta;
%     c_bin = R.C_rate > theta; % bin: binary
%     result.high_ratio = sum(sum(c_bin,1) > 0)/length(c_bin(1,:)); % ratio of high activity state
%     
%     % begining-ending points (A&B points) detection
%     Mnum = length(c_bin(:,1));
%     padding = false(Mnum,1);
%     beginning = [padding, c_bin(:,1:end-1) == 0 & c_bin(:,2:end) == 1];  % if ...01111..., detect first 1
%     ending = [c_bin(:,1:end-1) == 1 & c_bin(:,2:end) == 0, padding]; % if ...11110...., detect last 1
%     
%     
%     
%     
%     % up and down detection (for each cluster individually)
%     result.up_overlapping = 0;
%     AB_lumped = []; % populational up-state A&B points
%     AB_lumped_C_label = []; % cluster label for A&B points
%     AB_lumped_01 = []; % 0 for A point and 1 for B point
%     result.up_duration = [];
%     result.upA = [];
%     result.upB = [];
%     result.up_C_label = [];
%     result.down_duration = [];
%     
%     for i = 1:Mnum
%         A_temp = find(beginning(i,:)); % beginning point of up state
%         A_temp = A_temp(:)'; % row vector
%         B_temp = find(ending(i,:)); % ending point of up state
%         B_temp = B_temp(:)'; % row vector
%         AB_lumped = [AB_lumped, A_temp, B_temp]; % addding cluster up-states into populational ones
%         AB_lumped_C_label = [AB_lumped_C_label, ones(1,length(A_temp))*i, ones(1,length(B_temp))*i];
%         AB_lumped_01 = [AB_lumped_01, zeros(1,length(A_temp)), ones(1,length(B_temp))];% overlapping detection and warning
%         [AB_lumped, AB_ind] = sort(AB_lumped);             % A B a b, A a B b, A a b B,
%         AB_lumped_01 = AB_lumped_01(AB_ind);   % 0 1 0 1, 0 0 1 1, 0 0 1 1, the last two cases are overlapping
%         AB_lumped_C_label = AB_lumped_C_label(AB_ind);
%         
%         cluster_label = ones(size(A_temp))*i;
%         [up_duration, upA, upB, up_C_label] = duration_from_AB(A_temp,B_temp,dt, cluster_label);
%         [down_duration] = duration_from_AB(B_temp,A_temp,dt, cluster_label );
%         
%         result.up_duration    = [ result.up_duration up_duration(:)' ];
%         result.upA            = [ result.upA upA(:)' ];
%         result.upB            = [ result.upB upB(:)' ];
%         result.up_C_label     = [ result.up_C_label up_C_label(:)' ];
%         result.down_duration  = [ result.down_duration down_duration(:)' ];
%     end
%     
%     if ~isempty(AB_lumped)
%         if nnz(AB_lumped_01(1:end-1) == 0 & AB_lumped_01(2:end) == 0)
%             warning('overlapping in cluster rate up-states!');
%             result.up_overlapping = 1;
%         end
%     end
%         
%         
%     
% %     % up and down detection % old version, requires non-overlapping
% %     result.up_overlapping = 0;
% %     AB_lumped = []; % populational up-state A&B points
% %     AB_lumped_C_label = []; % cluster label for A&B points
% %     AB_lumped_01 = []; % 0 for A point and 1 for B point
% %     for i = 1:Mnum
% %         A_temp = find(beginning(i,:)); % beginning point of up state
% %         A_temp = A_temp(:)'; % row vector
% %         B_temp = find(ending(i,:)); % ending point of up state
% %         B_temp = B_temp(:)'; % row vector
% %         AB_lumped = [AB_lumped, A_temp, B_temp]; % addding cluster up-states into populational ones
% %         AB_lumped_C_label = [AB_lumped_C_label, ones(1,length(A_temp))*i, ones(1,length(B_temp))*i];
% %         AB_lumped_01 = [AB_lumped_01, zeros(1,length(A_temp)), ones(1,length(B_temp))];% overlapping detection and warning
% %         [AB_lumped, AB_ind] = sort(AB_lumped);             % A B a b, A a B b, A a b B,
% %         AB_lumped_01 = AB_lumped_01(AB_ind);   % 0 1 0 1, 0 0 1 1, 0 0 1 1, the last two cases are overlapping
% %         AB_lumped_C_label = AB_lumped_C_label(AB_ind);
% %         if ~isempty(AB_lumped)
% %             if nnz(AB_lumped_01(1:end-1) == 0 & AB_lumped_01(2:end) == 0)
% %                 warning('overlapping in cluster rate up-states!');
% %                 result.up_overlapping = 1;
% %                 break;
% %             end
% %         end
% %     end
% %     if result.up_overlapping == 0
% %         A = AB_lumped(AB_lumped_01==0);
% %         B = AB_lumped(AB_lumped_01==1);
% %         A_label = AB_lumped_C_label(AB_lumped_01==0);
% %         % up state duration
% %         [result.up_duration, result.upA, result.upB, result.up_C_label] = duration_from_AB(A,B,dt,A_label);
% %         % down state duration
% %         % Note that A point of up state is B point of down state and vice
% %         % versa. However, for down state, only duration is meaningful.
% %         [result.down_duration] = duration_from_AB(B,A,dt);
% %     else
% %         result.up_duration = [];
% %         result.upA = [];
% %         result.upB = [];
% %         result.up_C_label = [];
% %         result.down_duration = [];
% %     end
% 
% 
% 
% end
% 
% 
% 
% 
% function [duration, A_pair, B_pair, A_C_label] = duration_from_AB(A,B,dt,A_C_label)
% if nargin == 3
%     A_C_label = zeros(1,length(A));
% end
% 
% duration = [];
% if ~isempty(A) && ~isempty(B)
%     % pairing of point A's and point B's
%     A_pair = A;
%     B_pair = B;
%     if length(B_pair) > length(A_pair) 
%         % B A-B A-B A-B
%         B_pair(1) = [];
%     elseif length(B_pair) < length(A_pair) 
%         % A-B A-B A-B A 
%         A_pair(end) = [];
%         A_C_label(end) = [];
%     elseif B_pair(1) < A_pair(1) && length(B_pair) == length(A_pair) % the latter is always true but kept for readability
%         % B A-B A-B A-B A
%         B_pair(1) = [];
%         A_pair(end) = [];
%         A_C_label(end) = [];
%     end
%     % calculate durations
%     if ~isempty(A_pair) && ~isempty(B_pair)
%         duration = (B_pair - A_pair)*dt;
%     end
% else
%     disp('A and B points for up state are empty!');
%     A_pair = [];
%     B_pair = [];
%     A_C_label = [];
% end
% 
% end
