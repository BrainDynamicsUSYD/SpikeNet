function [ Result_cell ] = AnalyseYG( Result_cell )
%
%   Detailed explanation goes here
disp('AnalyseYG...');
tic;


% choices
re_arrange = 0;
get_corrcoef = 0;


Result_num = length(Result_cell);
for r_num = 1:Result_num
    
    
    % Get Firing rate
    Result_cell{r_num} = get_spike_count(Result_cell{r_num});
    
    % Get ISI distribution
    Result_cell{r_num} = get_ISI(Result_cell{r_num});
    
    
    % Re-arrange the neuron index
    % >> according to firing rate
    Result_cell{r_num} = rate_sort(Result_cell{r_num}, re_arrange);
    
    
    % Irregularity measures for the network based on CV_ISI
    Result_cell{r_num} = get_CV_ISI(Result_cell{r_num});
    
    % Correlation coefficient (zero time-lag)
    if get_corrcoef == 1
        Result_cell{r_num} = get_CC(Result_cell{r_num});
    end
    
    
    % Generate comments
    Result_cell{r_num} = generate_comments(Result_cell{r_num});
    
end
toc;
end

function [R] = get_spike_count(R)
fprintf('\t Getting spike count...\n');
% dump fields
dt = R.dt;
step_tot = R.step_tot;
N = R.N;
Num_pop = R.Num_pop;
num_spikes = R.num_spikes;
spike_hist = R.spike_hist;
%
spike_count = cell(Num_pop,1); % Initialise individual neuron average firing rate
Hz_overall = 0; % Network-wide average firing rate in Hz over simulation
for pop_ind = 1:Num_pop
    if nnz(num_spikes{pop_ind}) > 0
        spike_count{pop_ind} = full(sum(spike_hist{pop_ind},2));
        Hz_overall = Hz_overall + sum(spike_count{pop_ind})/(sum(N)*dt/1000*step_tot);
    end
end
% Record results
R.Analysis.Hz_overall = Hz_overall;
R.Analysis.spike_count = spike_count;
end

function [R] = get_ISI(R)
% dump fields
dt = R.dt;
N = R.N;
Num_pop = R.Num_pop;
num_spikes = R.num_spikes;
spike_hist = R.spike_hist;
spike_count = R.Analysis.spike_count;
%
fprintf('\t Getting ISI distribution...\n');
ISI_dist = cell(Num_pop,1); % ISI histogram
for pop_ind = 1:Num_pop
    if nnz(num_spikes{pop_ind}) > 0
        for i = 1:N(pop_ind)
            if spike_count{pop_ind}(i) >= 2
                spike_temp = find(spike_hist{pop_ind}(i,:)); % in time step in lieu of ms!!!
                Dt_temp = (spike_temp(2:end)-spike_temp(1:end-1))*dt; % in ms
                ISI_dist{pop_ind} = [ISI_dist{pop_ind} Dt_temp];
            end
        end
    end
end
% record results
R.Analysis.ISI_dist = ISI_dist;

end

function [R] = rate_sort(R, re_arrange)
if re_arrange == 1
    fprintf('\t Re-arranging neuron index according to individual firing rates...\n');
    for pop_ind = 1:R.Num_pop
        if nnz(R.num_spikes{pop_ind}) > 0
            % re-arrange spike_count accrodingly
            [R.Analysis.spike_count{pop_ind},rate_sort_ind] = sort(R.Analysis.spike_count{pop_ind});  
            % re-arrange spike_hist accrodingly
            R.spike_hist{pop_ind} = R.spike_hist{pop_ind}(rate_sort_ind,:);
            % reduced_spike_hist
            R.reduced.spike_hist{pop_ind} = R.reduced.spike_hist{pop_ind}(rate_sort_ind,:);
        end
    end
    R.Analysis.rate_sorted = 1;
    R.Analysis.rate_sort_ind = rate_sort_ind;
else
    R.Analysis.rate_sorted = 0;
end
end

function [R] = get_CV_ISI(R)
fprintf('\t Getting ISI coefficient of variation...\n');
minimum_spikes = 5; % 3? 5?
% dump fields
dt = R.dt;
N = R.N;
Num_pop = R.Num_pop;
num_spikes = R.num_spikes;
spike_hist = R.spike_hist;
spike_count = R.Analysis.spike_count;
%
CV_ISI = cell(Num_pop,1);% CV_ISI = STD[ISI]/MEAN[ISI], coefficient of variation (note that STD = sqrt(VAR) )
for pop_ind = 1:Num_pop
    if nnz(num_spikes{pop_ind}) > 0
        CV_ISI{pop_ind} = zeros(N(pop_ind),1);
        % warning for insufficient data
        if min(spike_count{pop_ind}) < minimum_spikes;
            warning('insufficient data for estimating CV_ISI for each neuron!');
        end
        % get CV_ISI for neurons with enough spikes
        for i = 1:N(pop_ind)
            spike_temp = (find(spike_hist{pop_ind}(i,:)))*dt; % in (ms) !
            if length(spike_temp) >= minimum_spikes
                ISI_temp = spike_temp(2:end)-spike_temp(1:end-1);
                CV_ISI{pop_ind}(i) = std(ISI_temp)/mean(ISI_temp);% CV
            end
        end
        
    end
end
% CV_ISI averaged over the network
CV_ISI_overall = 0;
for pop_ind = 1:Num_pop
    CV_ISI_overall = CV_ISI_overall+sum(CV_ISI{pop_ind})/sum(N);
end

% Record results
R.Analysis.CV_ISI = CV_ISI;
R.Analysis.CV_ISI_overall = CV_ISI_overall;

end

function [R] = get_CC(R)
fprintf('\t Getting correlation coefficient...\n');
corrcoef_sample_num = 10^4; % number of sampling pairs
corrcoef_matrix_elements = 10^2; % say, 100 neurons pair-wise gives 100 x 100 matrix, for each population
CC_kernel_width = 40; % ms, kernel length

% Dumpe fields
N = R.N;
spike_hist = R.reduced.spike_hist;
dt = R.reduced.dt;
Num_pop = R.Num_pop;
num_spikes = R.num_spikes;

% Define kernel
kernel_type = 'square';
CC_kernel = spike_train_kernel_YG(CC_kernel_width, dt, kernel_type);

%%%%%%%%%% sample pairs from entire network (all the populations)
NNZ_product = 1;
for pop_ind = 1:Num_pop
    NNZ_product = NNZ_product * nnz(num_spikes{pop_ind}); % if any one of them is zero
end

if NNZ_product > 0
    [popA, indA, popB, indB] = pairs_sample_from_network(N,corrcoef_sample_num);
    corrcoef_sample_num = length(popA); % Note that length(popA) <= corrcoef_sample_num !!!!!!!!!!!!!
    % Calculate pair-wise corrcoef
    CC_sample = zeros(1,corrcoef_sample_num);
    for i = 1:corrcoef_sample_num
        CC_sample(i) = CorrCoefYG(spike_hist{popA(i)}(indA(i),:), spike_hist{popB(i)}(indB(i),:), CC_kernel);
%         %%%%%%% display progress
%         if i > 1
%             fprintf(repmat('\b',1,10));
%         end
%         if i < corrcoef_sample_num
%             fprintf('%10g',i);
%         end
    end
    % Record results
    R.Analysis.CC_kernel_type = kernel_type;
    R.Analysis.CC_kernel_width = CC_kernel_width;
    R.Analysis.CC_network_sample = CC_sample;
    R.Analysis.CC_network_sample_mean = mean(CC_sample);
end

%%%%%%%%%% corrcoef matrix within each population
CC_pop_sample_matrix = cell(Num_pop,1);
for pop_ind = 1:Num_pop
    if nnz(num_spikes{pop_ind}) > 0
        if N(pop_ind) > corrcoef_matrix_elements
            sample = sort(randperm(N(pop_ind), corrcoef_matrix_elements));
        else
            sample = 1:N(pop_ind);
        end
        sample_size = length(sample);
        CC_pop_sample_matrix{pop_ind} = zeros(sample_size);
        for A = 1:sample_size-1
            for B = (A+1):sample_size
                CC_pop_sample_matrix{pop_ind}(A,B) = CorrCoefYG(spike_hist{pop_ind}(sample(A),:), spike_hist{pop_ind}(sample(B),:), CC_kernel);
                CC_pop_sample_matrix{pop_ind}(B,A) = CC_pop_sample_matrix{pop_ind}(A,B); % symmetry
                % display progress
                AB = (A-1)*sample_size+B;
                if AB > 1
                    fprintf(repmat('\b',1,10));
                end
                if AB < sample_size*sample_size
                    fprintf('%10g',AB);
                end
            end
        end
    end
end

% Record results
R.Analysis.CC_pop_sample_matrix = CC_pop_sample_matrix;

end

function [R] = generate_comments(R)

fprintf('\t Gennerating auto comments...\n');
% comment line length
cmt_llength = 100;
% dump field
ExplVar = R.ExplVar;
step_killed = R.step_killed;
Hz_overall = R.Analysis.Hz_overall;
% Comments
comments = ' ';
if step_killed >= 0
comments = sprintf('Runaway killed at step %g, ', step_killed);
end
comments = [comments sprintf('Mean firing rate %.4g Hz, ', Hz_overall)];
% Expl variables
fname_cell = fieldnames(ExplVar);
for f = 1:length(fname_cell)
    fn = fname_cell{f};
    if isnumeric(ExplVar.(fn))
        comments = [comments sprintf('%s = %.4g, ', fn, ExplVar.(fn)) ];
    else
        comments = [comments sprintf('%s = %s, ', fn, ExplVar.(fn)) ];
    end
end
% break comments into multiple lines if too long
line_num = ceil(length(comments) / cmt_llength);
if line_num > 1
    comments = [comments repmat(' ', 1, line_num*cmt_llength-length(comments)) ]; % pad with white space
    comment_lines = cell(0,1);
    for ll = 1:line_num
        comment_lines{ll} = comments((ll-1)*cmt_llength+1:ll*cmt_llength);
    end
    comments = comment_lines;
end
% record results
R.comments = comments;

end

%%% helper function
function [popA, indA, popB, indB] = pairs_sample_from_network(N, corrcoef_sample_num)
% e.g.
% N = [4000, 1000]; sample_num = 10^4
% popA = [1 2 1 1 1 2    1    2 ....] (10^4 elements)
% indA = [1 1 2 3 4 1000 4000 8 ...] (10^4 elements)
%
% Note that length(popA) <= corrcoef_sample_num !!!!!!!!!!!!!

N_tot = sum(N);
ind_tot2pop = []; % total index to pop id 1:Num_pop
ind_tot2ind = []; % total index to pop index 1:N(i)
for i = 1:length(N)
    ind_tot2pop = [ind_tot2pop i*ones(1,N(i))];
    ind_tot2ind = [ind_tot2ind 1:N(i)];
end
if N_tot^2 >= corrcoef_sample_num
    pairRaw = randperm(N_tot^2, corrcoef_sample_num); % pairRaw[.] = ind_A*ind_B
else
    pairRaw = 1:N_tot^2;
end
pairs = [floor((pairRaw-1)/N_tot)+1; ...
    mod(pairRaw-1, N_tot)+1];
pairs(:, (pairs(1,:) == pairs(2,:))) = []; % no self-comparision

popA = ind_tot2pop(pairs(1,:)); % population index of neuron A
indA = ind_tot2ind(pairs(1,:)); % neuron index of neuron A
popB = ind_tot2pop(pairs(2,:)); % population index of neuron B 
indB = ind_tot2ind(pairs(2,:)); % neuron index of neuron B
            
end