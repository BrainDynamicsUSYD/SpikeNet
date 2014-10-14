function [ Result_cell ] = AnalyseYG( Result_cell, varargin )
%
%   Detailed explanation goes here
disp('AnalyseYG...');
tic;


% choices
re_arrange = 0;
get_corrcoef = 0;

for i = 1:(length(varargin)/2)
    eval([varargin{i*2-1}, '=', num2str(varargin{i*2}) ]);
end


Result_num = length(Result_cell);
for r_num = 1:Result_num
    
    
    % Get Firing rate
    Result_cell{r_num} = get_rate(Result_cell{r_num});
    
%     % Get ISI distribution
%     Result_cell{r_num} = get_ISI(Result_cell{r_num});
    
    
%     % Re-arrange the neuron index
%     % >> according to firing rate
%     Result_cell{r_num} = rate_sort(Result_cell{r_num}, re_arrange);
    
    
%     % Irregularity measures for the network based on CV_ISI
%     Result_cell{r_num} = get_CV_ISI(Result_cell{r_num});
%     
    % Correlation coefficient (zero time-lag)
    if get_corrcoef == 1
        Result_cell{r_num} = get_CC_pop(Result_cell{r_num});
    end
    
    
    % Generate comments
    Result_cell{r_num} = generate_comments(Result_cell{r_num});
    
end
toc;
end



function [R] = get_ISI(R)
% dump fields
dt = R.dt;
N = R.N;
Num_pop = R.Num_pop;
num_spikes = R.num_spikes;
spike_hist = R.spike_hist;
rate = R.Analysis.rate;
%
fprintf('\t Getting ISI distribution...\n');
ISI_dist = cell(Num_pop,1); % ISI histogram
for pop_ind = 1:Num_pop
    if nnz(num_spikes{pop_ind}) > 0
        for i = 1:N(pop_ind)
            if rate{pop_ind}(i) >= 2
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
            % re-arrange rate accrodingly
            [R.Analysis.rate{pop_ind},rate_sort_ind] = sort(R.Analysis.rate{pop_ind});  
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
rate = R.Analysis.rate;
%
CV_ISI = cell(Num_pop,1);% CV_ISI = STD[ISI]/MEAN[ISI], coefficient of variation (note that STD = sqrt(VAR) )
for pop_ind = 1:Num_pop
    if nnz(num_spikes{pop_ind}) > 0
        CV_ISI{pop_ind} = zeros(N(pop_ind),1);
        % warning for insufficient data
        if min(rate{pop_ind}) < minimum_spikes;
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

function [R] = get_pop_CC_matrix(R)
fprintf('\t Getting correlation coefficient matrix for each population...\n');

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

% corrcoef matrix within each population
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



