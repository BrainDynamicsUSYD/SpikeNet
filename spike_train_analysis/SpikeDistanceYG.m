function [ Result_cell ]  = SpikeDistanceYG( Result_cell )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   This function is not good!!

disp('SpikeDistanceYG...');
tic;

minimum_spikes = 3; % minimum data size
sample_size = 50; % sampling the population


Result_num = length(Result_cell);

for r_num = 1:Result_num
    % Dump fields
    Num_pop = Result_cell{r_num}.Num_pop;
    dt = Result_cell{r_num}.dt;
    step_tot = Result_cell{r_num}.step_tot;
    N = Result_cell{r_num}.N;
    spike_hist_tot = Result_cell{r_num}.spike_hist_tot;
    spike_count = Result_cell{r_num}.spike_count;
    T = (1:step_tot)*dt;

    % Synchrony measurment
    SPIKE_dissimilarity = cell(Num_pop,1);
    for pop_ind = 1:Num_pop
        if min(spike_count{pop_ind}) >= minimum_spikes;
            % random sampling
            sample_ind = sort(randperm(N(pop_ind),sample_size));
            spike_mat = zeros(sample_size,1);
            for i = 1:sample_size
                spike_temp = (find(spike_hist_tot{pop_ind}(sample_ind(i),:)))*dt; % in (ms) !
                spike_mat(i,1:length(spike_temp)) = spike_temp;
            end
            
            
%             % Separate-and-conquer to avoid memory problem: toooooooo
%             % slow!!
%             plot = 0;
%             dissimilarity_profile = zeros(size(T));
%             distance_matrix = zeros(sample_size,sample_size);
%             for i = 1:sample_size-1
%                 for j = (i+1):sample_size
%                     result_temp = Kreuz_SPIKE_dissimilarity( spike_mat([i j],:), T, plot );
%                     dissimilarity_profile = dissimilarity_profile+result_temp.dissimilarity_profiles{1};
%                     distance_matrix(i,j) = result_temp.overall_dissimilarities;
%                     distance_matrix(j,i) = distance_matrix(i,j);
%                 end
%             end
%             dissimilarity_profile = dissimilarity_profile/(sample_size*(sample_size-1)/2);
%             overall_dissimilarity = sum(sum(triu(distance_matrix)))/(sample_size*(sample_size-1)/2);
%             SPIKE_dissimilarity{pop_ind}.overall_dissimilarity = overall_dissimilarity;
%             SPIKE_dissimilarity{pop_ind}.distance_matrix = distance_matrix;
%             SPIKE_dissimilarity{pop_ind}.dissimilarity_profile = dissimilarity_profile;

            % All-together: memory problem!!!!!!
            plot = 0;
            SPIKE_dissimilarity{pop_ind} = Kreuz_SPIKE_dissimilarity( spike_mat, T, plot );
            

        else
            disp('Insufficient data for estimating SPIKE dissimilarity!');
        end
    end
    
    % Write data
    Result_cell{r_num}.SPIKE_dissimilarity = SPIKE_dissimilarity;
    
end
    

toc;
end

