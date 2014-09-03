function [X0,Y0] = CirclePackingLayout(R0)
% Graph Layout using (non-optimal, non-compact) cicle packing.
%
% Circle radii should be (mapped from) total strength/degree of graph
% nodes.
%
% Layout circles are arranged in a counter-clockwise spiral fasion with
% decreasing circle sizes while spinning away from origin.
% 
% Performance: 0.2 sec with 10,000 nodes.

% Sorting
[R,Ind] = sort(R0,'descend'); % decreasing circle sizes

% Initialisation
N = length(R);
X = zeros(1,N);
Y = zeros(1,N);

ShowPlot = 0; % Default value
if ShowPlot == 1
    % Prepare for visualising results
    figure('NumberTitle','Off','Name','Circles');
    axis equal;
    hold on;
    theta = linspace(0,2*pi,100);
    dx = cos(theta);
    dy = sin(theta);
end

% The main while loop
i_current = 1;
while i_current <= N
    % Set/Reset base circle centre and radius
    if i_current == 1 % Set first base circle
        X_base = 0;
        Y_base = 0;
        R_base = R(1);
        % The first base circle is the largest circle
        theta = 0; % polar angle of current circle centres
        X(1) = X_base;
        Y(1) = Y_base;
        i_base = 1;
        i_current = 2;
    
    else % Reset base circle
        % Find new centre
        X_vec = [X_base,X(i_base+1:i_current-1)];
        Y_vec = [Y_base,Y(i_base+1:i_current-1)];
        R_vec = [R_base,R(i_base+1:i_current-1)];
        % Initial guess: weighted by radius
        X_base_new = dot(X_vec,R_vec)/sum(R_vec); % New X_base
        Y_base_new = dot(Y_vec,R_vec)/sum(R_vec); % New Y_base
        L = sqrt((X(i_base+1:i_current-1)-X_base_new).^2+(Y(i_base+1:i_current-1)-Y_base_new).^2)+R(i_base+1:i_current-1);
        R_base_new = max([sqrt((X_base-X_base_new)^2+(Y_base-Y_base_new)^2)+R_base, L]);
        
        
        % Apply random searching for Minimum Enclosing Circle as the new
        % base circle
        trial_num = 200;
        pert = 0.15; % max perturbation strength in percentage of radius
        for jj = 1:trial_num
            X_base_trial = X_base_new + (rand*2-1)*pert*R_base_new;
            Y_base_trial = Y_base_new + (rand*2-1)*pert*R_base_new;
            % Find radius given by the centre
            L = sqrt((X(i_base+1:i_current-1)-X_base_trial).^2+(Y(i_base+1:i_current-1)-Y_base_trial).^2)+R(i_base+1:i_current-1);
            R_base_trial = max([sqrt((X_base-X_base_trial)^2+(Y_base-Y_base_trial)^2)+R_base, L]);
            if R_base_trial < R_base_new
                % disp('Improved');
                X_base_new = X_base_trial;
                Y_base_new = Y_base_trial;
                R_base_new = R_base_trial;
            end
        end

        % set virtual base circle
        X_base = X_base_new;
        Y_base = Y_base_new;
        R_base = R_base_new;
        i_base = i_current-1;
        theta = 0;
    end
    if ShowPlot == 1
        plot(X_base+dx*R_base,Y_base+dy*R_base,'r--');
    end
    
    % Wrap patel circles around current base circle until full
    non_overlapping = true;
    while i_current <= N && non_overlapping == true
        % D_theta (between i_current-1 and i_current)
        if i_current - i_base == 1 % First patel circle
            D_theta = 0; % increment in polar angle
        else% Following patel circles
            aa = R_base*(R_base+R(i_current-1)+R(i_current));
            bb = R(i_current-1)*R(i_current);
            D_theta = acos((aa-bb)/(aa+bb)); % simply applying the cosine theorem
            % Note that due to descending size, D_theta <= pi/3 = 60 degree
        end
        % Check Non-overlapping condition
        % theta_min (between i_base+1 and i_current)
        aa = R_base*(R_base+R(i_base+1)+R(i_current));
        bb = R(i_base+1)*R(i_current);
        theta_min = acos((aa-bb)/(aa+bb));
        if 2*pi-(D_theta+theta) >= theta_min  % Non-overlapping
            theta = D_theta+theta; % increment theta
            X(i_current) = cos(theta)*(R_base+R(i_current))+X_base;
            Y(i_current) = sin(theta)*(R_base+R(i_current))+Y_base;
            if ShowPlot == 1
                plot(X(i_current)+dx*R(i_current),Y(i_current)+dy*R(i_current));
            end
            i_current = i_current + 1;
        else %  overlapping
            non_overlapping = false;
        end
    end
end

% Sorting back
[~,Ind_back] = sort(Ind,'ascend');
X0 = X(Ind_back);
Y0 = Y(Ind_back);
end