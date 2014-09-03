function visualise_potential_on_layout( p, V, V_range, cmap )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

max_V = max(V_range);
min_V = min(V_range);

% Continuous change of FaceColor
% map V to RGB values using colormap
V(V>max_V) = max_V;
V(V<min_V) = min_V;
C = (V-min_V)/(max_V-min_V);
Cind = ceil(C*length(cmap(:,1)));
Cind(Cind == 0) = 1;
Color = cmap(Cind,:); % continuous color mapping from V

% Continuous change of FaceAlpha
max_Alpha = 0.4;
min_Alpha = 0.04;
Alpha = (V-min_V)/(max_V-min_V)*(max_Alpha-min_Alpha)+min_Alpha;  % continuous alpha mapping from V

for i = 1:length(V)
    % % Discrete change of FaceAlpha
    %     if max_V-V(i) < 0.1*(max_V-min_V)
    %         set(p{i}, 'FaceColor', Color(i,:),'FaceAlpha',0.5);
    %     else
    %         set(p{i}, 'FaceColor', Color(i,:),'FaceAlpha',0.1);
    %     end
    
    % Continuous change of FaceAlpha
    set(p{i}, 'FaceColor', Color(i,:),'FaceAlpha',Alpha(i));
end

end

