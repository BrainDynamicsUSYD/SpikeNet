function visualise_potential_on_circle_packing_layout( p, V, V_range, cmap )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here



mode = 'state';

switch mode
    case 'continuous'
        % Continuous change of FaceColor
        % map V to RGB values using colormap
        max_V = max(V_range);
        min_V = min(V_range);
        V(V>max_V) = max_V;
        V(V<min_V) = min_V;
        C = (V-min_V)/(max_V-min_V); % [0,1]
        Cind = ceil(C*length(cmap(:,1)));
        Cind(Cind == 0) = 1;
        Color = cmap(Cind,:); % continuous color mapping from V
    case 'state'
        % State parameter: check that they are consistent with the model!
        V_th = -55;
        V_rt = -75;
        Color = ones(length(V),3); % white
        Color(V>=V_th,:) = repmat([1 0 0],sum(V>=V_th),1); % red
        Color(V==V_rt,:) = repmat([0 0 0],sum(V==V_rt),1); % black
end




% Update scatter color data
set(p,'CData',Color); % Toooooooooooooooooooooo slow!!!!


end

