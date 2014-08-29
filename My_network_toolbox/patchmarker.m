function [ p_output ] = patchmarker(   x,y,z,r,varargin)
% Plot marker as patch in 3D using small sphere to control transparency
%
% SYNTAX:
%     patchmarker(x,y,z,patchmarkersize)
%     patchmarker(x,y,z,patchmarkersize,'PropertyName',propertyvalue,...)
%     p = patchmarker(...)
%
%     if x, y, z are vector, p will be column vector (cell array) of
%     handles
%
% PROPERTIES: 
%     Accepts all parameter-values accepted by PATCH.
%     But: 'EdgeColor' is 'none', 'MarkerEdgeColor' is 'none' and 
%          'MarkerFaceColor' is 'none'. They shall not be changed.

N = 3; % keep it small but make small sphere used as 3D marker still look good, 3 seems to be the minimum
[xs, ys, zs] = sphere(N); % creat unit sphere

num = length(x);
if num > 1
    p_output = cell(num,1);
end

for i = 1:num
    % creat patch
    fvs = surf2patch(x(i)+xs*r(i),y(i)+ys*r(i),z(i)+zs*r(i)); % face and vertices struct
    p_temp = patch(fvs, 'FaceColor', 'r', 'FaceAlpha', 0.03,... % default FaceAlpha is 0.03
        'EdgeColor', 'none', 'MarkerEdgeColor','none', 'MarkerFaceColor', 'none', 'erasemode','normal');
    % read non-default settings
    for j = 1:length(varargin)/2
        set(p_temp,varargin{2*j-1},varargin{2*j})
    end
    % store handles for output
    if nargout == 0
        clear p_temp
    elseif num  == 1
        p_output = p_temp;
    else
        p_output{i} = p_temp;
    end
end

end

