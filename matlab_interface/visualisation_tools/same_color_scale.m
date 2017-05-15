function same_color_scale(ah)
% this function takes in a vector of axis handles and unify their color
% scales
%
% Usage:
%   for i = 1:6
%       ah(i) = subplot(2,3,i);
%   end
%   same_color_scale(ah);

a = cell2mat(get(ah,'CLIM'));
cl = [min(a(:,1)) max(a(:,2))];
set(ah,'CLIM',cl);

end