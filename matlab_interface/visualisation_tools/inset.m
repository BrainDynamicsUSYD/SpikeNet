function [ah1, ah2] = inset(inset_pos)
% when done, use axes(ah1) to switch back to the orginal axes
ah1 = gca;

if nargin == 0
    inset_pos = [0.7 0.7 0.2 0.2];
end
pos = get(gca,'position');
new_pos = [pos(1) + inset_pos(1)*pos(3),...
           pos(2) + inset_pos(2)*pos(4),...
           inset_pos(3)*pos(3),...
           inset_pos(4)*pos(4)];
       
ah2 = axes('Position',new_pos, 'layer','top');
set(ah1, 'layer','bottom', 'color', 'none');

end
