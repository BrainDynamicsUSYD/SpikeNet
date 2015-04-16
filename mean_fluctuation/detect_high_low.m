

function [low_du_tot, high_du_tot] =  detect_high_low( y_in, thre, high_du_min_steps, dt, cutoff  )


% if nargin == 1
%     thre = 0;
% end
% 
% if nargin <= 2
%     high_du_min_steps = 0;
% end

low_du_tot = [];
high_du_tot = [];
[Mnum, T] = size(y_in);
if Mnum > T
    y_in = y_in';
end
[Mnum, ~] = size(y_in);

for M = 1:Mnum
    seq = y_in >= thre;
    seq = persistence_requirement( seq, high_du_min_steps );
    [~, high_du, low_du, ~, ~] = seq_postprocess(seq, dt, cutoff);
    low_du_tot = [low_du_tot low_du];
    high_du_tot = [high_du_tot high_du];
end




% % % % previous solution

% if nargin == 1
%     thre = 0;
% end
% 
% low_du_tot = [];
% high_du_tot = [];
% [Mnum, T] = size(y_in);
% if Mnum > T
%     y_in = y_in';
% end
% [Mnum, T] = size(y_in);
% 
% 
% for M = 1:Mnum
%     y = y_in(M,:);
%     
%     high_du = [];
%     low_du = [];
%     state_pre = (y(1) >= thre) ; % 1 for high, 0 for low
%     if state_pre == 1
%         high_du = 1;
%     else
%         low_du = 1;
%     end
% 
%     for i = 2:T
%         
%         % detect
%         state_now = (y(i) >= thre);
%         
%         if state_now == 1 && state_now == state_pre
%             high_du(end) = high_du(end) + 1;
%         elseif state_now == 0 && state_now == state_pre
%             low_du(end) = low_du(end) + 1;
%         elseif state_now == 1 && state_now ~= state_pre
%             high_du = [high_du 1];
%         elseif state_now == 0 && state_now ~= state_pre
%             low_du = [low_du 1];
%         end
%         state_pre = state_now;
%         
%     end
%     
%     if length(high_du) > 2
%         high_du(1) = [];
%         high_du(end) = [];
%     else
%         high_du = [];
%     end
%     if length(low_du) > 2
%         low_du(1) = [];
%         low_du(end) = [];
%     else
%         low_du = [];
%     end
%     
% 
%     low_du_tot = [low_du_tot low_du];
%     high_du_tot = [high_du_tot high_du];
% end

end