clc;
clear;
close all;



    
% generate Matrix rain code

% Say we have 
%   (1) the screen grid size = h-by-w
%   (2) the number of time steps = t
%   (3) the number of random charactors = C

h = 40;
w = 20;
t = 50;
c_lib = 'A':'Z';

c_length = length(c_lib);

vidObj = VideoWriter('test_raining_code.avi');
open(vidObj);
                    
matrix_green = [1/3 1 1]; % HSV

head_len = 3;
tail_len_range = [5 10];

life_range = [15 30];

rain_pos = [];
rain_pos_rev = [];
h_pos = [];
w_pos = [];
rain_tail_len = [];

v_pos = []; % different speed of the tick
v_range = [0.3 0.6]; % less than 1
v_tick = []; % when the tick runs over an integer, the rain falls


c_map = randi([1, c_length], 1, h*w); % charactor map, a vector (use vec2ind)
c_change_rate = 0.05; % charactor map randomly update rate

hsv_v_map_rand = unifrnd(0.6, 1, 1, w*h); % hsv value map
    

new_rate = 0.1; % rate of generating new rain

figure(1);
set(gcf,'color','k')%, 'visible','off')
% whitebg(gcf,'k')
axis off;
axis([0 w+1 0 h+1]); %hold on;
set(gca,'YDir','reverse')
set(gca,'color','k')
hold on;

if life_range(1) < head_len + max(tail_len_range)
    error('life_range(1) < head_len + tail_len!')
end

tic
for i = 1:t
    
    % generate new rains
    new_rain = setdiff(1:w, unique(w_pos(h_pos < h/2)));
    new_rain(rand(1,length(new_rain)) > new_rate) = []; 
    w_pos = [w_pos new_rain]; %#ok<*AGROW>
    h_pos = [h_pos ones(1,length(new_rain))]; 
    v_pos = [v_pos unifrnd(v_range(1), v_range(2), 1,length(new_rain))];
    v_tick = [v_tick zeros(1, length(new_rain))];
    
    new_rain_length = randi(life_range, 1,length(new_rain));
    rain_pos = [rain_pos new_rain_length]; % counting from N to 1 (from head to tail)
    rain_pos_rev = [rain_pos_rev ones(1,length(new_rain))];% counting from 1 to N (from head to tail)
    rain_tail_len = [rain_tail_len randi(tail_len_range, 1,length(new_rain))]; % random tail length

    % falling down
    f_fall =  h_pos == 1 & (rain_pos > 1) & ((floor(v_tick+v_pos)-floor(v_tick)) == 1);
    
    w_pos_new = {w_pos, w_pos(f_fall)};
    rain_pos_new = {rain_pos, rain_pos(f_fall)-1};
    h_pos_new = h_pos;
    h_pos_new((floor(v_tick+v_pos)-floor(v_tick)) == 1) = h_pos_new((floor(v_tick+v_pos)-floor(v_tick)) == 1) + 1;
    h_pos_new = {h_pos_new, ones(1, sum(f_fall))};
    v_pos_new = {v_pos, v_pos(f_fall)};
    v_tick_new = {v_tick, v_tick(f_fall)};
    
    rain_pos_rev_new = {rain_pos_rev, rain_pos_rev(f_fall)+1};
    rain_tain_len_new = {rain_tail_len, rain_tail_len(f_fall)};
    
    w_pos = cell2mat(w_pos_new);
    rain_pos = cell2mat(rain_pos_new);
    h_pos = cell2mat(h_pos_new); 
    v_pos = cell2mat(v_pos_new);
    v_tick = cell2mat(v_tick_new)+ v_pos;
    rain_pos_rev = cell2mat(rain_pos_rev_new);
    rain_tail_len = cell2mat(rain_tain_len_new);

    % delete those falling out of screen
    w_pos(h_pos > h) = [];
    rain_pos(h_pos > h) = [];
    v_pos(h_pos > h) = [];
    v_tick(h_pos > h) = [];
    rain_pos_rev(h_pos > h) = [];
    rain_tail_len(h_pos > h) = [];
    h_pos(h_pos > h) = [];

    
    
    % update character-to-position map
    c_change = rand(1,h*w) < c_change_rate; 
    c_map( c_change ) = randi([1, c_length], 1, sum(c_change)); 
    c_pos = c_lib(c_map(sub2ind([h w], h_pos, w_pos)));
    
    % generate brightness
    hsv_v_pos = hsv_v_map_rand(sub2ind([h w], h_pos, w_pos));
    % head effect
    hsv_s_pos_head = ones(size(rain_pos_rev));
    is_rain_head = rain_pos_rev <= head_len;
    hsv_s_pos_head(is_rain_head) = 1- 0.8*exp(- (rain_pos_rev(is_rain_head)-1)/head_len);
    % tail effect
    hsv_v_pos_tail = ones(size(rain_pos_rev));
    is_rain_tail = rain_pos <= rain_tail_len;
    hsv_v_pos_tail(is_rain_tail) = 1-exp(-(rain_pos(is_rain_tail)-1)./rain_tail_len(is_rain_tail)*2);
    
 
    
    % draw
    cla;
    for j = 1:length(w_pos)
        color_tmp = matrix_green;
        color_tmp(3) = color_tmp(3)*hsv_v_pos(j)*hsv_v_pos_tail(j);
        color_tmp(2) = hsv_s_pos_head(j);
        color_tmp = hsv2rgb(color_tmp);
        if rain_pos_rev(j) == 1
            text( w_pos(j), h_pos(j), c_pos(j), 'color', color_tmp , ...
                'FontWeight', 'bold', 'HorizontalAlignment','center','VerticalAlignment','top',...
                'FontName','Courier')
        else
            text( w_pos(j), h_pos(j), c_pos(j), 'color', color_tmp , ...
                'FontWeight', 'bold', 'HorizontalAlignment','center','VerticalAlignment','top',...
                'FontName','Courier')
        end
        
    end
    drawnow;
    writeVideo(vidObj, getframe(gcf));

    i/t
end
toc;
close(vidObj);






