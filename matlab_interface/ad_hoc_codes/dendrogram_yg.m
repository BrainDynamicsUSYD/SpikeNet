function dendrogram_yg(level, n1, n2)
% if given node n1 and n2, the path between them will be shown by bold
% line

height = 1;
linewidth = 3;


set(gcf,'color','w');
hold on;

n = 2^(level-1);
for i = 1:level
    X = (1+2^(i-1))/2:2^(i-1):n;
    Xc{i} = X;
    Y = ones(1,n*2^(-(i-1)))*i*height;
    Yc{i} = Y;
    draw_level(X,Y, height);
    
    text( 0, i*height, num2str(i),'VerticalAlignment','middle' );
end

for i = 1:n
    text( i, 0.6, num2str(i), 'HorizontalAlignment','center');
end

axis off;
xlim([0 n+1])
ylim([0 (level+1)*height])


% show path
if nargin == 3
    
    path1 = [];
    path2 = [];
    for i = 1:level
        C1 = ceil( ( n1/2^(i-1)) );
        C2 = ceil( ( n2/2^(i-1)) );
        if C1 ~= C2
            path1 = [path1; Xc{i}(C1) Yc{i}(C1); Xc{i}(C1) Yc{i}(C1)+height ];
            path2 = [path2; Xc{i}(C2) Yc{i}(C2); Xc{i}(C2) Yc{i}(C2)+height ];
        else
            path = [path1; flipud(path2)];
            for p = 1:length(path(:,1))-2
                line(path(p:p+1, 1), path(p:p+1, 2), 'LineWidth',linewidth)
            end
            vectarrow(path(end-1, :), path(end, :), linewidth)
            break;
        end
    end
end


end



function draw_level(X,Y,height)
% X = [node1 node2 ... nodeN]
% (node1,node2) will be a branch

if nargin == 2
    height = 1;
end

for i = 1:(length(X)/2)
    draw_branch([X(i*2-1) Y(i*2-1)], [X(i*2) Y(i*2)], height)
end

if length(X) == 1
    tip = 0.2;
    draw_branch([X Y], [X Y], height*tip)
end
   
end

function draw_branch(leaf1, leaf2, height)
% leaf = [x y]
if nargin == 2
    height = 1;
end

X1 = leaf1(1);
Y1 = leaf1(2);
X2 = leaf2(1);
Y2 = leaf2(2);

H = height;

line([X1;X1], [Y1;Y1+H]);
line([X2;X2], [Y2;Y2+H]);
line([X1;X2], [Y1;Y2]+H);
    
end


function vectarrow(p0,p1, linewidth)
%   vectarrow(p0,p1) plots a line vector with arrow pointing from point p0
%   to point p1.
%
%       2D vector
%       p0 = [1 2];     % Coordinate of the first point p0
%       p1 = [4 5];     % Coordinate of the second point p1
%       vectarrow(p0,p1)


alpha = 0.2;  % Size of arrow head relative to the length of the vector
beta = 0.3;  % Width of the base of the arrow head relative to the length

x0 = p0(1);
y0 = p0(2);
x1 = p1(1);
y1 = p1(2);
x2 = x0 + (1-alpha+eps)*(x1-x0);
y2 = y0 + (1-alpha+eps)*(y1-y0);
plot([x0;x2],[y0;y2], 'LineWidth',linewidth);   % Draw a line between p0 and p1

p = p1-p0;


hu = [x1-alpha*(p(1)+beta*(p(2)+eps)); x1; x1-alpha*(p(1)-beta*(p(2)+eps))];
hv = [y1-alpha*(p(2)-beta*(p(1)+eps)); y1; y1-alpha*(p(2)+beta*(p(1)+eps))];

hold on
fill(hu(:),hv(:), 'b');

end