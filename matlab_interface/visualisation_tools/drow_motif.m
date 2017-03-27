function drow_motif( links, centre, r)

r = 1;
centre = [0 0];

node_size = 20;
arrow_size = 10;


x = centre(1) + r*cosd([-30 90 210]);
y = centre(2) + r*sind([-30 90 210]);
plot(x,y,'.','MarkerSize', node_size);
hold on;

for i = 1:length(links(:,1))
    arrow([x(links(i,1)) y(links(i,1))],[x(links(i,2)) y(links(i,2))],arrow_size)
end


end