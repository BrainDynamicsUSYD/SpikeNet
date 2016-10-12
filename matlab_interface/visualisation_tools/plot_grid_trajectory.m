function  plot_grid_trajectory( R )
%UNTITLED Summary of this function goes here
%   periodic boundary condition is removed

dir = R.grid.raw.jump_dir;
dist = R.grid.raw.jump_dist;
x = cumsum(dist.*cos(dir));
y = cumsum(dist.*sin(dir));
plot(x,y);

end

