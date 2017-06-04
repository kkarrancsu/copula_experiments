% domain experiments for the cos
M = 500;
x = rand(M,1)*10-15;
y = rand(M,1)*10-15;
z = y.*sin(x) - x.*cos(y);
% scatter3(x,y,z)

regions = cos_regions_3d(x,y,z)
