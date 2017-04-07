function [dist_M, dist_W,dist_V,dist_H] = fhindeptest2(x,y)

n = length(x);

u = pobs(x);
v = pobs(y);
[u,I] = sort(u);
v = v(I);
UV = [u v];

dist_M = point_to_line(UV, repmat([0,0],n,1), repmat([1,1],n,1));
dist_W = point_to_line(UV, repmat([0,1],n,1), repmat([1,0],n,1));   % could also be done w/ a fliplr
dist_V = point_to_line(UV, repmat([0.5,0],n,1), repmat([0.5,1],n,1));
% dist_H = point_to_line(UV, repmat([0,0.5],n,1), repmat([1,0.5],n,1));
dist_H = point_to_line(fliplr(UV), repmat([0.5,0],n,1), repmat([0.5,1],n,1));

end

