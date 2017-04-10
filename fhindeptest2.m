function [dist_M, dist_W,dist_V,dist_H] = fhindeptest2(x,y)

n = length(x);

u = pobs(x);
v = pobs(y);
[u,I] = sort(u);
v = v(I);
UV = [u v];
UV_shuffled = shuffle_pobs(UV);

dist_M = point_to_line(UV, repmat([0,0],n,1), repmat([1,1],n,1));
dist_W = point_to_line(UV, repmat([0,1],n,1), repmat([1,0],n,1));   % could also be done w/ a fliplr
dist_V = point_to_line(UV_shuffled, repmat([0,0],n,1), repmat([1,1],n,1));
dist_H = point_to_line(UV_shuffled, repmat([0,1],n,1), repmat([1,0],n,1));

end

