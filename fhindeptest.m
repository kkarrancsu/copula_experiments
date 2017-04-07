function [dist_M, dist_W] = fhindeptest(x,y)

u = pobs(x);
v = pobs(y);
[u,I] = sort(u);
v = v(I);
UV = [u v];

% compute distance between M(u,v) - Frechet-Hoeffding upper bound
v_M = u;
UV_M = [u v_M];

dist_M = sqrt(sum((UV-UV_M).^2,2));

% compute distance between W(u,v) - Frechet-Hoeffding lower-bound
v_W = 1-u;
UV_W = [u v_W];
dist_W = sqrt(sum((UV-UV_W).^2,2));

end