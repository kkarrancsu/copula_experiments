function [u,v] = pobs_sorted(x,y)

u = pobs(x);
v = pobs(y);

[u,I] = sort(u);
v = v(I);

end