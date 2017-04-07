function [dist_M, dist_W,dist_V,dist_H] = fhindeptest(x,y,distType)

u = pobs(x);
v = pobs(y);
[u,I] = sort(u);
v = v(I);
UV = [u v];

% compute distance between M(u,v) - Frechet-Hoeffding upper bound
v_M = u;
UV_M = [u v_M];
if(strcmpi(distType,'euclidean'))
    dist_M = computeEuclideanDist(UV,UV_M);
elseif(strcmpi(distType,'signed'))
    dist_M = computeSignedDist(UV,UV_M);
end

% compute distance between W(u,v) - Frechet-Hoeffding lower-bound
v_W = 1-u;
UV_W = [u v_W];
if(strcmpi(distType,'euclidean'))
    dist_W = computeEuclideanDist(UV,UV_W);
elseif(strcmpi(distType,'signed'))
    dist_W = computeSignedDist(UV,UV_W);
end

uu = 0.5*ones(length(x),1);
UV_vertical = [u v];
if(strcmpi(distType,'euclidean'))
    dist_V = computeEuclideanDist(UV,UV_vertical);
elseif(strcmpi(distType,'signed'))
    dist_V = computeSignedDist(UV,UV_vertical);
end

vv = uu;
UV_horizontal = [u vv];
if(strcmpi(distType,'euclidean'))
    dist_H = computeEuclideanDist(UV,UV_horizontal);
elseif(strcmpi(distType,'signed'))
    dist_H = computeSignedDist(UV,UV_horizontal);
end


end

function [dist] = computeEuclideanDist(X,Y)
dist = sqrt(sum((X-Y).^2,2));
end

function [dist] = computeSignedDist(X,Y)
dist = sqrt(sum((X-Y).^2,2)).*sign(X(:,2)-Y(:,2));
end