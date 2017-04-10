function [U] = shuffle_pobs(u,seed)

if(nargin>1)
    rng(seed);
end

n = size(u,1);
d = size(u,2);
U = zeros(size(u));
for ii=1:d
    x = u(:,ii);
    I = randperm(n);
    U(:,ii) = x(I);
end

end