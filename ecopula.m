function [e_cop_points] = ecopula(u,v)
% computes the empirical copula, only for the points provided, rather than
% the entire grid

n = length(u); e_cop_points = zeros(1,n);
for ii=1:n
    e_cop_points(ii) = sum(u(ii)>=u & v(ii)>=v)/(n+1);
end

end