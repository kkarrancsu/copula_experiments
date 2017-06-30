function [z] = ktau_zscore(tau,n)
% computes the z-score for kendall's tau value, given n points

z = tau./( (2*(2*n+5))./(9*n.*(n-1)) );

end