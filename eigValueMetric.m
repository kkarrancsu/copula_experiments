function [metric] = eigValueMetric(U)
R = nancov(U);
zz = eig(R);
metric = zz(2)/zz(1);
end