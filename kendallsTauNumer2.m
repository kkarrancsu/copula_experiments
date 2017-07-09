function [K] = kendallsTauNumer2(X,Y)

uu = repmat(X,1,length(X));
vv = repmat(Y,1,length(Y));

n = size(uu,1);
opMat = -1*ones(size(uu));
opMat(1:n+1:n*n) = 1;

U = uu.*opMat;
V = vv.*opMat;
cmpMatrixU = repmat(diag(U)',size(uu,1),1);
cmpMatrixV = repmat(diag(V)',size(vv,1),1);
sumMatU = U+cmpMatrixU;
sgnMatU = tril(sign(sumMatU),-1);
sumMatV = V+cmpMatrixV;
sgnMatV = tril(sign(sumMatV),-1);

K = sum(sum(sgnMatU.*sgnMatV));

end