function [cosValue, RR] = cos2d_v4(x,y)
%COS2D_v4 - computes the 2-D cos metric

general_debug_print = 0;
lambda_debug_print = 0;

n = length(x);
d = 2;

u = pobs(x);
v = pobs(y);
[u,I] = sort(u);
v = v(I);

U = [u v];

E_copula = zeros(n,1);
Pi_copula = u.*v;
for ii=1:n
    E_copula(ii) = sum(u(ii)>=u(1:n) & v(ii)>=v(1:n))/(n+1);
end

idxVec = 1:n;
RR = [idxVec' U E_copula Pi_copula ones(n,1)*-999]; % the last column will be filled in 
                                                    % by the computed domain

E_COPULA_IDX = 4;
PI_COPULA_IDX = 5;
DOMAIN_IDX = 6;
                                                    
domainIdx = 1;
zz = RR(:,E_COPULA_IDX)>RR(:,PI_COPULA_IDX);
zz_diff = diff(zz);
I = find(zz_diff~=0);
startIdx = 1;
for ii=I'
    endIdx = ii;
    if(endIdx==startIdx)
        continue;
    end
    RR(startIdx:endIdx,DOMAIN_IDX) = domainIdx;
    domainIdx = domainIdx + 1;
    startIdx = endIdx+1;
end
if(startIdx<=n)
    RR(startIdx:n,DOMAIN_IDX) = domainIdx;
end

% compute lambda_lo for the first domain
pi1 = RR(1,PI_COPULA_IDX);
e1 = RR(1,E_COPULA_IDX); 
if(is_close(e1,pi1))
    lambda_lo = 0;
else
    lambda_lo = 1;
end

% compute the metric for each domain
numDomains = domainIdx;
domainMetrics = zeros(2,numDomains);
for ii=1:numDomains
    % get the data associated w/ this domain
    II = RR(:,6)==ii;
    domainData = RR(II,:);
    domainLen = length(find(II==1));
    minAbsoluteIdx = min(RR(II,1));
    maxAbsoluteIdx = max(RR(II,1));
    if(general_debug_print)
        fprintf('[cos] -- Domain -> %d:%d\n',minAbsoluteIdx,maxAbsoluteIdx);
    end

    iMax = domainLen;  % because ECopula is sorted
    
%     maxECopVal = RR(maxAbsoluteIdx-1:maxAbsoluteIdx+1,E_COPULA_IDX);
    maxECopVal = RR(maxAbsoluteIdx,E_COPULA_IDX);
    piMaxVal = domainData(iMax,PI_COPULA_IDX);
    WmaxVal = max(1-d+sum(domainData(iMax,2:3)),0);
    MmaxVal = min(domainData(iMax,2:3));
    
    lambda_hi = computeRelativeDistance(WmaxVal,piMaxVal,MmaxVal,maxECopVal,n,domainLen,lambda_debug_print);
    
    gamma = (lambda_lo+lambda_hi)/2;
    
    if(general_debug_print)
        fprintf('[cos] -- lambda_min=%0.02f lambda_max=%0.02f\n',lambda_lo,lambda_hi);
    end

    lambda_lo = lambda_hi;
    
    domainMetrics(1,ii) = gamma;
    domainMetrics(2,ii) = domainLen;
    
end

cosValue = sum(domainMetrics(1,:).*domainMetrics(2,:))/(n+numDomains-1);

end

function [lambda] = computeRelativeDistance(wVal,piVal,mVal,cValVec,n,domainLen,lambda_debug_print)

% cVal_n1 = cValVec(1); cVal_p1 = cValVec(3);
% cVal = cValVec(2);
cVal = cValVec(1);

if(is_close(cVal,piVal))
    if(lambda_debug_print)
        fprintf('[cos] -- >> c\n');
    end
    lambda = 0;
else
    if(cVal>piVal)
        if(lambda_debug_print)
            fprintf('[cos] -- >> a\n');
        end
        lambda = (cVal-piVal)/(mVal-piVal);
    else
        if(lambda_debug_print)
            fprintf('[cos] -- >> b\n');
        end
        lambda = (piVal-cVal)/(piVal-wVal);
    end
end

% test for local-minimum
% if(is_close(cVal,mVal) || is_close(cVal,wVal) || is_close(piVal,mVal) || ...
%    (is_close(piVal,wVal) && domainLen>4) || (cVal==1/(n+1) && domainLen>4) || ...
%    is_close(cVal_n1,cVal_p1) || (is_close(cVal_n1,cVal) && domainLen>4) )
% if(is_close(cVal,mVal) || is_close(cVal,wVal) || is_close(piVal,mVal) || ...
%    (is_close(piVal,wVal) && domainLen>4) || (cVal==1/(n+1) && domainLen>4) )
%     lambda = 1;
% end
if(is_close(cVal,mVal) || is_close(cVal,wVal) || (cVal==1/(n+1) && domainLen>4))
    if(lambda_debug_print)
        fprintf('[cos] -- >> d\n');
    end
    lambda = 1;
end

if(is_close(piVal,mVal) || (is_close(piVal,wVal) && domainLen>4))
    if(lambda_debug_print)
        fprintf('[cos] -- >> f\n');
    end
    lambda = 1;
end

% REDUNDANT TO C above ... shouldn't hit it twice?
% if(is_close(cVal,piVal))
%     fprintf('[cos] -- >> g\n');
%     lambda = 1;
% end

end

function [tf] = is_close(a,b)
tf = 0;
if(abs(a-b)<10*eps)
    tf = 1;
end
end