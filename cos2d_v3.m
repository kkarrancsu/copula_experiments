function [cosValue, RR] = cos2d_v3(x,y)
%COS2D_v3 - computes the 2-D cos metric
n = length(x);
d = 2;

X = [x y]; 
U = pobs(X);
U = sortrows(U);
u = U(:,1); v = U(:,2);

E_copula = zeros(n,1);
Pi_copula = u.*v;
for ii=1:n
%     E_copula(ii) = sum(u(ii)>=u(1:n) & v(ii)>=v(1:n))/(n+1);
    E_copula(ii) = sum(u(ii)>=u(1:n) & v(ii)>=v(1:n))/n;
end

idxVec = 1:n;
RR = [idxVec' U E_copula Pi_copula ones(n,1)*-999]; % the last column will be filled in 
                                                    % by the computed domain

domainIdx = 1;
zz = RR(:,4)>=RR(:,5);
zz_diff = diff(zz);
I = find(zz_diff~=0);
startIdx = 1;
for ii=I'
    endIdx = ii;
    if(endIdx==startIdx)
        continue;
    end
    RR(startIdx:endIdx,6) = domainIdx;
    domainIdx = domainIdx + 1;
    startIdx = endIdx+1;
end
if(startIdx<=n)
    RR(startIdx:n,6) = domainIdx;
end

% compute the metric for each domain
numDomains = domainIdx;
domainMetrics = zeros(2,numDomains);
for ii=1:numDomains
    % get the data associated w/ this domain
    II = RR(:,6)==ii;
    domainData = RR(II,:);
    domainLen = length(find(II==1));
    numAvg = min(domainLen,3);
    
    % compute relative distance function
    if(mean(domainData(1:numAvg,4))>=mean(domainData(end-(numAvg-1):end,4)))
        maxECopVal = domainData(1,4);
        iMax = 1;
        minECopVal = domainData(end,4);
        iMin = size(domainData,1);
    else
        minECopVal = domainData(1,4);
        iMin = 1;
        maxECopVal = domainData(end,4);
        iMax = size(domainData,1);
    end
    
    piMinVal = domainData(iMin,end-1);
    piMaxVal = domainData(iMax,end-1);
    if(minECopVal>=piMinVal)
        MminVal = min(domainData(iMin,2:3));
        lambdaMin = (minECopVal-piMinVal)/(MminVal-piMinVal);
    else
        WminVal = max(1-d+sum(domainData(iMin,2:3)),0);
        lambdaMin = (minECopVal-piMinVal)/(WminVal-piMinVal);
    end
    if(maxECopVal>=piMaxVal)
        MmaxVal = min(domainData(iMax,2:3));
        lambdaMax = (maxECopVal-piMaxVal)/(MmaxVal-piMaxVal);
    else
        WmaxVal = max(1-d+sum(domainData(iMax,2:3)),0);
        lambdaMax = (maxECopVal-piMaxVal)/(WmaxVal-piMaxVal);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check if we are at a local minimum or global minimum, according to
    % the following rules: (from Mohsen's original paper/algorithm)
    % Calculate the absolute difference between the three consecutive values of Cn(u(i),vj) 
    % centered at u_i_min (respectively at u_i_max) and decide that the central point is a local optimum if 
    % (i) both absolute differences are smaller than or equal to 1/n and 
    % (ii) there are more than four points within the two adjacent domains, 𝔇𝑖 and 𝔇𝑖+1
    gamma = (lambdaMin+lambdaMax)/2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     gamma = abs(corr(domainData(:,2),domainData(:,3),'type','kendall'));
    domainMetrics(1,ii) = gamma;
    domainMetrics(2,ii) = domainLen;
end

cosValue = nanmean(domainMetrics(1,:));
% cosValue = sum(domainMetrics(1,:).*domainMetrics(2,:))/(n+numDomains-1);

end