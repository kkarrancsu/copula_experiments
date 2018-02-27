function [cosValue, RR] = cos2d_v3(x,y)
%COS2D_v3 - computes the 2-D cos metric
n = length(x);
d = 2;

MIN_TOL = 0.001;

% X = [x y]; 
% U = pobs(X);
% U = sortrows(U);
% u = U(:,1); v = U(:,2);
% [u,v] = pobs_sorted_cc(x,y); u = u*n/(n+1); v = v*n/(n+1);

u = pobs(x);
v = pobs(y);
[u,I] = sort(u);
v = v(I);

U = [u v];

E_copula = zeros(n,1);
Pi_copula = u.*v;
for ii=1:n
    E_copula(ii) = sum(u(ii)>=u(1:n) & v(ii)>=v(1:n))/(n+1);
%     E_copula(ii) = sum(u(ii)>=u(1:n) & v(ii)>=v(1:n))/n;
end

idxVec = 1:n;
RR = [idxVec' U E_copula Pi_copula ones(n,1)*-999]; % the last column will be filled in 
                                                    % by the computed domain

E_COPULA_IDX = 4;
PI_COPULA_IDX = 5;
DOMAIN_IDX = 6;
                                                    
domainIdx = 1;
zz = RR(:,E_COPULA_IDX)>=RR(:,PI_COPULA_IDX);
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

% compute the metric for each domain
numDomains = domainIdx;
domainMetrics = zeros(2,numDomains);
for ii=1:numDomains
    % get the data associated w/ this domain
    II = RR(:,6)==ii;
    domainData = RR(II,:);
    domainLen = length(find(II==1));
%     numAvg = min(domainLen,3);
    
    % compute relative distance function
%     if(mean(domainData(1:numAvg,4))>=mean(domainData(end-(numAvg-1):end,4)))
%         maxECopVal = domainData(1,4);
%         iMax = 1;
%         minECopVal = domainData(end,4);
%         iMin = size(domainData,1);
%     else
%         minECopVal = domainData(1,4);
%         iMin = 1;
%         maxECopVal = domainData(end,4);
%         iMax = size(domainData,1);
%     end
    [minECopVal,iMin] = min(domainData(:,E_COPULA_IDX));
    [maxECopVal,iMax] = max(domainData(:,E_COPULA_IDX));
    
    piMinVal = domainData(iMin,PI_COPULA_IDX);
    piMaxVal = domainData(iMax,PI_COPULA_IDX);
    MminVal = min(domainData(iMin,2:3));
    WminVal = max(1-d+sum(domainData(iMin,2:3)),0);
    minNumerator = (minECopVal-piMinVal);
    if(minECopVal>=piMinVal)
        minDenominator = (MminVal-piMinVal);
    else
        minDenominator = (WminVal-piMinVal); 
    end
    if(abs(minDenominator)<MIN_TOL && abs(minNumerator)<MIN_TOL)
        lambdaMin = nan;
    elseif(abs(minDenominator)<MIN_TOL || abs(minNumerator)<MIN_TOL)
        lambdaMin = 0;
    else
        lambdaMin = min(minNumerator/minDenominator,1);
    end
    
    MmaxVal = min(domainData(iMax,2:3));
    WmaxVal = max(1-d+sum(domainData(iMax,2:3)),0);
    maxNumerator = (maxECopVal-piMaxVal);
    if(maxECopVal>=piMaxVal)
        maxDenominator = (MmaxVal-piMaxVal);
    else
        maxDenominator = (WmaxVal-piMaxVal); 
    end
    
    if(abs(maxDenominator)<MIN_TOL && abs(maxNumerator)<MIN_TOL)
        lambdaMax = nan;
    elseif(abs(maxDenominator)<MIN_TOL || abs(maxNumerator)<MIN_TOL)
        lambdaMax = 0;
    else
        lambdaMax = min(maxNumerator/maxDenominator,1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check if we are at a local optimum or global optimum, according to
    % the following rules: (from Mohsen's original paper/algorithm)
    % Calculate the absolute difference between the three consecutive values of Cn(u(i),vj) 
    % centered at u_i_min (respectively at u_i_max) and decide that the central point is a local optimum if 
    % (i) both absolute differences are smaller than or equal to 1/n and 
    % (ii) there are more than four points within the two adjacent domains,
    % Di and Di+1
    if(isnan(lambdaMin) && isnan(lambdaMax))
        gamma = 0;
    elseif(isnan(lambdaMin) && ~isnan(lambdaMax))
        gamma = lambdaMax;
    elseif(~isnan(lambdaMin) && isnan(lambdaMax))
        gamma = lambdaMin;
    else
        gamma = (lambdaMin+lambdaMax)/2;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    domainMetrics(1,ii) = gamma;
    domainMetrics(2,ii) = domainLen;
    
%     fprintf('**************************************\n');
%     fprintf('[u_i_min,v_i_min]=[%0.04f,%0.04f] [u_i_max,v_i_max]=[%0.04f,%0.04f]\n',...
%         domainData(iMin,2),domainData(iMin,3),...
%         domainData(iMax,2),domainData(iMax,3));
%     fprintf('i_min >> C[u,v]=%0.04f Pi[u,v]=%0.04f M[u,v]=%0.04f W[u,v]=%0.04f\n',...
%         minECopVal,piMinVal,MminVal,WminVal);
%     fprintf('i_max >> C[u,v]=%0.04f Pi[u,v]=%0.04f M[u,v]=%0.04f W[u,v]=%0.04f\n',...
%         maxECopVal,piMaxVal,MmaxVal,WmaxVal);
%     fprintf('lambda_min=%0.04f lambda_max=%0.04f gamma=%0.04f\n',...
%         lambdaMin,lambdaMax,gamma);
%     fprintf('**************************************\n\n');
    
end
% domainMetrics
cosValue = sum(domainMetrics(1,:).*domainMetrics(2,:))/(n+numDomains-1);

end