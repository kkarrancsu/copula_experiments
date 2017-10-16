function [cosValue, RR] = cos2d_v2(x,y)
%COS3D - computes the 2-D cos metric
n = length(x);
d = 2;

X = [x y]; 
U = pobs(X);
U = sortrows(U);
u = U(:,1); v = U(:,2);

E_copula = zeros(n,1);
for ii=1:n
%     E_copula(ii) = sum(u(ii)>=u(1:n) & v(ii)>=v(1:n))/(n+1);
    E_copula(ii) = sum(u(ii)>=u(1:n) & v(ii)>=v(1:n))/n;
end

% plot(u,E_copula)
% plot(E_copula)

idxVec = 1:n;
RR = [idxVec' U E_copula ones(n,1)*-999]; % the last column will be filled in by the computed
                                              % domain

domainIdx = 1;
while(sum(RR(:,end)>0)<n)  % condition checks to make sure we've assigned a domain to every 
                         % possible point
    % find the first point which doesn't have a domain associated with it
    % yet
    curIdx = find(RR(:,end)==-999,1);
    % assign it the next available domain
    RR(curIdx,5) = domainIdx;
    eCopulaVec = RR(curIdx,end-1);
    eCopulaSorted = 1;
    while(eCopulaSorted && (sum(RR(:,end)>0)<n))
        % get all points which haven't been searched yet.
        I = RR(:,end)==-999;
        SS = RR(I,:);
        % find all points and rank in order of euclidean distance to the point
        % available
        D = pdist2(RR(curIdx,2:3),SS(:,2:3));
        [~,ii] = min(D);
        
        rrIdx = SS(ii,1);
        eCopVal = RR(rrIdx,end-1);
        tryVec = [eCopulaVec eCopVal];
        if(issorted(tryVec,'monotonic'))
            eCopulaVec = tryVec;
            RR(rrIdx,end) = domainIdx;
            curIdx = find(RR(:,end)==-999,1);
        else
            eCopulaSorted = 0;
        end
    end
    domainIdx = domainIdx + 1;
end

% compute the metric for each domain
numDomains = domainIdx-1;
domainMetrics = zeros(2,numDomains);
for ii=1:numDomains
    % get the data associated w/ this domain
    II = RR(:,5)==ii;
    domainData = RR(II,:);
    domainLen = length(find(II==1));
    
    % compute relative distance function
    [minECopVal, iMin] = min(domainData(:,end-1));
    [maxECopVal, iMax] = max(domainData(:,end-1));
    piMinVal = prod(domainData(iMin,2:3));
    piMaxVal = prod(domainData(iMax,2:3));
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
    isLocalOptimum = 0;
    
    if(domainLen>=3)
        % perform the tests
        iiMin = min(iMin,iMax);
        iiMax = max(iMin,iMax);
        if((iiMax-iiMin)>=2)  % make sure we have enough points to conduct the tests below
            loEcopDiffCheckVec = diff(RR(iiMin:iiMin+2,4))<=(1/n);
            hiEcopDiffCheckVec = diff(RR(iiMax-2:iiMax,4))<=(1/n);
            if(sum(loEcopDiffCheckVec)==2 && sum(hiEcopDiffCheckVec)==2)
                isLocalOptimum = 1;
            end
            if(isLocalOptimum && domainLen<=4)
                if(ii<numDomains)
                    % count # of samples in next domain
                    II2 = RR(:,5)==(ii+1);
                    domain2Len = length(find(II2==1));
                    if((domainLen+domain2Len)<=4)
                        isLocalOptimum = 0;
                    end
                end
            end
        end
    else
        isLocalOptimum = 0;
    end
    if(isLocalOptimum)
        gamma = 1;
    else
        gamma = (lambdaMin+lambdaMax)/2;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    domainMetrics(1,ii) = gamma;
    domainMetrics(2,ii) = domainLen;
end

% cosValue = mean(domainMetrics(1,:));
cosValue = sum(domainMetrics(1,:).*domainMetrics(2,:))/(n+numDomains-1);

end