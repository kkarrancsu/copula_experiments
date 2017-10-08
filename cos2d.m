function [cosValue, RR] = cos2d(x,y)
%COS3D - computes the 2-D cos metric
n = length(x);

X = [x y]; 
U = pobs(X);
U = sortrows(U);
u = U(:,1); v = U(:,2);

E_copula = zeros(n,1);
for ii=1:n
    E_copula(ii) = sum(u(ii)>=u(1:n) & v(ii)>=v(1:n))/(n+1);
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
    while(eCopulaSorted)
        % get all points which haven't been searched yet.
        I = RR(:,end)==-999;
        SS = RR(I,:);
        % find all points and rank in order of euclidean distance to the point
        % available
        D = pdist2(RR(curIdx,2:3),SS(:,2:3));
        [~,II] = sort(D);
        
        % search this list until we get to where our E_copula value is no
        % longer sorted
        noDirectionGood = 1;
        for ii=II
            SS(ii,:)
            rrIdx = SS(ii,1);
            eCopVal = RR(rrIdx,end-1);
            tryVec1 = [eCopulaVec eCopVal];
            tryVec2 = [eCopVal eCopulaVec];
            if(issorted(tryVec1))
                eCopulaVec = tryVec1;
                RR(rrIdx,end) = domainIdx;
                noDirectionGood = 0;
            elseif(issorted(tryVec2))
                eCopulaVec = tryVec2;
                RR(rrIdx,end) = domainIdx;
                noDirectionGood = 0;
            end
        end

        % repeat until we have no direction in which the E_copula is sorted
        if(noDirectionGood || (sum(RR(:,end)>0)==n) )
            eCopulaSorted = 0;
        end
    end
    domainIdx = domainIdx + 1;
end

% compute the metric for each domain
numDomains = domainIdx-1
domainMetrics = zeros(1,numDomains);
for ii=1:numDomains
    % get the data associated w/ this domain
    II = RR(:,5)==ii;
    domainData = RR(II,:);
    
    % compute relative distance function
    [minECopVal, iMin] = min(domainData(:,4));
    [maxECopVal, iMax] = max(domainData(:,4));
    piMinVal = prod(domainData(iMin,2:3));
    piMaxVal = prod(domainData(iMax,2:3));
    if(minECopVal>=piMinVal)
        MminVal = min(domainData(iMin,2:3));
        lambdaMin = (minECopVal-piMinVal)/(MminVal-piMinVal);
    else
        WminVal = max(1-3+sum(domainData(iMin,2:3)),0);
        lambdaMin = (minECopVal-piMinVal)/(WminVal-piMinVal);
    end
    if(maxECopVal>=piMaxVal)
        MmaxVal = min(domainData(iMax,2:3));
        lambdaMax = (maxECopVal-piMaxVal)/(MmaxVal-piMaxVal);
    else
        WmaxVal = max(1-3+sum(domainData(iMax,2:3)),0);
        lambdaMax = (maxECopVal-piMaxVal)/(WmaxVal-piMaxVal);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TODO: check if we are at a local minimum or global minimum
    gamma = (lambdaMin+lambdaMax)/2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    domainMetrics(ii) = gamma;
end

cosValue = mean(domainMetrics);

end