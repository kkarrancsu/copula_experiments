function [metric] = cim_v2(x,y,minscanincr)
% WARNING -- ONLY DOES FULL SCAN!

if(nargin<3)
    minscanincr = 0.025;
end

[u,v] = pobs_sorted(x,y,1);
M = length(u);

signature = cim_region_finder(x, y, minscanincr);
zscVec = zeros(1,length(signature.keys()));
cimVec = zeros(1,length(signature.keys()));
jj = 1;
for key=signature.keys()
    scanincr = key{1};
    zz = signature(scanincr);
    tau = zz.metricVec;
    n = zz.numPtsVec;
    den = (2*(2*n+5))./(9*n.*(n-1));
    zScoreProfile = tau./den;
    
    regionIdxs = [findruns(zScoreProfile) length(zScoreProfile)];
    
    % compute the CIM metric
    if(~isempty(regionIdxs))
        minU = 0-eps; maxU = regionIdxs(1)*scanincr;
        pts = getBoundPoints(u,v,minU,maxU);
        nn = size(pts,1);
        cimTotalVal = 0; 
        zTotalVal = 0;
        cimVal = abs(corr(pts(:,1),pts(:,2),'type','kendall'));
        zVal = cimVal/( (2*(2*nn+5))/(9*nn*(nn-1)) );
        w = nn/M;
        cimTotalVal = cimTotalVal + w*cimVal;
        zTotalVal = zTotalVal + w*zVal;
        for ii=2:length(regionIdxs)
            minU = regionIdxs(ii-1)*scanincr;
            maxU = regionIdxs(ii)*scanincr;
            pts = getBoundPoints(u,v,minU,maxU);
            nn = length(pts);
            cimVal = abs(corr(pts(:,1),pts(:,2),'type','kendall'));
            zVal = cimVal/( (2*(2*nn+5))/(9*nn*(nn-1)) );
            w = nn/M;
            cimTotalVal = cimTotalVal + w*cimVal;
            zTotalVal = zTotalVal + w*zVal;
        end
    else
        cimTotalVal = abs(corr(x,y,'type','kendall'));
        nn = n(end);
        zTotalVal = cimTotalVal/( (2*(2*nn+5))/(9*nn*(nn-1)) );
    end
    
    zscVec(jj) = zTotalVal;
    cimVec(jj) = cimTotalVal;
    jj = jj + 1;
end

% [~,maxI] = max(zscVec);
% metric = cimVec(maxI);
metric = max(cimVec);

end


function [matchPts] = getBoundPoints(ax1pts, ax2pts, ax1min, ax1max)
    ax1_match = find(ax1pts>ax1min & ax1pts<=ax1max);
    matchPts = [ax1pts(ax1_match) ax2pts(ax1_match)];
end