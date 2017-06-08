function [signature] = cim_region_finder(x, y, minscanincr)

% minscanincr should be 0.025
% mode -- if 0, we perform a sliding window (non-aggregating)
% mode -- if 1, we perform a sliding window (aggregating)

% convert X and Y to pseudo-observations, and scale to be between 0-1
% M = length(x);
% u = pobs(x)*(M+1)/M;
% v = pobs(y)*(M+1)/M;
% sort for computational efficiency
[u,v] = pobs_sorted(x,y,1);     % the 1 is a scale flag

scanincr = 1;
signature = containers.Map('KeyType','double','ValueType','any');
while(scanincr>=minscanincr)
    [metricVec,numPtsVec] = scanForDep(u,v,scanincr);
    data = struct;
    data.metricVec = metricVec;
    data.numPtsVec = numPtsVec;
    signature(scanincr) = data;
    scanincr = scanincr/2;
end

end

function [metricVec, numPtsVec] = scanForDep(ax1pts, ax2pts, scanincr)
%scanForDep - scans for dependencies across the first axis (if you would
%like to scan across the second axis, simply swap the input arguments to 
%this function).
M = length(ax1pts);
ax1min = 0; ax1max = scanincr; ax2min = 0; ax2max = 1;
numElem = 1/scanincr + 1;
metricVec = zeros(1,numElem);
numPtsVec = zeros(1,numElem);
elemIdx = 2;
while ax1max<=1
    % find all the points which are contained within this cover rectangle
    matchPts = getPointsWithinBounds(ax1pts, ax2pts, ax1min, ax1max, ax2min, ax2max,M);
    numPts = size(matchPts,1);
    if(numPts>=2)   % make sure we have enough points to compute the metric
        % compute the concordance
        metricRectangle = abs(corr(matchPts(:,1),matchPts(:,2),'type','kendall'));
        metricVec(elemIdx) = metricRectangle;
        numPtsVec(elemIdx) = numPts;
        elemIdx = elemIdx + 1;
    end
    ax1max = ax1max + scanincr;
    
    % if we didn't exactly match the bounds
    if((ax1max-scanincr)>1)
        matchPts = getPointsWithinBounds(ax1pts, ax2pts, ax1min, 1, ax2min, ax2max, M);
        numPts = size(matchPts,1);
        if(numPts>=2)   % make sure we have enough points to compute the metric
            % compute the concordance
            metricRectangle = abs(taukl( matchPts(:,1),matchPts(:,2) ));
            metricVec(elemIdx) = metricRectangle;
            numPtsVec(elemIdx) = numPts;
            elemIdx = elemIdx + 1;
        end
    end
end

metricVec(elemIdx:end) = [];
numPtsVec(elemIdx:end) = [];

end

function [matchPts, matchIdxs] = getPointsWithinBounds(ax1pts, ax2pts, ax1min, ax1max, ax2min, ax2max, M)
%     ax1_match = find(ax1pts>ax1min & ax1pts<=ax1max);
    ax1MinIdx = ax1min*M+1; ax1MaxIdx = ax1max*M;
    ax1_match = ax1MinIdx:ax1MaxIdx;
    if(ax2min==0 && ax2max==1)
        matchIdxs = ax1_match;
    else
        ax2_match = find(ax2pts>ax2min & ax2pts<=ax2max);
        matchIdxs = intersect(ax1_match,ax2_match);
    end
    matchPts = [ax1pts(matchIdxs) ax2pts(matchIdxs)];
end