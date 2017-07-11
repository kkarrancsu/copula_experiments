function [metric] = cim_v8a_rev3acc(x, y, minScanIncr)
% in this version, this is the same as cim_v8a_cc.m, except the scanForDep
% function is copied from cim_v8a.m

% convert X and Y to pseudo-observations, and scale to be between 0-1
[u,v] = pobs_sorted_cc(x,y);
MAX_NUM_RECT = ceil(length(x)/2);

axisCfgs = [1 2];
ax2minmaxCfgs = { {[0,1]}, {[0,0.5],[0.5,1]} };

% perform a scan pattern while varying U with V full-range, then swap the U-V axes
vecLen = length(axisCfgs)*length(ax2minmaxCfgs);
numScans = ceil(log2(1/minScanIncr))+1;
% scanPattern = [1,0.5,0.25,0.125,0.0625,0.03125,0.015625];
% numScans = length(scanPattern);

metricCell = zeros(numScans,MAX_NUM_RECT); numPtsCell = zeros(numScans,MAX_NUM_RECT);
numRectanglesCreatedVec = zeros(numScans);

% pre-allocate to max-length for matlab coder speed purposes
metricVecAggr = cell(1,2);
numPtsVecAggr = cell(1,2);
numRectanglesVecAggr = cell(1,2); 

% assign cell elements empty stuff to make matlab-coder happy :x
for ii=1:2
    metricVecAggr{ii} = metricCell;
    numPtsVecAggr{ii} = numPtsCell;
    numRectanglesVecAggr{ii} = numRectanglesCreatedVec;
end

metrics = zeros(1,vecLen); 
rectangleAggrIdx = 1;
for axisCfg=axisCfgs
    for ax2minmaxCfgsIdx=1:length(ax2minmaxCfgs)
        ax2minmaxCfg = ax2minmaxCfgs{ax2minmaxCfgsIdx};
        
        ax2minmaxCfgLen = length(ax2minmaxCfg);
        
        for ax2mmCfgIdx=1:ax2minmaxCfgLen
            ax2mmCfg = ax2minmaxCfg{ax2mmCfgIdx};
            ax2min = ax2mmCfg(1);
            ax2max = ax2mmCfg(2);

            scanincr = 1;
            for zz=1:numScans
%                 scanincr = scanPattern(zz);
                switch(axisCfg)
                    case 1
                        ax1pts = u; ax2pts = v;
                    otherwise  % changed from case 2 to otherwise for matlab coder
                        ax1pts = v; ax2pts = u;
                end

                [metricVecTmp, numPtsVecTmp, numRectanglesCreated] = ...
                    scanForDep(ax1pts,ax2pts,ax2min,ax2max,scanincr);
                
                % pad, so that we can conform to api :(
                metricVecTmpPadded = [metricVecTmp zeros(1,MAX_NUM_RECT-length(metricVecTmp))];
                numPtsVecTmpPadded = [numPtsVecTmp zeros(1,MAX_NUM_RECT-length(numPtsVecTmp))];
                
                metricCell(zz,:) = metricVecTmpPadded;  
                numPtsCell(zz,:) = numPtsVecTmpPadded;
                numRectanglesCreatedVec(zz) = numRectanglesCreated;
                
                scanincr = scanincr/2;
            end
            metricVecAggr{ax2mmCfgIdx} = metricCell;
            numPtsVecAggr{ax2mmCfgIdx} = numPtsCell;
            numRectanglesVecAggr{ax2mmCfgIdx} = numRectanglesCreatedVec;
        end
        % compute the metric for this.  putting stuff outside the 2nd
        % for-loop allows us to combine the results for {[0,1]} and
        % {[0,0.5][0.5,1]} easily.  metricVecAggr should have the
        % results for when ax2 is {[0,1]} and compute a metric, then it
        % should have the results for {[0,0.5],[0.5,1]} and compute a
        % metric for it.  at the end of processing, we compute a
        % maximum.
        m = computeMetricFromAggregates(metricVecAggr, numPtsVecAggr, ax2minmaxCfgLen, numRectanglesVecAggr);
        metrics(rectangleAggrIdx) = m;
        rectangleAggrIdx = rectangleAggrIdx + 1;
    end
end

metric = max(metrics);

end

function [metric] = computeMetricFromAggregates(metricVecAggr, numPtsVecAggr, ax2minmaxCfgLen, numRectanglesAggr)
coder.inline('always');
metrics = zeros(2,ax2minmaxCfgLen);

for jj=1:ax2minmaxCfgLen
    groupMetrics = metricVecAggr{jj};
    groupVecLens = numPtsVecAggr{jj};
    groupRectangles = numRectanglesAggr{jj};
    
    weightedMetric = -999;
    numPts = -999;
    for ii=1:length(groupRectangles)
        % compute the metric for each group
        numRects = groupRectangles(ii);
        gMetric = groupMetrics(ii,1:numRects);
        gVecLen = groupVecLens(ii,1:numRects);
        
        % compute the weighted metric
        weightedMetricCompute = sum( gMetric.*gVecLen/(sum(gVecLen)) );
        numPtsCompute = sum(gVecLen);
        
        % take the max
        if(weightedMetricCompute>weightedMetric)
            weightedMetric = weightedMetricCompute;
            numPts = numPtsCompute;
        end
    end
    metrics(1,jj) = weightedMetric;
    metrics(2,jj) = numPts;
end

% combine group metrics
metric = sum( metrics(2,:)/sum(metrics(2,:)).*metrics(1,:) );

end

% the only modification to this function was to return the rectanglesIdx as
% the 3rd output argument, rather than the rectangles themselves, b/c we
% are only caring about the metric right now, not the actual regions.
function [metricVec, numPtsVec, rectanglesIdx] = scanForDep(ax1pts, ax2pts, ax2min, ax2max, scanincr)
%scanForDep - scans for dependencies across the first axis (if you would
%like to scan across the second axis, simply swap the input arguments to 
%this function).

ax1min = 0; ax1max = scanincr;
newRectangle = 1;
metricVec = [];
numPtsVec = [];
rectangles = []; rectanglesIdx = 1;
numStdDev = 4;
while ax1max<=1
    % find all the points which are contained within this cover rectangle
    matchPts = getPointsWithinBounds(ax1pts, ax2pts, ax1min, ax1max, ax2min, ax2max);
    
    numPts = size(matchPts,1);
    if(numPts>=2)   % make sure we have enough points to compute the metric
        % compute the concordance
        metricRectangle = abs(taukl( matchPts(:,1),matchPts(:,2)));
        stdTau = ((1-metricRectangle)*sqrt( (2*(2*numPts+5))/(9*numPts*(numPts-1)) ) )*numStdDev;
        if(newRectangle)
            newRectangle = 0;
        else
            % compare to the previous concordance, if there is a change by the
            % threshold amount, rewind the axes of the cover rectangle and 
            if( (metricRectangle < (metricRectanglePrev-stdTau)) )
                metricVec = [metricVec metricRectanglePrev];
                numPtsVec = [numPtsVec numPtsPrev];
                rectangles(:,rectanglesIdx) = [ax1min ax1max-scanincr ax2min ax2max]; rectanglesIdx = rectanglesIdx + 1;
                % start the new cover rectangle
                ax1min = ax1max - scanincr;
                ax1max = ax1min;        % it will be incremented below
                newRectangle = 1;
            end
        end

        metricRectanglePrev = metricRectangle;
        numPtsPrev = numPts;
    end
    ax1max = ax1max + scanincr;
    
    if(ax1max>1)
        if(exist('metricRectanglePrev','var'))
            metricVec = [metricVec metricRectanglePrev];
            numPtsVec = [numPtsVec length(matchPts)];
            rectangles(:,rectanglesIdx) = [ax1min 1 ax2min ax2max]; rectanglesIdx = rectanglesIdx + 1;
        end
    end
end

end

function [matchPts, matchIdxs] = getPointsWithinBounds(ax1pts, ax2pts, ax1min, ax1max, ax2min, ax2max)
coder.inline('always');
ax1_match = find(ax1pts>ax1min & ax1pts<=ax1max);
ax2_match = find(ax2pts>ax2min & ax2pts<=ax2max);
matchIdxs = intersect(ax1_match,ax2_match);
matchPts = [ax1pts(matchIdxs) ax2pts(matchIdxs)];
end