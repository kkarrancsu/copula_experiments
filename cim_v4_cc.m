function [metric] = cim_v4_cc(x, y)
%CIM - Copula Index for Detecting Dependence and Monotonicity between
%Stochastic Signals.  See associated paper... to be published and preprint
%located here: 
% Inputs:
%  x - the x variable
%  y - the y variable
%  varargin{1} - minscanincr - the minimum scanning increment.  Large
%                              values will filter out high frequency
%                              dependencies, small values decrease the
%                              statistical power of the dependency metric
%  varargin{2} - diffthresh  - the threshold at which a change in
%                              concordance amount is detected.  Larger
%                              values are more robust to noise, but tend to
%                              miss high frequency changes.
%  varargin{3} - alpha       - the value used to determine significance
%                              level of a box's concordance level
% Outputs:
%  metric - the calculated dependency metric between x and y
%  resid  - the residual between the estimated concordance boxes and the
%           observed statistical variables.  Each concordance box's
%           residuals are provided separately
%  residAssocIdxs - the indices of the independent variable associated with
%                   each residual point, this is used by rscdm for residual
%                   alignment.
%
%**************************************************************************
%*                                                                        *
%* Copyright (C) 2017  Kiran Karra <kiran.karra@gmail.com>                *
%*                                                                        *
%* This program is free software: you can redistribute it and/or modify   *
%* it under the terms of the GNU General Public License as published by   *
%* the Free Software Foundation, either version 3 of the License, or      *
%* (at your option) any later version.                                    *
%*                                                                        *
%* This program is distributed in the hope that it will be useful,        *
%* but WITHOUT ANY WARRANTY; without even the implied warranty of         *
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
%* GNU General Public License for more details.                           *
%*                                                                        *
%* You should have received a copy of the GNU General Public License      *
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
%*                                                                        *
%**************************************************************************

% convert X and Y to pseudo-observations, and scale to be between 0-1
[u,v] = pobs_sorted_cc(x,y);
MAX_NUM_RECT = ceil(length(x)/2);

axisCfgs = [1 2];
ax2minmaxCfgs = { {[0,1]}, {[0,0.5],[0.5,1]} };

% perform a scan pattern while varying U with V full-range, then swap the U-V axes
vecLen = length(axisCfgs)*length(ax2minmaxCfgs);
scanPattern = [1,0.5,0.25,0.125,0.0625,0.03125,0.015625];
scanVecLen = length(scanPattern);

metricCell = zeros(scanVecLen,MAX_NUM_RECT); numPtsCell = zeros(scanVecLen,MAX_NUM_RECT);
numRectanglesCreatedVec = zeros(scanVecLen);

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

            for zz=1:scanVecLen
                scanincr = scanPattern(zz);
                switch(axisCfg)
                    case 1
                        ax1pts = u; ax2pts = v;
                    otherwise  % changed from case 2 to otherwise for matlab coder
                        ax1pts = v; ax2pts = u;
                end

                [metricVecTmp, numPtsVecTmp, numRectanglesCreated] = ...
                    scanForDep(ax1pts,ax2pts,ax2min,ax2max,scanincr,MAX_NUM_RECT);
                
                metricCell(zz,:) = metricVecTmp;
                numPtsCell(zz,:) = numPtsVecTmp;
                numRectanglesCreatedVec(zz) = numRectanglesCreated;
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

function [metricVec, numPtsVec, rectanglesIdx] = scanForDep(ax1pts, ax2pts, ax2min, ax2max, scanincr, maxNumRect)
%scanForDep - scans for dependencies across the first axis (if you would
%like to scan across the second axis, simply swap the input arguments to 
%this function).

ax1min = 0; ax1max = scanincr;
newRectangle = 1;

metricVec = zeros(1,maxNumRect);
numPtsVec = zeros(1,maxNumRect);
rectanglesIdx = 1;

metricRectanglePrev = -999;
numPtsPrev = 1;  % should get overwritten
while ax1max<=1
    % find all the points which are contained within this cover rectangle
    matchPts = getPointsWithinBounds(ax1pts, ax2pts, ax1min, ax1max, ax2min, ax2max);
    
    numPts = size(matchPts,1);
    if(numPts>=2)   % make sure we have enough points to compute the metric
        % compute the concordance
        metricRectangle = abs(taukl_cc( matchPts(:,1),matchPts(:,2)));
        zsc  = metricRectangle./sqrt( (2*(2*numPts+5))./(9*numPts.*(numPts-1)) );   
        if(newRectangle)
            newRectangle = 0;
        else
            % compare to the previous concordance, if there is a change by the
            % threshold amount, rewind the axes of the cover rectangle and 
            percentageChange = (metricRectangle-metricRectanglePrev)/metricRectanglePrev;
            diffthreshAdaptive = (1/zsc);
            if(percentageChange<(-1*diffthreshAdaptive))
                metricVec(rectanglesIdx) = metricRectanglePrev;
                numPtsVec(rectanglesIdx) = numPtsPrev;
                rectanglesIdx = rectanglesIdx + 1;
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
        if(metricRectanglePrev>=0)
            metricVec(rectanglesIdx) = metricRectanglePrev;
            numPtsVec(rectanglesIdx) = length(matchPts);
%             rectanglesIdx = rectanglesIdx + 1;
        end
    end
end

end

function [matchPts, matchIdxs] = getPointsWithinBounds(ax1pts, ax2pts, ax1min, ax1max, ax2min, ax2max)
    ax1_match = find(ax1pts>ax1min & ax1pts<=ax1max);
    ax2_match = find(ax2pts>ax2min & ax2pts<=ax2max);
    matchIdxs = intersect(ax1_match,ax2_match);
    matchPts = [ax1pts(matchIdxs) ax2pts(matchIdxs)];
end