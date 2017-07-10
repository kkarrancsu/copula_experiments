function [metric, rectangleCellOut] = ...
    cim_v8a(x, y, varargin)
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

% default values
minscanincr = 0.015625;

% overwrite defaults with user-inputted values
nVarargin = length(varargin);
switch nVarargin
    case 0
    otherwise
        minscanincr = varargin{1};
end

% convert X and Y to pseudo-observations, and scale to be between 0-1
M = length(x);
% u = pobs(x)*(M+1)/M;
% v = pobs(y)*(M+1)/M;
[u,v] = pobs_sorted(x,y,1);

axisCfgs = [1 2];
ax2minmaxCfgs = { {[0,1]}, {[0,0.5],[0.5,1]} };

% perform a scan pattern while varying U with V full-range, then swap the U-V axes
metrics = []; maxIICell = {};
rectangleAggr = {};  rectangleAggrIdx = 1;
for axisCfg=axisCfgs
    for ax2minmaxCfgsIdx=1:length(ax2minmaxCfgs)
        ax2minmaxCfg = ax2minmaxCfgs{ax2minmaxCfgsIdx};
        metricVecAggr = cell(1,length(ax2minmaxCfg));
        numPtsVecAggr = cell(1,length(ax2minmaxCfg));
        rectangleCellAggr = cell(1,length(ax2minmaxCfg));
        for ax2mmCfgIdx=1:length(ax2minmaxCfg)
            ax2mmCfg = ax2minmaxCfg{ax2mmCfgIdx};
            ax2min = ax2mmCfg(1);
            ax2max = ax2mmCfg(2);

            metricCell = {}; numPtsCell = {};
            rectanglesCell = {}; rectanglesCellIdx = 1;
            scanincr = 1;
            while(scanincr>=minscanincr)
                switch(axisCfg)
                    case 1
                        ax1pts = u; ax2pts = v;
                    case 2
                        ax1pts = v; ax2pts = u;
                end

                [metricVecTmp, numPtsVecTmp, rectangles] = ...
                    scanForDep(ax1pts,ax2pts,ax2min,ax2max,scanincr);
                
                metricCell{rectanglesCellIdx} = metricVecTmp;
                numPtsCell{rectanglesCellIdx} = numPtsVecTmp;
                rectanglesCell{rectanglesCellIdx} = rectangles; 
                
                rectanglesCellIdx = rectanglesCellIdx + 1;

                scanincr = scanincr/2;
            end
            metricVecAggr{ax2mmCfgIdx} = metricCell;
            numPtsVecAggr{ax2mmCfgIdx} = numPtsCell;
            rectangleCellAggr{ax2mmCfgIdx} = rectanglesCell;
        end
        % compute the metric for this.  putting stuff outside the 2nd
        % for-loop allows us to combine the results for {[0,1]} and
        % {[0,0.5][0.5,1]} easily.  metricVecAggr should have the
        % results for when ax2 is {[0,1]} and compute a metric, then it
        % should have the results for {[0,0.5],[0.5,1]} and compute a
        % metric for it.  at the end of processing, we compute a
        % maximum.
        [m,maxIIVec] = computeMetricFromAggregates(metricVecAggr, numPtsVecAggr);
        metrics = [metrics m]; maxIICell{rectangleAggrIdx} = maxIIVec;
        rectangleAggr{rectangleAggrIdx} = rectangleCellAggr;
        rectangleAggrIdx = rectangleAggrIdx + 1;
    end
end

[metric, metricMaxIdx] = max(metrics);
if(nargout>1)
    idx1 = metricMaxIdx;
    tmp = maxIICell{idx1};
    if(mod(metricMaxIdx,2)~=0)  % means full v-scale scan was best option
        idx2 = 1; idx3 = tmp(1);
        rectangleCellOut = rectangleAggr{idx1}{idx2}{idx3};
    else
        % WARNING -- need to change this along w/ the scan-patterns
        idx3_1 = tmp(1); idx3_2 = tmp(2);
        % disambiguate this to determine the number of regions!
        regionMat1 = rectangleAggr{idx1}{1}{idx3_1};
        regionMat2 = rectangleAggr{idx1}{2}{idx3_2};
        % TODO: in the future, think about fancy merging of the regions to
        % get a more accurate monotonicity count, for now, flatten the
        % region matrices and return as a first order approximation of the
        % monotonic regions
        rectangleCellOut = [regionMat1 regionMat2];
    end
end

end

function [metric, maxIIVec] = computeMetricFromAggregates(metricVecAggr, numPtsVecAggr)

metrics = zeros(2,length(metricVecAggr));
maxIIVec = [];

for jj=1:length(metricVecAggr)
    groupMetrics = metricVecAggr{jj};
    groupVecLens = numPtsVecAggr{jj};
    
    weightedMetric = -999;
    numPts = -999;
    for ii=1:length(groupMetrics)
        % compute the metric for each group
        gMetric = groupMetrics{ii};
        gVecLen = groupVecLens{ii};
        
        % compute the weighted metric
        weightedMetricCompute = sum( gMetric.*gVecLen/(sum(gVecLen)) );
        numPtsCompute = sum(gVecLen);
        
        % take the max
        if(weightedMetricCompute>weightedMetric)
            weightedMetric = weightedMetricCompute;
            numPts = numPtsCompute;
            maxIIVal = ii;
        end
    end
    maxIIVec = [maxIIVec maxIIVal];
    metrics(1,jj) = weightedMetric;
    metrics(2,jj) = numPts;
end

% combine group metrics
metric = sum( metrics(2,:)/sum(metrics(2,:)).*metrics(1,:) );

end

function [metricVec, numPtsVec, rectangles] = scanForDep(ax1pts, ax2pts, ax2min, ax2max, scanincr)
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
    ax1_match = find(ax1pts>ax1min & ax1pts<=ax1max);
    ax2_match = find(ax2pts>ax2min & ax2pts<=ax2max);
    matchIdxs = intersect(ax1_match,ax2_match);
    matchPts = [ax1pts(matchIdxs) ax2pts(matchIdxs)];
end