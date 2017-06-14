function [borderPts] = rd3(ax1pts, ax2pts, scanincr, diffthresh)
%scanForDep - scans for dependencies across the first axis (if you would
%like to scan across the second axis, simply swap the input arguments to 
%this function).

ax1min = 0; ax1max = scanincr; ax2min = 0; ax2max = 1;
newRectangle = 1;
borderPts = [];
while ax1max<=1
    % find all the points which are contained within this cover rectangle
    matchPts = getPointsWithinBounds(ax1pts, ax2pts, ax1min, ax1max, ax2min, ax2max);
    
    numPts = size(matchPts,1);
    if(numPts>=2)   % make sure we have enough points to compute the metric
        % compute the concordance
        metricRectangle = abs(taukl( matchPts(:,1),matchPts(:,2) ));
        zsc  = metricRectangle./sqrt( (2*(2*numPts+5))./(9*numPts.*(numPts-1)) );   

        if(newRectangle)
            newRectangle = 0;
        else
            % compare to the previous concordance, if there is a change by the
            % threshold amount, rewind the axes of the cover rectangle and 
            percentageChange = (metricRectangle-metricRectanglePrev)/metricRectanglePrev*100;
            diffthreshAdaptive = (1/zsc)*diffthresh;
            if(percentageChange<(-1*diffthreshAdaptive))
                borderPts = [borderPts ax1max-scanincr];
                % start the new cover rectangle
                ax1min = ax1max - scanincr;
                ax1max = ax1min;        % it will be incremented below
                newRectangle = 1;
            end
        end

        metricRectanglePrev = metricRectangle;
    end
    ax1max = ax1max + scanincr;
    
    if(ax1max>1)
        if(exist('metricRectanglePrev','var'))
            borderPts = [borderPts 1];
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