function [signature] = probabilistic_region_finder(x, y)

% sort for computational efficiency
[u,v] = pobs_sorted(x,y,1);     % the 1 is a scale flag
scanincrVec = 10:10:length(u);
signature = containers.Map('KeyType','double','ValueType','any');

for ii=1:length(scanincrVec)
    scanincr = scanincrVec(ii);
    signature(scanincr) = scanForDep(u,v,scanincr);
end

end

function [distStats] = scanForDep(ax1pts, ax2pts, scanincr)
% scanincr is the # of points!

M = length(ax1pts);
ax1StartIdx = 1; ax1EndIdx = scanincr;

numLoopIter = ceil(M/scanincr);
distStats = zeros(11,numLoopIter);   % (1,:)  -> theoretical max overlap between distribution(R1,R2)
                                     % (2,:)  -> actual overlap between distribution(R1,R2)
                                     % (3,:)  -> theoretical max overlap between distribution(R3,R2)
                                     % (4,:)  -> actual overlap between distribution(R3,R2)
                                     % (5,:)  -> theoretical max overlap between distribution(R1,R3)
                                     % (6,:)  -> actual overlap between distribution(R1,R3)
                                     % (7,:)  -> mu1
                                     % (8,:)  -> sigma1
                                     % (9,:)  -> mu2
                                     % (10,:) -> sigma2
                                     % (11,:) -> mu3
                                     % (12,:) -> sigma3

for ii=1:numLoopIter
    % compute distributions for R1,R2, and R3
    [mu1,sigma1] = getGroupStatistics(ax1pts,ax2pts,ax1StartIdx,ax1EndIdx);
    [mu2,sigma2] = getGroupStatistics(ax1pts,ax2pts,ax1StartIdx+scanincr,ax1EndIdx+scanincr);
    [mu3,sigma3] = getGroupStatistics(ax1pts,ax2pts,ax1StartIdx,ax1EndIdx+scanincr);

    distStats(1,ii) = computeOvlpProb(mu1,sigma1,mu1,sigma2);
    distStats(2,ii) = computeOvlpProb(mu1,sigma1,mu2,sigma2);
    distStats(3,ii) = computeOvlpProb(mu3,sigma3,mu3,sigma2);
    distStats(4,ii) = computeOvlpProb(mu3,sigma3,mu2,sigma2);
    distStats(5,ii) = computeOvlpProb(mu1,sigma1,mu1,sigma3);
    distStats(6,ii) = computeOvlpProb(mu1,sigma1,mu3,sigma3);
    distStats(7,ii)  = mu1; distStats(8,ii)  = sigma1;
    distStats(9,ii)  = mu2; distStats(10,ii) = sigma2;
    distStats(11,ii) = mu3; distStats(12,ii) = sigma3;
end

end

function [mu,sigma] = getGroupStatistics(ax1pts,ax2pts,ax1StartIdx,ax1EndIdx)
    if(ax1StartIdx>length(ax1pts))
        % means we've exceeded the bounds
        mu = 0; sigma = 1;
    else
        if(ax1EndIdx>length(ax1pts))
            ax1EndIdx = length(ax1pts);
        end
        uu = ax1pts(ax1StartIdx:ax1EndIdx); vv = ax2pts(ax1StartIdx:ax1EndIdx);
        n = ax1EndIdx-ax1StartIdx+1;
        mu = corr(uu,vv,'type','kendall'); sigma = getVariance(n);
    end
end

function [sigma] = getVariance(numPts)
    % TODO: replace with the real variance -- the below is a 
    % calculation under the null hypothesis that X indep. Y!
    sigma = sqrt( (2*(2*numPts+5))./(9*numPts.*(numPts-1)) );
end