function [avgMetric, avgNumPts] = getSignatureAvg(signature, mIdx, lIdx, numMCSims, scanincrs)
    for scanincrIdx=1:length(scanincrs)
        scanincr = scanincrs(scanincrIdx);
        zz = signature{mIdx,lIdx,1}(scanincr);
        metricVec = zz.metricVec;
        numPtsVec = zz.numPtsVec;
        for ii=2:numMCSims
            zz = signature{mIdx,lIdx,ii}(scanincr);
            metricVec = metricVec + zz.metricVec;
            numPtsVec = numPtsVec + zz.numPtsVec;
        end
        avgMetricLoc = metricVec/numMCSims;
        avgNumPtsLoc = numPtsVec/numMCSims;
        
        avgMetric{scanincrIdx} = avgMetricLoc;
        avgNumPts{scanincrIdx} = avgNumPtsLoc;
    end
end