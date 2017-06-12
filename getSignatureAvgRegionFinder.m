function [avgResults] = getSignatureAvgRegionFinder(signature, mIdx, lIdx, numMCSims, scanincrs)
    avgResults = cell(1,length(scanincrs));
    for scanincrIdx=1:length(scanincrs)
        scanincr = scanincrs(scanincrIdx);
        zz = signature{mIdx,lIdx,1}(scanincr);
        for ii=2:numMCSims
            zz = zz + signature{mIdx,lIdx,ii}(scanincr);
        end
        zz = zz./numMCSims;
        avgResults{scanincrIdx} = zz;
    end
end