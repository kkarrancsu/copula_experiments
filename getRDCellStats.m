function [numRegionsMat, regionValuesCell] = getRDCellStats(rdCell, MVec, noiseVec, cornerPts, scanincrs, ...
                                                           MvalOfInterest, noiseValsOfInterest, cornerPtOfInterest,  ...
                                                           scanincrOfInterest, numSimsToCompute)

MIdx = find(MvalOfInterest==MVec);
noiseIdxs = zeros(1,length(noiseValsOfInterest));
for ii=1:length(noiseIdxs)
    noiseIdxs(ii) = find(noiseValsOfInterest(ii)==noiseVec);
end
cornerPtIdx = find(cornerPtOfInterest==cornerPts);
scanincrIdx = find(scanincrOfInterest==scanincrs);

numRegionsMat = zeros(numSimsToCompute,length(noiseValsOfInterest));
regionValuesCell = cell(1,length(noiseValsOfInterest));

for ii=1:length(noiseIdxs)
    noiseIdx = noiseIdxs(ii);
    rawDataAccum = [];
    for jj=1:numSimsToCompute
        rawData = rdCell{MIdx, noiseIdx, cornerPtIdx, scanincrIdx, jj};
        numRegionsMat(jj,ii) = length(rawData);
        rawDataAccum = [rawDataAccum rawData];
    end
    regionValuesCell{ii} = rawDataAccum;
end

end