function [] = nonmonotonic_search(data, alpha, monotonic_thresh)

[R,RectanglesCell] = paircim_v4(data);
sz = size(R);
M = size(data,1);  % used for p-value computations

resultsCell = {}; resultsCellIdx = 1;
for ii=1:sz(1)
    for jj=ii+1:sz(2)
        % check if dependency is "significant"
        cimVal = R(ii,jj);
        pval = cimv4_pval(cimVal,M);
        if(pval<alpha)
            % see if we are within the monotonic_thresh to prevent
            % over-fitting
            tauklVal = abs(taukl(data(:,ii),data(:,jj)));
            percentageDiff = abs(cimVal-tauklVal)/tauklVal;
            if(percentageDiff>monotonic_thresh)
                % means we didn't overfit and we now store the regions and
                % the dependency indices
                rco = RectanglesCell{ii,jj};
                % count the number of regions
                x = struct;
                x.feature1 = ii;
                x.feature2 = jj;
                x.depVal = cimVal;
                x.numRegions = size(rco,2);
                resultsCell{resultsCellIdx} = x;
                resultsCellIdx = resultsCellIdx + 1;
            end
        end
    end
end

end