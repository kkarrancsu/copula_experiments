function [regionIdxs] = findruns(zScoreProfile)

signDiffVec = sign([diff(zScoreProfile)]);
regionIdxs = [];
for ii=2:length(signDiffVec)
    if(signDiffVec(ii)~=signDiffVec(ii-1))
        regionIdxs = [regionIdxs ii-1];
    end
end

end
