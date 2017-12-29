function [outCell,x1OvlpIdxs,x2OvlpIdxs,x1UniqueIdxs,x2UniqueIdxs] = separateOverlaps_binary(x,y)
% y - the output variable, should be binary

outCell = cell(1,2);
% first  index contains [x,y] which isn't overlapping
% second index contains [x,y] which is overlapping

uniquePredictorVals = unique(y);

idxSet1 = (y==uniquePredictorVals(1));
idxSet2 = (y==uniquePredictorVals(2));
x1 = x(idxSet1); y1 = y(idxSet1);
x2 = x(idxSet2); y2 = y(idxSet2);

x1OvlpIdxs = (x1>=min(x2)); x1UniqueIdxs = (x1<min(x2));
x2OvlpIdxs = (x2<=max(x1)); x2UniqueIdxs = (x2>max(x1));

XY1 = [[x1(x1OvlpIdxs) y1(x1OvlpIdxs)]; ...
       [x2(x2OvlpIdxs) y2(x2OvlpIdxs)]];
XY2 = [[x1(x1UniqueIdxs) y1(x1UniqueIdxs)]; ...
       [x2(x2UniqueIdxs) y2(x2UniqueIdxs)]];

outCell{1} = XY2;  % put the unique stuff (non-overlapping) first
outCell{2} = XY1;  % put the overlapping points second

end