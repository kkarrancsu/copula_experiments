% domain experiments for the cos
M = 500;
x = rand(M,1)*10-15;
y = rand(M,1)*10-15;
z = y.*sin(x) - x.*cos(y);
% scatter3(x,y,z)

regions = cos_regions_3d(x,y,z)

%% How to create a lexicographical ordering?
clear;
clc;
close all;
dbstop if error;

M = 500;
x = rand(M,1)*8-4;
% y = rand(M,1)*8-4;
% z = -x.^2-y.^2+6;
y = x;
% z = x.^2 + y.^2;
z = sin(x) + cos(y);

[cosValue, RR] = cos3d(x,y,z);

% % do some initial verification of the RR matrix
numDomains = length(unique(RR(:,6)))
for domainNum=1:numDomains
    % get all data for RR in this domain
    II = RR(:,6)==domainNum;
    fprintf('************************************************************\n');
    SS = RR(II,:);
    eCopVec = SS(2:end,5);
    [eCopVecSorted,eCopSortII] = sort(eCopVec);
    distVec = pdist2(SS(1,2:4),SS(2:end,2:4));
    distVecSorted = distVec(eCopSortII');
    [eCopVecSorted distVecSorted']
    issorted(distVecSorted)
    fprintf('************************************************************\n');
end