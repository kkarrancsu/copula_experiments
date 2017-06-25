%% Test Tau prediction with region shift

clear;
clc;

M = 200;
x = rand(M,1);
y = 4*(x-0.5).^2;

[x,I] = sort(x);
y = y(I);

numPtsIncr = 10;
MExploreVec = M/2:numPtsIncr:M;

ii = 0;
for MM=MExploreVec
    numR2 = MM-M/2;
    x_subset = x(1:MM); y_subset = y(1:MM);
    tauActual = corr(x_subset,y_subset,'type','kendall');
    
    xx = fliplr(x(M/2:M/2+numR2)); yy = y(M/2:M/2+numR2);
    scatter(x_subset,y_subset,'b'); hold on;
    xx_subset = x(1:M/2); yy_subset = y(1:M/2);
    u = [xx_subset; xx]; v = [yy_subset; yy];
    scatter(xx,yy,'r');
    pause;
    tauPredict = corr(u,v,'type','kendall');
    
    fprintf('tauActual=%0.02f tauPredict=%0.02f\n', ...
        tauActual, tauPredict);
end