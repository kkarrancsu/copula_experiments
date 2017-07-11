%% debug cim_v8a_rev3cc

clear;
clc;
close all;

rng(1);

M = 500;
minscanincr = 0.015625;
numMCSim = 1000;

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');
for ii=1:numMCSim
    dispstat(sprintf('%d/%d',ii, numMCSim),'timestamp');
    
    x = rand(M,1);
    y = (x > 0.5);
    m = cim_v8a_rev3cc(x,y,minscanincr);
end