%% Script which characterizes the CIM algorithm in many different ways

%% Characterize the CIM Algorithm Power
clear;
clc;
close all;
dbstop if error;

if(ispc)
    load('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\power_all_cimv4.mat');
elseif(ismac)
    load('/Users/Kiran/ownCloud/PhD/sim_results/independence/power_all_cimv4.mat');
else
    load('/home/kiran/ownCloud/PhD/sim_results/independence/power_all_cimv4.mat');
end

% select for which M to plot
M = 500;
% find the index we need to extract
MIdx = find(M_vec==M);

num_noise_test_min = 1;
num_noise_test_max = 20;
noiseVec = num_noise_test_min:num_noise_test_max;
powerMat = zeros(6,8,length(noiseVec));
powerMat(1,:,:) = rsdmPower(:,noiseVec,MIdx);
powerMat(2,:,:) = cosPower(:,noiseVec,MIdx);
powerMat(3,:,:) = rdcPower(:,noiseVec,MIdx);
powerMat(4,:,:) = ticePower(:,noiseVec,MIdx);
powerMat(5,:,:) = dcorrPower(:,noiseVec,MIdx);
powerMat(6,:,:) = ccorrPower(:,noiseVec,MIdx);
noiseVec = (num_noise_test_min:num_noise_test_max)/10;

labels = {'CIM', 'CoS', 'RDC', 'TICe', 'dCor', 'cCor'};
plotStyle = 1;
plotPower(powerMat, M, labels, noiseVec, num_noise_test_min, num_noise_test_max, plotStyle)

%% test the power sensitivity of the CIM algorithm

clear;
clc;
close all;

cimVersion = 6;
scanincrsToTest = [0.25, 0.125, .0625, .03125, .015625];
switch cimVersion
    case 1
        cimfunc = @cim;
        fnameStr = 'cim';
    case 3
        cimfunc = @cim_v3;
        fnameStr = 'cimv3';
    case 4
        cimfunc = @cim_v4;
        fnameStr = 'cimv4';
    case 5
        cimfunc = @cim_v5;
        fnameStr = 'cimv5';
    case 6
        cimfunc = @cim_v6;
        fnameStr = 'cimv6';
    case 7
        cimfunc = @cim_v7;
        fnameStr = 'cimv7';
    otherwise
end
MVecToTest = 100:100:1000;
for M=MVecToTest
    powerCurve = cim_power_sensitivity(cimfunc,M,scanincrsToTest);
    % save the results
    if(ispc)
        save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\%s_powerSensitivity_M_%d.mat', fnameStr, M));
    elseif(ismac)
        save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/%s_powerSensitivity_M_%d.mat', fnameStr, M));
    else
        save(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/%s_powerSensitivity_M_%d.mat', fnameStr, M));
    end
end

%% plot the power sensitivity
clear;
clc;
close all;

M = 900;
cimVersion = 4;

% load the correct data
switch cimVersion
    case 1
        fnameStr = 'cim';
    case 3
        fnameStr = 'cimv3';
    case 4
        fnameStr = 'cimv4';
    case 5
        fnameStr = 'cimv5';
    case 6
        fnameStr = 'cimv6';
    case 7
        fnameStr = 'cimv7';
    otherwise
end
if(ispc)
    load(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\%s_powerSensitivity_M_%d.mat', fnameStr, M));
elseif(ismac)
    load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/%s_powerSensitivity_M_%d.mat', fnameStr, M));
else
    load(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/%s_powerSensitivity_M_%d.mat', fnameStr, M));
end

% TODO: use the plotPower function to plot the power for different scanning
% intervals to show sensitivity
labels = {'0.25', '0.125', '0.0625', '0.03125', '0.015625'};
plotStyle = 1;
num_noise_test_min = 1;
num_noise_test_max = 20;
noiseVec = num_noise_test_min:num_noise_test_max;
plotPower(powerCurve, M, labels, noiseVec, num_noise_test_min, num_noise_test_max, plotStyle)

%% Plot the power-sensitivity in a different way
clear;
clc;
close all;
dbstop if error;

cimVersion = 4;
MVecToPlot = 100:100:1000;
num_noise_test_min = 1;
num_noise_test_max = 20+1;
plotPowerSensitivity_withinM(cimVersion,MVecToPlot,num_noise_test_min,num_noise_test_max);
figtitle(sprintf('Algorithm Power Sensitivity (M=%d - %d)',min(MVecToPlot),max(MVecToPlot)));

%% Run the algorithm sensitivity analysis

clear;
clc;
close all;

cimVersion = 4;
scanincrsToTest = [0.25, 0.125, .0625, .03125, .015625];
switch cimVersion
    case 1
        cimfunc = @cim;
        fnameStr = 'cim';
    case 3
        cimfunc = @cim_v3;
        fnameStr = 'cimv3';
    case 4
        cimfunc = @cim_v4;
        fnameStr = 'cimv4';
    case 5
        cimfunc = @cim_v5;
        fnameStr = 'cimv5';
    case 6
        cimfunc = @cim_v6;
        fnameStr = 'cimv6';
    case 7
        cimfunc = @cim_v7;
        fnameStr = 'cimv7';
    otherwise
end
MVecToTest = 100:100:1000;
for M=MVecToTest
    algoSensitivityData = cim_algo_sensitivity(cimfunc,M,scanincrsToTest);    
    % save the results
    if(ispc)
        save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\%s_algoSensitivity_M_%d.mat', fnameStr, M));
    elseif(ismac)
        save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/%s_algoSensitivity_M_%d.mat', fnameStr, M));
    else
        save(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/%s_algoSensitivity_M_%d.mat', fnameStr, M));
    end
end

%% Plot the algorithm sensitivity analysis
clear;
clc;
close all;
dbstop if error;

M = 1000;
cimVersion = 4;

% load the correct data
switch cimVersion
    case 1
        fnameStr = 'cim';
    case 3
        fnameStr = 'cimv3';
    case 4
        fnameStr = 'cimv4';
    case 5
        fnameStr = 'cimv5';
    case 6
        fnameStr = 'cimv6';
    case 7
        fnameStr = 'cimv7';
    otherwise
end
if(ispc)
    load(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\%s_algoSensitivity_M_%d.mat', fnameStr, M));
elseif(ismac)
    load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/%s_algoSensitivity_M_%d.mat', fnameStr, M));
else
    load(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/%s_algoSensitivity_M_%d.mat', fnameStr, M));
end

num_noise_test_min = 1;
num_noise_test_max = 20;
noiseVec = num_noise_test_min:num_noise_test_max;

plotAlgoSensitivity(algoSensitivityData, scanincrsToTest, noiseVec, M)

%% Plot the Algorithm sensitivity analysis including M
clear;
clc;
close all;
dbstop if error;

cimVersion = 4;
MVecToPlot = 100:100:1000;
% plotAlgoSensitivity_acrossM(cimVersion,MVecToPlot);
num_noise_test_min = 1;
num_noise_test_max = 20+1;

plotAlgoSensitivity_withinM(cimVersion,MVecToPlot,num_noise_test_min,num_noise_test_max);
figtitle(sprintf('Algorithm Sensitivity (M=%d - %d)',min(MVecToPlot),max(MVecToPlot)));
%% simulate the difference between computing monotonic, detecting the region, and theoretical CIM value
clear;
clc;

xMin = 0;
xMax = 1;
num_noise = 30;                    % The number of different noise levels used
noise = 3;                         % A constant to determine the amount of noise

MVecResults = [100:100:1000 2500 5000 10000];
num_noise_test_min = 0;
num_noise_test_max = 30;
noiseVec = num_noise_test_min:num_noise_test_max;
numMCSim = 200;

minscanincrVal = 0.015625;

cimVersion = 4;

switch(cimVersion)
    case 1
        fnameStr = 'cimv1';
        cimfunc = @cim;
    case 3
        fnameStr = 'cimv3';
        cimfunc = @cim_v3;
    case 4
        fnameStr = 'cimv4';
        cimfunc = @cim_v4;
    case 5
        fnameStr = 'cimv5';
        cimfunc = @cim_v5;
    case 6
        fnameStr = 'cimv6';
        cimfunc = @cim_v6;
    case 7
        fnameStr = 'cimv7';
        cimfunc = @cim_v7;
    otherwise
        error('Unrecognized CIM version!');
end
        
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');

for M=MVecResults
    % Test Empirical copula distance as a distance measure under various
    % amounts of noise
    linearDep     = zeros(3,length(noiseVec));
    quadraticDep  = zeros(3,length(noiseVec));
    cubicDep      = zeros(3,length(noiseVec));
    sinusoidalDep = zeros(3,length(noiseVec));
    hiFreqSinDep  = zeros(3,length(noiseVec));
    fourthRootDep = zeros(3,length(noiseVec));
    circleDep     = zeros(3,length(noiseVec));
    stepDep       = zeros(3,length(noiseVec));
    indep         = zeros(3,length(noiseVec));
                    % 1 - using tau directly
                    % 2 - segmenting into regions
                  
    for l=noiseVec
        dispstat(sprintf('Computing for noise level=%d >> M=%d',l,M),'keepthis', 'timestamp');
        t1 = 0; c1 = 0; c11 = 0;
        t2 = 0; c2 = 0; c22 = 0;
        t3 = 0; c3 = 0; c33 = 0;
        t4 = 0; c4 = 0; c44 = 0;
        t6 = 0; c6 = 0; c66 = 0;
        t5 = 0; c5 = 0; c55 = 0;
        t7 = 0; c7 = 0; c77 = 0;
        t8 = 0; c8 = 0; c88 = 0;
        t9 = 0; c9 = 0; c99 = 0;
        parfor mcSimNum=1:numMCSim
            x = rand(M,1)*(xMax-xMin)+xMin;
            y1 = x + noise*(l/num_noise)*randn(M,1);
            y2 = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
            y3 = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
            y4 = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
            y5 = sin(16*pi*x) + noise*(l/num_noise)*randn(M,1);
            y6 = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
            y7=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
            y8 = (x > 0.5) + noise*5*l/num_noise*randn(M,1);
            y9 = rand(M,1)*(xMax-xMin)+xMin;

            U1 = pobs([x y1]);
            U2 = pobs([x y2]);
            U3 = pobs([x y3]);
            U4 = pobs([x y4]);
            U5 = pobs([x y5]);
            U6 = pobs([x y6]);
            U7 = pobs([x y7]);
            U8 = pobs([x y8]);
            U9 = pobs([x y9]);
            
            t1 = t1 + abs(corr(x,y1,'type','kendall'));
            c1 = c1 + abs(corr(x,y1,'type','kendall'));
            c11 = c11 + cimfunc(x,y1,minscanincrVal);
            
            t2 = t2 + abs(corr(x,y2,'type','kendall'));
            r1 = inBoundedPts(U2(:,1),U2(:,2),0,0.5,0,1); 
            r2 = inBoundedPts(U2(:,1),U2(:,2),0.5,1,0,1);
            c2 = c2 + ( abs(corr(r1(:,1),r1(:,2),'type','kendall'))*0.5 + ...
                        abs(corr(r2(:,1),r2(:,2),'type','kendall'))*0.5 );
            c22 = c22 + cimfunc(x,y2,minscanincrVal);
            
            t3 = t3 + abs(corr(x,y3,'type','kendall'));
            r1 = inBoundedPts(U3(:,1),U3(:,2),0,0.1138,0,1); 
            r2 = inBoundedPts(U3(:,1),U3(:,2),0.1138,0.5768,0,1);
            r3 = inBoundedPts(U3(:,1),U3(:,2),0.5768,1,0,1);
            c3 = c3 + ( abs(corr(r1(:,1),r1(:,2),'type','kendall'))*0.1138 + ...
                        abs(corr(r2(:,1),r2(:,2),'type','kendall'))*0.4368 + ...
                        abs(corr(r3(:,1),r3(:,2),'type','kendall'))*0.4232);
            c33 = c33 + cimfunc(x,y3,minscanincrVal);
                    
            t4 = t4 + abs(corr(x,y4,'type','kendall'));
            r1 = inBoundedPts(U4(:,1),U4(:,2),0,0.1138,0,1); 
            r2 = inBoundedPts(U4(:,1),U4(:,2),0.1138,0.3413,0,1);
            r3 = inBoundedPts(U4(:,1),U4(:,2),0.3413,.6208,0,1);
            r4 = inBoundedPts(U4(:,1),U4(:,2),0.6208,.8483,0,1);
            r5 = inBoundedPts(U4(:,1),U4(:,2),0.8483,1,0,1);
            c4 = c4 + ( abs(corr(r1(:,1),r1(:,2),'type','kendall'))*0.1138 + ...
                        abs(corr(r2(:,1),r2(:,2),'type','kendall'))*0.2275 + ...
                        abs(corr(r3(:,1),r3(:,2),'type','kendall'))*0.2795 + ...
                        abs(corr(r4(:,1),r4(:,2),'type','kendall'))*0.2275 + ...
                        abs(corr(r5(:,1),r5(:,2),'type','kendall'))*0.1517);
            c44 = c44 + cimfunc(x,y4,minscanincrVal);
            
            t5 = t5 + abs(corr(x,y5,'type','kendall'));
            zz1 = 0; zz2 = 0.02900;   r1 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w1 = zz2-zz1;
            zz1 = zz2; zz2 = 0.09158; r2 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w2 = zz2-zz1;
            zz1 = zz2; zz2 = 0.1528;  r3 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w3 = zz2-zz1;
            zz1 = zz2; zz2 = 0.2108;  r4 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w4 = zz2-zz1;
            zz1 = zz2; zz2 = 0.2715;  r5 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w5 = zz2-zz1;
            zz1 = zz2; zz2 = 0.3339;  r6 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w6 = zz2-zz1;
            zz1 = zz2; zz2 = 0.3999;  r7 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w7 = zz2-zz1;
            zz1 = zz2; zz2 = 0.4613;  r8 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w8 = zz2-zz1;
            zz1 = zz2; zz2 = 0.5227;  r9 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w9 = zz2-zz1;
            zz1 = zz2; zz2 = 0.5893;  r10 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w10 = zz2-zz1;
            zz1 = zz2; zz2 = 0.6547;  r11 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w11 = zz2-zz1;
            zz1 = zz2; zz2 = 0.7135;  r12 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w12 = zz2-zz1;
            zz1 = zz2; zz2 = 0.7734;  r13 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w13 = zz2-zz1;
            zz1 = zz2; zz2 = 0.8438;  r14 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w14 = zz2-zz1;
            zz1 = zz2; zz2 = 0.9074;  r15 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w15 = zz2-zz1;
            zz1 = zz2; zz2 = 0.9694;  r16 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w16 = zz2-zz1;
            zz1 = zz2; zz2 = 1;       r17 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w17 = zz2-zz1;
            c5 = c5 + ( abs(corr(r1(:,1),r1(:,2),'type','kendall'))*w1 + ...
                        abs(corr(r2(:,1),r2(:,2),'type','kendall'))*w2 + ...
                        abs(corr(r3(:,1),r3(:,2),'type','kendall'))*w3 + ...
                        abs(corr(r4(:,1),r4(:,2),'type','kendall'))*w4 + ...
                        abs(corr(r5(:,1),r5(:,2),'type','kendall'))*w5 + ...
                        abs(corr(r6(:,1),r6(:,2),'type','kendall'))*w6 + ...
                        abs(corr(r7(:,1),r7(:,2),'type','kendall'))*w7 + ...
                        abs(corr(r8(:,1),r8(:,2),'type','kendall'))*w8 + ...
                        abs(corr(r9(:,1),r9(:,2),'type','kendall'))*w9 + ...
                        abs(corr(r10(:,1),r10(:,2),'type','kendall'))*w10 + ...
                        abs(corr(r11(:,1),r11(:,2),'type','kendall'))*w11 + ...
                        abs(corr(r12(:,1),r12(:,2),'type','kendall'))*w12 + ...
                        abs(corr(r13(:,1),r13(:,2),'type','kendall'))*w13 + ...
                        abs(corr(r14(:,1),r14(:,2),'type','kendall'))*w14 + ...
                        abs(corr(r15(:,1),r15(:,2),'type','kendall'))*w15 + ...
                        abs(corr(r16(:,1),r16(:,2),'type','kendall'))*w16 + ...
                        abs(corr(r17(:,1),r17(:,2),'type','kendall'))*w17);
            c55 = c55 + cimfunc(x,y5,minscanincrVal);
            
            t6 = t6 + abs(corr(x,y6,'type','kendall'));
            c6 = c6 + abs(corr(x,y6,'type','kendall'));
            c66 = c66 + cimfunc(x,y6,minscanincrVal);
            
            t7 = t7 + abs(corr(x,y7,'type','kendall'));
            r1 = inBoundedPts(U7(:,1),U7(:,2),0,0.5,0,0.5); 
            r2 = inBoundedPts(U7(:,1),U7(:,2),0.5,1,0,0.5);
            r3 = inBoundedPts(U7(:,1),U7(:,2),0,0.5,0.5,1);
            r4 = inBoundedPts(U7(:,1),U7(:,2),0.5,1,0.5,1);
            c7 = c7 + ( abs(corr(r1(:,1),r1(:,2),'type','kendall'))*0.25 + ...
                        abs(corr(r2(:,1),r2(:,2),'type','kendall'))*0.25 + ...
                        abs(corr(r3(:,1),r3(:,2),'type','kendall'))*0.25 + ...
                        abs(corr(r4(:,1),r4(:,2),'type','kendall'))*0.25);
            c77 = c77 + cimfunc(x,y7,minscanincrVal);
            
            t8 = t8 + abs(corr(x,y8,'type','kendall'));
            c8 = c8 + abs(corr(x,y8,'type','kendall'));
            c88 = c88 + cimfunc(x,y8,minscanincrVal);
            
            t9 = t9 + abs(corr(x,y9,'type','kendall'));
            c9 = c9 + abs(corr(x,y9,'type','kendall'));
            c99 = c99 + cimfunc(x,y9,minscanincrVal);
        end
        linearDep(1,l+1) = t1/numMCSim; linearDep(2,l+1) = c1/numMCSim; linearDep(3,l+1) = c11/numMCSim;
        quadraticDep(1,l+1) = t2/numMCSim; quadraticDep(2,l+1) = c2/numMCSim; quadraticDep(3,l+1) = c22/numMCSim;
        cubicDep(1,l+1) = t3/numMCSim; cubicDep(2,l+1) = c3/numMCSim; cubicDep(3,l+1) = c33/numMCSim;
        sinusoidalDep(1,l+1) = t4/numMCSim; sinusoidalDep(2,l+1) = c4/numMCSim; sinusoidalDep(3,l+1) = c44/numMCSim;
        hiFreqSinDep(1,l+1) = t5/numMCSim; hiFreqSinDep(2,l+1) = c5/numMCSim; hiFreqSinDep(3,l+1) = c55/numMCSim;
        fourthRootDep(1,l+1) = t6/numMCSim; fourthRootDep(2,l+1) = c6/numMCSim; fourthRootDep(3,l+1) = c66/numMCSim;
        circleDep(1,l+1) = t7/numMCSim; circleDep(2,l+1) = c7/numMCSim; circleDep(3,l+1) = c77/numMCSim;
        stepDep(1,l+1) = t8/numMCSim; stepDep(2,l+1) = c8/numMCSim; stepDep(3,l+1) = c88/numMCSim;
        indep(1,l+1) = t9/numMCSim; indep(2,l+1) = c9/numMCSim; indep(3,l+1) = c99/numMCSim;
    end

    % save the data
    if(ispc)
        save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\regionDetection_%s_M_%d.mat', fnameStr, M));
    elseif(ismac)
        save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/clustering/regionDetection_%s_M_%d.mat', fnameStr, M));
    else
        save(sprintf('/home/kiran/ownCloud/PhD/sim_results/clustering/regionDetection_%s_M_%d.mat', fnameStr, M));
    end
end
%% Plot the results
clear;
clc;

cim_version = 4;
switch cim_version
    case 3
        fnameStr = 'cimv3';
    case 4
        fnameStr = 'cimv4';
    case 5
        fnameStr = 'cimv5';
    case 6
        fnameStr = 'cimv6';
    case 7
        fnameStr = 'cimv7';
    otherwise
        error('no data exists!');
end

MVecDataAvailable = [100:100:1000 2500];  % add 5000 and 10000 to this as they come online
tol = 0.015;
noiseMinToAnalyze = 0; noiseMaxToAnalyze = 20;
noiseVecToAnalyze = noiseMinToAnalyze:noiseMaxToAnalyze;  

% do a +1 to account for matlab indexing
noiseVecToAnalyze = noiseVecToAnalyze + 1;

% aggregatefunc = @sum;
% aggregatefunc = @max;
aggregatefunc = @mean;

numDep = 8;
depTestVec = ones(1,numDep);
MVecResults = zeros(1,numDep);
for MIdx=1:length(MVecDataAvailable)
    M = MVecDataAvailable(MIdx);
    if(ispc)
        load(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\regionDetection_%s_M_%d.mat', fnameStr, M));
    elseif(ismac)
        load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/clustering/regionDetection_%s_M_%d.mat', fnameStr, M));
    else
        load(sprintf('/home/kiran/ownCloud/PhD/sim_results/clustering/regionDetection_%s_M_%d.mat', fnameStr, M));
    end
    
    % compute error between theoretical CIM and actual CIM for each of the
    % dependencies, if the total error is within the sseTol value, then we
    % store that as the M Value for which the error is acceptable.
    
    if(depTestVec(1))
%         linearDepErr = (linearDep(2,noiseToTest)-linearDep(3,noiseToTest)).^2;
        linearDepErr = abs(linearDep(2,noiseVecToAnalyze)-linearDep(3,noiseVecToAnalyze));
        if(aggregatefunc(linearDepErr)<=tol)
            MVecResults(1) = M;
            depTestVec(1) = 0;
            linearDepToPlot = linearDep;
        end
    end
    
    if(depTestVec(2))
%         quadraticDepErr = (quadraticDep(2,noiseToTest)-quadraticDep(3,noiseToTest)).^2;
        quadraticDepErr = abs(quadraticDep(2,noiseVecToAnalyze)-quadraticDep(3,noiseVecToAnalyze));
        if(aggregatefunc(quadraticDepErr)<=tol)
            MVecResults(2) = M;
            depTestVec(2) = 0;
            quadraticDepToPlot = quadraticDep;
        end
    end
    
    if(depTestVec(3))
%         cubicDepErr = (cubicDep(2,noiseToTest)-cubicDep(3,noiseToTest)).^2;
        cubicDepErr = abs(cubicDep(2,noiseVecToAnalyze)-cubicDep(3,noiseVecToAnalyze));
        if(aggregatefunc(cubicDepErr)<=tol)
            MVecResults(3) = M;
            depTestVec(3) = 0;
            cubicDepToPlot = cubicDep;
        end
    end
    
    if(depTestVec(4))
%         sinusoidalDepErr = (sinusoidalDep(2,noiseToTest)-sinusoidalDep(3,noiseToTest)).^2;
        sinusoidalDepErr = abs(sinusoidalDep(2,noiseVecToAnalyze)-sinusoidalDep(3,noiseVecToAnalyze));
        if(aggregatefunc(sinusoidalDepErr)<=tol)
            MVecResults(4) = M;
            depTestVec(4) = 0;
            sinusoidalDepToPlot = sinusoidalDep;
        end
    end
    
    if(depTestVec(5))
%         hiFreqSinDepErr = (hiFreqSinDep(2,noiseToTest)-hiFreqSinDep(3,noiseToTest)).^2;
%         hiFreqSinDep(2,1) = 1;
        hiFreqSinDepErr = abs(hiFreqSinDep(2,noiseVecToAnalyze)-hiFreqSinDep(3,noiseVecToAnalyze));
        if(aggregatefunc(hiFreqSinDepErr)<=tol)
            MVecResults(5) = M;
            depTestVec(5) = 0;
            hiFreqSinDepToPlot = hiFreqSinDep;
        end
    end
    
    if(depTestVec(6))
%         fourthRootDepErr = (fourthRootDep(2,noiseToTest)-fourthRootDep(3,noiseToTest)).^2;
        fourthRootDepErr = abs(fourthRootDep(2,noiseVecToAnalyze)-fourthRootDep(3,noiseVecToAnalyze));
        if(aggregatefunc(fourthRootDepErr)<=tol)
            MVecResults(6) = M;
            depTestVec(6) = 0;
            fourthRootDepToPlot = fourthRootDep;
        end
    end
    
    if(depTestVec(7))
%         circleDepErr = (circleDep(2,noiseToTest)-circleDep(3,noiseToTest)).^2;
        circleDepErr = abs(circleDep(2,noiseVecToAnalyze)-circleDep(3,noiseVecToAnalyze));
        if(aggregatefunc(circleDepErr)<=tol)
            MVecResults(7) = M;
            depTestVec(7) = 0;
            circleDepToPlot = circleDep;
        end
    end
    
    if(depTestVec(8))
%         stepDepErr = (stepDep(2,noiseToTest)-stepDep(3,noiseToTest)).^2;
        stepDep(2,1) = 1; % due to a bug in the way we simulated we need to force
                          % this value to 1
        stepDepErr = abs(stepDep(2,noiseVecToAnalyze)-stepDep(3,noiseVecToAnalyze));
        if(aggregatefunc(stepDepErr)<=tol)
            MVecResults(8) = M;
            depTestVec(8) = 0;
            stepDepToPlot = stepDep;
        end
    end
end

% use max M for any remainders that didn't meet the SSE
for depTestValIdx=1:length(depTestVec)
    if(depTestVec(depTestValIdx)==1)
        M = max(MVecDataAvailable);
        if(ispc)
            load(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\regionDetection_%s_M_%d.mat', fnameStr, M));
        elseif(ismac)
            load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/clustering/regionDetection_%s_M_%d.mat', fnameStr, M));
        else
            load(sprintf('/home/kiran/ownCloud/PhD/sim_results/clustering/regionDetection_%s_M_%d.mat', fnameStr, M));
        end
        MVecResults(depTestValIdx) = M;
        switch depTestValIdx
            case 1
                linearDepToPlot = linearDep;
            case 2
                quadraticDepToPlot = quadraticDep;
            case 3
                cubicDepToPlot = cubicDep;
            case 4
                sinusoidalDepToPlot = sinusoidalDep;
            case 5
                hiFreqSinDepToPlot = hiFreqSinDep;
%                 hiFreqSinDepToPlot(2,1) = 1;
            case 6
                fourthRootDepToPlot = fourthRootDep;
            case 7
                circleDepToPlot = circleDep;
            case 8
                stepDepToPlot = stepDep;
                stepDepToPlot(2,1) = 1;
            otherwise
                error('UNK!');
        end
    end
end

noiseVecPlot = noiseVecToAnalyze-1;

subplot(2,4,1);
% hh1 = plot(noiseToTest,linearDepToPlot(1,noiseToTest),'o-.', ...
%      noiseToTest,linearDepToPlot(2,noiseToTest),'+-.', ...
%      noiseToTest,linearDepToPlot(3,noiseToTest),'d-.');
hh1 = plot(noiseVecPlot/10,linearDepToPlot(2,noiseVecToAnalyze),'+-.', ...
     noiseVecPlot/10,linearDepToPlot(3,noiseVecToAnalyze),'d-.');
grid on;
xlabel('noise');
% legend('\tau','CIM Theoretical',legendStr);
h = legend('CIM','$$\widehat{CIM}$$');
set(h,'Interpreter','latex')
title({'Linear Dependency', sprintf('min(M)=%d',MVecResults(1))});
hh1(1).LineWidth = 1.5; 
hh1(2).LineWidth = 1.5; 
% hh1(3).LineWidth = 1.5; 

subplot(2,4,2);
hh1 = plot(noiseVecPlot/10,quadraticDepToPlot(2,noiseVecToAnalyze),'+-.', ...
     noiseVecPlot/10,quadraticDepToPlot(3,noiseVecToAnalyze),'d-.');
grid on;
xlabel('noise');
title({'Quadratic Dependency', sprintf('min(M)=%d',MVecResults(2))});
hh1(1).LineWidth = 1.5; 
hh1(2).LineWidth = 1.5; 
% hh1(3).LineWidth = 1.5; 

subplot(2,4,3);
hh1 = plot(noiseVecPlot/10,cubicDepToPlot(2,noiseVecToAnalyze),'+-.', ...
     noiseVecPlot/10,cubicDepToPlot(3,noiseVecToAnalyze),'d-.');
grid on;
xlabel('noise');
title({'Cubic Dependency', sprintf('min(M)=%d',MVecResults(3))});
hh1(1).LineWidth = 1.5; 
hh1(2).LineWidth = 1.5; 
% hh1(3).LineWidth = 1.5; 

subplot(2,4,4);
hh1 = plot(noiseVecPlot/10,sinusoidalDepToPlot(2,noiseVecToAnalyze),'+-.', ...
     noiseVecPlot/10,sinusoidalDepToPlot(3,noiseVecToAnalyze),'d-.');
grid on;
xlabel('noise');
title({'LF-Sin Dependency', sprintf('min(M)=%d',MVecResults(4))});
hh1(1).LineWidth = 1.5; 
hh1(2).LineWidth = 1.5; 
% hh1(3).LineWidth = 1.5; 

subplot(2,4,5);
hh1 = plot(noiseVecPlot/10,hiFreqSinDepToPlot(2,noiseVecToAnalyze),'+-.', ...
     noiseVecPlot/10,hiFreqSinDepToPlot(3,noiseVecToAnalyze),'d-.');
grid on;
xlabel('noise');
title({'HF-Sin Dependency', sprintf('min(M)=%d',MVecResults(5))});
hh1(1).LineWidth = 1.5; 
hh1(2).LineWidth = 1.5; 
% hh1(3).LineWidth = 1.5; 

subplot(2,4,6);
hh1 = plot(noiseVecPlot/10,fourthRootDepToPlot(2,noiseVecToAnalyze),'+-.', ...
     noiseVecPlot/10,fourthRootDepToPlot(3,noiseVecToAnalyze),'d-.');
grid on;
xlabel('noise');
title({'Fourth-Root Dependency', sprintf('min(M)=%d',MVecResults(6))});
hh1(1).LineWidth = 1.5; 
hh1(2).LineWidth = 1.5; 
% hh1(3).LineWidth = 1.5; 

subplot(2,4,7);
hh1 = plot(noiseVecPlot/10,circleDepToPlot(2,noiseVecToAnalyze),'+-.', ...
     noiseVecPlot/10,circleDepToPlot(3,noiseVecToAnalyze),'d-.');
grid on;
xlabel('noise');
title({'Circular Dependency', sprintf('min(M)=%d',MVecResults(7))});
hh1(1).LineWidth = 1.5; 
hh1(2).LineWidth = 1.5; 
% hh1(3).LineWidth = 1.5; 

subplot(2,4,8);
hh1 = plot(noiseVecPlot/10,stepDepToPlot(2,noiseVecToAnalyze),'+-.', ...
     noiseVecPlot/10,stepDepToPlot(3,noiseVecToAnalyze),'d-.');
grid on;
xlabel('noise');
title({'Step-Function Dependency', sprintf('min(M)=%d',MVecResults(8))});
hh1(1).LineWidth = 1.5; 
hh1(2).LineWidth = 1.5; 
% hh1(3).LineWidth = 1.5; 

% subplot(3,3,9);
% hh1 = plot(noiseVec,indep(1,:),'o-.', ...
%      noiseVec,indep(2,:),'+-.', ...
%      noiseVec,indep(3,:),'d-.');
% grid on;
% xlabel('noise');
% title('Independence');
% hh1(1).LineWidth = 1.5; 
% hh1(2).LineWidth = 1.5; 
% hh1(3).LineWidth = 1.5; 

%% HI-FREQ SIN VALUE VERIFICATION

clear;
clc;

xMin = 0;
xMax = 1;
num_noise = 30;                    % The number of different noise levels used
noise = 3;                         % A constant to determine the amount of noise

MVecResults = [100:100:1000 2500 5000 10000];
num_noise_test_min = 0;
num_noise_test_max = 20;
noiseVec = num_noise_test_min:num_noise_test_max;
numMCSim = 200;

minscanincrVal = 0.015625;

cimVersion = 4;

switch(cimVersion)
    case 1
        fnameStr = 'cimv1';
        cimfunc = @cim;
    case 3
        fnameStr = 'cimv3';
        cimfunc = @cim_v3;
    case 4
        fnameStr = 'cimv4';
        cimfunc = @cim_v4;
    case 5
        fnameStr = 'cimv5';
        cimfunc = @cim_v5;
    case 6
        fnameStr = 'cimv6';
        cimfunc = @cim_v6;
    case 7
        fnameStr = 'cimv7';
        cimfunc = @cim_v7;
    otherwise
        error('Unrecognized CIM version!');
end
        
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');

for M=MVecResults
    % Test Empirical copula distance as a distance measure under various
    % amounts of noise
    hiFreqSinDep  = zeros(2,length(noiseVec));
                  
    for l=noiseVec
        dispstat(sprintf('Computing for noise level=%d >> M=%d',l,M),'keepthis', 'timestamp');
        t5 = 0; c5 = 0; c55 = 0;
        parfor mcSimNum=1:numMCSim
            x = rand(M,1)*(xMax-xMin)+xMin;
            y5 = sin(16*pi*x) + noise*(l/num_noise)*randn(M,1);

            U5 = pobs([x y5]);
            t5 = t5 + abs(corr(x,y5,'type','kendall'));
            zz1 = 0; zz2 = 0.02900;   r1 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w1 = zz2-zz1;
            zz1 = zz2; zz2 = 0.09158; r2 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w2 = zz2-zz1;
            zz1 = zz2; zz2 = 0.1528;  r3 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w3 = zz2-zz1;
            zz1 = zz2; zz2 = 0.2108;  r4 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w4 = zz2-zz1;
            zz1 = zz2; zz2 = 0.2715;  r5 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w5 = zz2-zz1;
            zz1 = zz2; zz2 = 0.3339;  r6 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w6 = zz2-zz1;
            zz1 = zz2; zz2 = 0.3999;  r7 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w7 = zz2-zz1;
            zz1 = zz2; zz2 = 0.4613;  r8 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w8 = zz2-zz1;
            zz1 = zz2; zz2 = 0.5227;  r9 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w9 = zz2-zz1;
            zz1 = zz2; zz2 = 0.5893;  r10 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w10 = zz2-zz1;
            zz1 = zz2; zz2 = 0.6547;  r11 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w11 = zz2-zz1;
            zz1 = zz2; zz2 = 0.7135;  r12 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w12 = zz2-zz1;
            zz1 = zz2; zz2 = 0.7734;  r13 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w13 = zz2-zz1;
            zz1 = zz2; zz2 = 0.8438;  r14 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w14 = zz2-zz1;
            zz1 = zz2; zz2 = 0.9074;  r15 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w15 = zz2-zz1;
            zz1 = zz2; zz2 = 0.9694;  r16 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w16 = zz2-zz1;
            zz1 = zz2; zz2 = 1;       r17 = inBoundedPts(U5(:,1),U5(:,2),zz1,zz2,0,1); w17 = zz2-zz1;
            c5 = c5 + ( abs(corr(r1(:,1),r1(:,2),'type','kendall'))*w1 + ...
                        abs(corr(r2(:,1),r2(:,2),'type','kendall'))*w2 + ...
                        abs(corr(r3(:,1),r3(:,2),'type','kendall'))*w3 + ...
                        abs(corr(r4(:,1),r4(:,2),'type','kendall'))*w4 + ...
                        abs(corr(r5(:,1),r5(:,2),'type','kendall'))*w5 + ...
                        abs(corr(r6(:,1),r6(:,2),'type','kendall'))*w6 + ...
                        abs(corr(r7(:,1),r7(:,2),'type','kendall'))*w7 + ...
                        abs(corr(r8(:,1),r8(:,2),'type','kendall'))*w8 + ...
                        abs(corr(r9(:,1),r9(:,2),'type','kendall'))*w9 + ...
                        abs(corr(r10(:,1),r10(:,2),'type','kendall'))*w10 + ...
                        abs(corr(r11(:,1),r11(:,2),'type','kendall'))*w11 + ...
                        abs(corr(r12(:,1),r12(:,2),'type','kendall'))*w12 + ...
                        abs(corr(r13(:,1),r13(:,2),'type','kendall'))*w13 + ...
                        abs(corr(r14(:,1),r14(:,2),'type','kendall'))*w14 + ...
                        abs(corr(r15(:,1),r15(:,2),'type','kendall'))*w15 + ...
                        abs(corr(r16(:,1),r16(:,2),'type','kendall'))*w16 + ...
                        abs(corr(r17(:,1),r17(:,2),'type','kendall'))*w17);
        end
        hiFreqSinDep(1,l+1) = t5/numMCSim; hiFreqSinDep(2,l+1) = c5/numMCSim;
    end

    % save the data
    if(ispc)
        save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\regionTheoretical_hiFreqSin_M_%d.mat', M));
    elseif(ismac)
        save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/clustering/regionTheoretical_hiFreqSin_M_%d.mat', M));
    else
        save(sprintf('/home/kiran/ownCloud/PhD/sim_results/clustering/regionTheoretical_hiFreqSin_M_%d', M));
    end
end

%% Plot the HI-FREQ sinusoidal stuff to see if our simulations match up with intuition
clear;
clc;

M = 1000;
if(ispc)
    load(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\regionTheoretical_hiFreqSin_M_%d.mat', M));
elseif(ismac)
    load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/clustering/regionTheoretical_hiFreqSin_M_%d.mat', M));
else
    load(sprintf('/home/kiran/ownCloud/PhD/sim_results/clustering/regionTheoretical_hiFreqSin_M_%d', M));
end

plot(noiseVec, hiFreqSinDep(2,:)); grid on;
xlabel('noise'); ylabel('CIM');