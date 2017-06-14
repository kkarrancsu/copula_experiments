%% Test the overlap probability computation
clear;
clc;

mu1 = 5.28; sigma1 = 0.91;
mu2 = 8.45; sigma2 = 1.36;

ovlp1 = computeOvlpProb(mu1,sigma1,mu1,sigma1);  % should be exactly 0
ovlp2 = computeOvlpProb(mu1,sigma1,mu2,sigma2);  % should be close to 0.158
ovlp3 = computeOvlpProb(2,2,1,1);  % should be close to 0.6
ovlp4 = computeOvlpProb(0,2,1,2);  % should be close to 0.8

fprintf('ovlp1=%0.02f ovlp2=%0.02f ovlp3=%0.02f ovlp4=%0.02f\n', ...
    ovlp1, ovlp2, ovlp3, ovlp4);

%% test the probabilistic region finder
clear;
clc;
dbstop if error;

numMCSims = 200;
% MVec = 50:50:500;
MVec = [100,250,500];  % only compute what we will plot, for now
xMin = 0;
xMax = 1;
num_noise = 30;                    % The number of different noise levels used
noise = 3; 
num_noise_test_min = 0;
num_noise_test_max = 30;
% noiseVec = num_noise_test_min:num_noise_test_max;
noiseVec = [0,10,20];   % only compute what we will plot, for now

linearResultsCell = cell(length(MVec),length(noiseVec),numMCSims);
quadraticResultsCell = cell(length(MVec),length(noiseVec),numMCSims);

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');
for mIdx=1:length(MVec)
    M = MVec(mIdx);
    for noiseIdx=1:length(noiseVec)
        l = noiseVec(noiseIdx);
        dispstat(sprintf('Computing for M=%d noise=%d', M, l),'keepthis', 'timestamp');
        % simulate data under the null w/ correct marginals
        parfor ii=1:numMCSims
            x = rand(M,1)*(xMax-xMin)+xMin;
            y1 = x + noise*(l/num_noise)*randn(M,1);
            y2 = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
            
            linearCell = probabilistic_region_finder(x,y1);
            linearResultsCell{mIdx,noiseIdx,ii} = linearCell;
            quadraticCell = probabilistic_region_finder(x,y2);
            quadraticResultsCell{mIdx,noiseIdx,ii} = quadraticCell;
        end
    end
end

% store the results
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\probabilistic_region_finder.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/clustering/probabilistic_region_finder.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/clustering/probabilistic_region_finder.mat');
end

%% Plot results for the probabilistic region finder
clear;
clc;
close all;
dbstop if error;

if(ispc)
    load('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\probabilistic_region_finder.mat');
elseif(ismac)
    load('/Users/Kiran/ownCloud/PhD/sim_results/clustering/probabilistic_region_finder.mat');
else
    load('/home/kiran/ownCloud/PhD/sim_results/clustering/probabilistic_region_finder.mat');
end

% which sample sizes to plot
MVecToPlot = [100, 500];
noiseLevelsToPlot = [0,10];
scanincrsToPlot = [20];
subplotCfg = [length(MVecToPlot), length(noiseLevelsToPlot)];
numPlots = prod(subplotCfg);
lineMarkers = {'+-.','o-.','*-.','d-.','x-.','s.-'};

numDeps = 2;
figureVec = zeros(1,numDeps*2);
for ii=1:numDeps*2
    figureVec(ii) = figure(ii);
end

plotIdx = 1;
numToAvg = numMCSims;
% numToAvg = 1;
for MIdx=1:length(MVecToPlot)
    % find which index this goes into the resultsMat
    M = MVecToPlot(MIdx);
    mIdx = find(M==MVec);
    
    for noiseLevelIdx=1:length(noiseLevelsToPlot)
        legendCell1 = {}; legendCell1Idx = 1;
        legendCell2 = {}; legendCell2Idx = 1;
        
        noiseLevel = noiseLevelsToPlot(noiseLevelIdx);
        lIdx = find(noiseLevel==noiseVec);
        % get the data        
        avgResults_linear = ...
            getSignatureAvgRegionFinder(linearResultsCell,mIdx,lIdx,numToAvg,scanincrsToPlot);
        avgResults_parabola = ...
            getSignatureAvgRegionFinder(quadraticResultsCell,mIdx,lIdx,numToAvg,scanincrsToPlot);
        
        set(0,'CurrentFigure',figureVec(1));
        B = [subplotCfg plotIdx]; B = mat2cell(B,1,ones(1,numel(B)));
        subplot(B{:});
        for ii=1:length(scanincrsToPlot)
            scanincr = scanincrsToPlot(ii);
            y = avgResults_linear{ii};
            y1 = y(7,:); x1 = linspace(scanincr/M,1,length(y1)); y1_err = y(8,:);
            y2 = y(9,:); x2 = linspace(scanincr/M,1,length(y2)); y2_err = y(10,:);
            y3 = y(11,:); x3 = linspace(scanincr/M,1,length(y3)); y3_err = y(12,:);
            errorbar(x1,y1,y1_err,lineMarkers{ii},'LineWidth',3);
            hold on; grid on;
            errorbar(x2,y2,y2_err,lineMarkers{ii},'LineWidth',3);
            errorbar(x3,y3,y3_err,lineMarkers{ii},'LineWidth',3);
            legendCell1{legendCell1Idx} = sprintf('\\mu(R1) - \\Delta i=%d', scanincrsToPlot(ii));
            legendCell1Idx = legendCell1Idx + 1;
            legendCell1{legendCell1Idx} = sprintf('\\mu(R2) - \\Delta i=%d', scanincrsToPlot(ii));
            legendCell1Idx = legendCell1Idx + 1;
            legendCell1{legendCell1Idx} = sprintf('\\mu(R3) - \\Delta i=%d', scanincrsToPlot(ii));
            legendCell1Idx = legendCell1Idx + 1;
        end
        title(sprintf('M=%d \\sigma^2=%d',M,noiseLevel));
        if(MIdx==1 && noiseLevelIdx==1)
            legend(legendCell1,'Location','SouthWest');
        end
        
        set(0,'CurrentFigure',figureVec(1+numDeps));
        B = [subplotCfg plotIdx]; B = mat2cell(B,1,ones(1,numel(B)));
        subplot(B{:});
        for ii=1:length(scanincrsToPlot)
            scanincr = scanincrsToPlot(ii);
            y = avgResults_linear{ii};
            y1 = y(1,:); x1 = linspace(scanincr/M,1,length(y1));
            y2 = y(3,:); x2 = linspace(scanincr/M,1,length(y2));
            y3 = y(5,:); x3 = linspace(scanincr/M,1,length(y3));
            plot(x1,y1,lineMarkers{ii},'LineWidth',3);
            hold on; grid on;
            plot(x2,y2,lineMarkers{ii},'LineWidth',3);
            plot(x3,y3,lineMarkers{ii},'LineWidth',3);
            legendCell2{legendCell2Idx} = sprintf('\\Delta(R1,R2) - \\Delta i=%d', scanincrsToPlot(ii));
            legendCell2Idx = legendCell2Idx + 1;
            legendCell2{legendCell2Idx} = sprintf('\\Delta(R3,R2) - \\Delta i=%d', scanincrsToPlot(ii));
            legendCell2Idx = legendCell2Idx + 1;
            legendCell2{legendCell2Idx} = sprintf('\\Delta(R1,R3) - \\Delta i=%d', scanincrsToPlot(ii));
            legendCell2Idx = legendCell2Idx + 1;
        end
        title(sprintf('M=%d \\sigma^2=%d',M,noiseLevel));
        if(MIdx==1 && noiseLevelIdx==1)
            legend(legendCell2,'Location','SouthWest');
        end
        
        set(0,'CurrentFigure',figureVec(2));
        subplot(B{:});
        % compute and plot the z-score from the metric & numPts
        for ii=1:length(scanincrsToPlot)
            scanincr = scanincrsToPlot(ii);
            y = avgResults_parabola{ii};
            y1 = y(7,:); x1 = linspace(scanincr/M,1,length(y1)); y1_err = y(8,:);
            y2 = y(9,:); x2 = linspace(scanincr/M,1,length(y2)); y2_err = y(10,:);
            y3 = y(11,:); x3 = linspace(scanincr/M,1,length(y3)); y3_err = y(12,:);
            errorbar(x1,y1,y1_err,lineMarkers{ii},'LineWidth',3);
            hold on; grid on;
            errorbar(x2,y2,y2_err,lineMarkers{ii},'LineWidth',3);
            errorbar(x3,y3,y3_err,lineMarkers{ii},'LineWidth',3);
            hold on; grid on;
        end
        title(sprintf('M=%d \\sigma^2=%d',M,noiseLevel));
        if(MIdx==1 && noiseLevelIdx==1)
            legend(legendCell1,'Location','SouthWest');
        end
        
        set(0,'CurrentFigure',figureVec(2+numDeps));
        subplot(B{:});
        % compute and plot the z-score from the metric & numPts
        for ii=1:length(scanincrsToPlot)
            scanincr = scanincrsToPlot(ii);
            y = avgResults_parabola{ii};
            y1 = y(1,:); x1 = linspace(scanincr/M,1,length(y1));
            y2 = y(3,:); x2 = linspace(scanincr/M,1,length(y2));
            y3 = y(5,:); x3 = linspace(scanincr/M,1,length(y3));
            plot(x1,y1,lineMarkers{ii},'LineWidth',3);
            hold on; grid on;
            plot(x2,y2,lineMarkers{ii},'LineWidth',3);
            plot(x3,y3,lineMarkers{ii},'LineWidth',3);
            hold on; grid on;
        end
        title(sprintf('M=%d \\sigma^2=%d',M,noiseLevel));
        if(MIdx==1 && noiseLevelIdx==1)
            legend(legendCell2,'Location','SouthWest');
        end
        
        plotIdx = plotIdx + 1;
    end
    set(0,'CurrentFigure',figureVec(1)); figtitle('Linear');
    set(0,'CurrentFigure',figureVec(1+numDeps)); figtitle('Linear');
    set(0,'CurrentFigure',figureVec(2)); figtitle('Parabolic');
    set(0,'CurrentFigure',figureVec(2+numDeps)); figtitle('Parabolic');
end

%% Test different region detection mechanisms by creating different data in the
% copula space

clear;
clc;
close all;

numMCSims = 500;
MVec = [50 100 200:100:1000];
noiseVec = 1:30;
% simulate different change rates
cornerPts = 0.1:0.1:0.9;
scanincrs = [1 0.5 0.25 0.125 0.0625 0.03125 0.015625];

rd1Cell_up = cell(length(MVec),length(noiseVec),length(cornerPts),length(scanincrs),numMCSims);
rd3Cell_up = cell(length(MVec),length(noiseVec),length(cornerPts),length(scanincrs),numMCSims);
rd4Cell_up = cell(length(MVec),length(noiseVec),length(cornerPts),length(scanincrs),numMCSims);
rd5Cell_up = cell(length(MVec),length(noiseVec),length(cornerPts),length(scanincrs),numMCSims);
rd6Cell_up = cell(length(MVec),length(noiseVec),length(cornerPts),length(scanincrs),numMCSims);
rd7Cell_up = cell(length(MVec),length(noiseVec),length(cornerPts),length(scanincrs),numMCSims);

rd1Cell_down = cell(length(MVec),length(noiseVec),length(cornerPts),length(scanincrs),numMCSims);
rd3Cell_down = cell(length(MVec),length(noiseVec),length(cornerPts),length(scanincrs),numMCSims);
rd4Cell_down = cell(length(MVec),length(noiseVec),length(cornerPts),length(scanincrs),numMCSims);
rd5Cell_down = cell(length(MVec),length(noiseVec),length(cornerPts),length(scanincrs),numMCSims);
rd6Cell_down = cell(length(MVec),length(noiseVec),length(cornerPts),length(scanincrs),numMCSims);
rd7Cell_down = cell(length(MVec),length(noiseVec),length(cornerPts),length(scanincrs),numMCSims);

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');
for mIdx=1:length(MVec)
    M = MVec(mIdx);
    for noiseIdx=1:length(noiseVec)
        l = noiseVec(noiseIdx);
        dispstat(sprintf('Computing for M=%d noise=%d', M, l),'keepthis', 'timestamp');
        for cornerPtIdx=1:length(cornerPts)
            cornerPt = cornerPts(cornerPtIdx);
            for scanincrIdx=1:length(scanincrs)
                scanincr = scanincrs(scanincrIdx);
                % simulate data under the null w/ correct marginals
                parfor ii=1:numMCSims
                    maxU = min(1,2*cornerPt);
                    u = linspace(0,maxU,M)';
                    vv1 = linspace(1,0,M/2); vv2 = linspace(0,1,M/2);
                    v1 = [vv1 vv2]'; v2 = [vv2 vv1]';  % test both directions, down --> up and down --> up
                    
                    bp1_1 = rd1(u,v1,scanincr,120,0.05);    bp1_2 = rd1(u,v2,scanincr,120,0.05);
                    bp3_1 = rd3(u,v1,scanincr,120);         bp3_2 = rd3(u,v2,scanincr,120);
                    bp4_1 = rd4(u,v1,scanincr);             bp4_2 = rd4(u,v2,scanincr);
                    bp5_1 = rd5(u,v1,scanincr);             bp5_2 = rd5(u,v2,scanincr);
                    bp6_1 = rd6(u,v1,scanincr);             bp6_2 = rd6(u,v2,scanincr);
                    bp7_1 = rd7(u,v1,scanincr);             bp7_2 = rd7(u,v2,scanincr);
                    
                    rd1Cell_up{mIdx,noiseIdx,cornerPtIdx,scanincrIdx,ii} = bp1_1;
                    rd1Cell_down{mIdx,noiseIdx,cornerPtIdx,scanincrIdx,ii} = bp1_2;
                    
                    rd3Cell_up{mIdx,noiseIdx,cornerPtIdx,scanincrIdx,ii} = bp3_1;
                    rd3Cell_down{mIdx,noiseIdx,cornerPtIdx,scanincrIdx,ii} = bp3_2;
                    
                    rd4Cell_up{mIdx,noiseIdx,cornerPtIdx,scanincrIdx,ii} = bp4_1;
                    rd4Cell_down{mIdx,noiseIdx,cornerPtIdx,scanincrIdx,ii} = bp4_2;
                    
                    rd5Cell_up{mIdx,noiseIdx,cornerPtIdx,scanincrIdx,ii} = bp5_1;
                    rd5Cell_down{mIdx,noiseIdx,cornerPtIdx,scanincrIdx,ii} = bp5_2;
                    
                    rd6Cell_up{mIdx,noiseIdx,cornerPtIdx,scanincrIdx,ii} = bp6_1;
                    rd6Cell_down{mIdx,noiseIdx,cornerPtIdx,scanincrIdx,ii} = bp6_2;
                    
                    rd7Cell_up{mIdx,noiseIdx,cornerPtIdx,scanincrIdx,ii} = bp7_1;
                    rd7Cell_down{mIdx,noiseIdx,cornerPtIdx,scanincrIdx,ii} = bp7_2;
                end
            end
        end
    end
end

% store the results
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\region_detection_tests.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/clustering/region_detection_tests.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/clustering/region_detection_tests.mat');
end

%% Plot & Interpret the results

clear;
clc;
close all;

% load the results
if(ispc)
    load('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\region_detection_tests.mat');
elseif(ismac)
    load('/Users/Kiran/ownCloud/PhD/sim_results/clustering/region_detection_tests.mat');
else
    load('/home/kiran/ownCloud/PhD/sim_results/clustering/region_detection_tests.mat');
end
