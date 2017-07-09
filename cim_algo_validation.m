%% Some experiments to validate the CIM algorithm

%% check kendall's tau's distribution for functionally dependent but P(conc)=P(disc) relationships
clear;
clc;

numMCSims = 500;
MVec = 50:50:500;
xMin = 0;
xMax = 1;
num_noise = 30;                    % The number of different noise levels used
noise = 3; 
num_noise_test_min = 0;
num_noise_test_max = 30;
noiseVec = num_noise_test_min:num_noise_test_max;

resultsMat = zeros(length(MVec),length(noiseVec),numMCSims);

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
            y2 = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
            % compute tau so we can look at the distribution of it, to see if
            % it matches the N(0,2(2n+5)/(9n(n+1))) as proved by Kendall
            tauVal = corr(x,y2,'type','kendall');
            resultsMat(mIdx,noiseIdx,ii) = tauVal;
        end
    end
end

% store the results
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\tau_distribution.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/clustering/tau_distribution.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/clustering/tau_distribution.mat');
end

%% analyze the tau distribution
clear;
clc;
 
 % store the results
if(ispc)
    load('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\tau_distribution.mat');
elseif(ismac)
    load('/Users/Kiran/ownCloud/PhD/sim_results/clustering/tau_distribution.mat');
else
    load('/home/kiran/ownCloud/PhD/sim_results/clustering/tau_distribution.mat');
end

% which sample sizes to plot
MVecToPlot = [50, 250, 500];
noiseLevelsToPlot = [0,5,10,20];

% plot the empirical distribution of Kendall's Tau, and overlay the
% theoretical curve
numPlots = length(MVecToPlot);
p = numSubplots(numPlots);
for plotIdx=1:numPlots
    % get the data
    M = MVecToPlot(plotIdx);
    % find which index this goes into the resultsMat
    mIdx = find(M==MVec);
    
    % find the indices for the noiseLevelsToPlot in the noiseVec
    legendCell = cell(1,length(noiseLevelsToPlot));
    legendCellIdx = 1;
    minXi = 1000; maxXi = -1000;
    for noiseLevel=noiseLevelsToPlot
        lIdx = find(noiseLevel==noiseVec);
        % get the data
        tauData = squeeze(resultsMat(mIdx,lIdx,:));
        
        B = [p plotIdx]; B = mat2cell(B,1,ones(1,numel(B)));
        subplot(B{:});
        [f,xi] = ksdensity(tauData);
        plot(xi,f,'-');
        hold on;
        
        if(min(xi)<minXi)
            minXi = min(xi);
        end
        if(max(xi)>maxXi)
            maxXi = max(xi);
        end
        
        % add to the legendCell
        legendCell{legendCellIdx} = sprintf('\\sigma^2=%d',noiseLevel);
        legendCellIdx = legendCellIdx + 1;
    end
    % overlay teh theoretical normal distribution of what we expect to
    % see
    y2 = normpdf(xi,0,(2*(2*M+5))/(9*M*(M-1)));
    plot(xi,y2,'-.');
    legendCell{legendCellIdx} = 'theoretical';
    
    title(sprintf('M=%d',M));
    grid on; xlabel('\tau');
    legend(legendCell);
    
end

%% test the CIM region finder
clear;
clc;

numMCSims = 250;
MVec = [50, 100, 500];
xMin = 0;
xMax = 1;
num_noise = 30;                    % The number of different noise levels used
noise = 3; 
% noiseVec = num_noise_test_min:num_noise_test_max;
noiseVec = [0,10,20];
minscanincr = 0.015;  % ends up that .015625 is the last one that gets run

resultsMat_linear   = cell(length(MVec),length(noiseVec),numMCSims);
resultsMat_parabola = cell(length(MVec),length(noiseVec),numMCSims);
resultsMat_cubic    = cell(length(MVec),length(noiseVec),numMCSims);
resultsMat_sinu     = cell(length(MVec),length(noiseVec),numMCSims);
resultsMat_hfsinu   = cell(length(MVec),length(noiseVec),numMCSims);
resultsMat_fr       = cell(length(MVec),length(noiseVec),numMCSims);
resultsMat_step     = cell(length(MVec),length(noiseVec),numMCSims);

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
            y3 = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
            y4 = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
            y5 = sin(16*pi*x) + noise*(l/num_noise)*randn(M,1);
            y6 = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
            y8 = (x > 0.5) + noise*5*l/num_noise*randn(M,1);
            
            signature_na_linear   = cim_region_finder(x,y1,minscanincr);
            signature_na_parabola = cim_region_finder(x,y2,minscanincr);
            signature_na_cubic    = cim_region_finder(x,y3,minscanincr);
            signature_na_sinu     = cim_region_finder(x,y4,minscanincr);
            signature_na_hfsinu   = cim_region_finder(x,y5,minscanincr);
            signature_na_fr       = cim_region_finder(x,y6,minscanincr);
            signature_na_step     = cim_region_finder(x,y8,minscanincr);
            
            resultsMat_linear{mIdx,noiseIdx,ii}   = signature_na_linear;
            resultsMat_parabola{mIdx,noiseIdx,ii} = signature_na_parabola;
            resultsMat_cubic{mIdx,noiseIdx,ii} = signature_na_cubic;
            resultsMat_sinu{mIdx,noiseIdx,ii} = signature_na_sinu;
            resultsMat_hfsinu{mIdx,noiseIdx,ii} = signature_na_hfsinu;
            resultsMat_fr{mIdx,noiseIdx,ii} = signature_na_fr;
            resultsMat_step{mIdx,noiseIdx,ii} = signature_na_step;
        end
    end
end
dispstat('Simulation complete!','keepthis', 'timestamp');
% store the results
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\region_finder.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/clustering/region_finder.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/clustering/region_finder.mat');
end

%% analyze the signature
clear;
clc;
close all;
dbstop if error;
 
 % store the results
if(ispc)
    load('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\region_finder.mat');
elseif(ismac)
    load('/Users/Kiran/ownCloud/PhD/sim_results/clustering/region_finder.mat');
else
    load('/home/kiran/ownCloud/PhD/sim_results/clustering/region_finder.mat');
end

% which sample sizes to plot
MVecToPlot = [50, 100, 500];
noiseLevelsToPlot = [0,10,20];
% scanincrsToPlot = [0.5, 0.25, 0.125, .0625, .03125];
scanincrsToPlot = [0.5, .25, 0.125, .0625, .03125];
subplotCfg = [length(MVecToPlot), length(noiseLevelsToPlot)];
numPlots = prod(subplotCfg);
lineMarkers = {'+-.','o-.','*-.','d-.','x-.','s.-'};

numDeps = 7;
figureVec = zeros(1,numDeps);
for ii=1:numDeps
    figureVec(ii) = figure(ii);
end

plotIdx = 1;
% numToAvg = numMCSims;
numToAvg = 1;
for MIdx=1:length(MVecToPlot)
    % find which index this goes into the resultsMat
    M = MVecToPlot(MIdx);
    mIdx = find(M==MVec);
    
    for noiseLevelIdx=1:length(noiseLevelsToPlot)
        legendCell1 = {}; legendCell1Idx = 1;
        
        noiseLevel = noiseLevelsToPlot(noiseLevelIdx);
        lIdx = find(noiseLevel==noiseVec);
        % get the data
        [avgMetric_linear, avgNumPts_linear] = ...
            getSignatureAvg(resultsMat_linear,mIdx,lIdx,numToAvg,scanincrsToPlot);
        [avgMetric_parabola, avgNumPts_parabola] = ...
            getSignatureAvg(resultsMat_parabola,mIdx,lIdx,numToAvg,scanincrsToPlot);
        [avgMetric_cubic, avgNumPts_cubic] = ...
            getSignatureAvg(resultsMat_cubic,mIdx,lIdx,numToAvg,scanincrsToPlot);
        [avgMetric_sinu, avgNumPts_sinu] = ...
            getSignatureAvg(resultsMat_sinu,mIdx,lIdx,numToAvg,scanincrsToPlot);
        [avgMetric_hfsinu, avgNumPts_hfsinu] = ...
            getSignatureAvg(resultsMat_hfsinu,mIdx,lIdx,numToAvg,scanincrsToPlot);
        [avgMetric_fr, avgNumPts_fr] = ...
            getSignatureAvg(resultsMat_fr,mIdx,lIdx,numToAvg,scanincrsToPlot);
        [avgMetric_step, avgNumPts_step] = ...
            getSignatureAvg(resultsMat_step,mIdx,lIdx,numToAvg,scanincrsToPlot);
        
        set(0,'CurrentFigure',figureVec(1));
        B = [subplotCfg plotIdx]; B = mat2cell(B,1,ones(1,numel(B)));
        subplot(B{:});
        for ii=1:length(scanincrsToPlot)
            scanincr = scanincrsToPlot(ii);
            n = avgNumPts_linear{ii};
            den = (2*(2*n+5))./(9*n.*(n-1));
            y = avgMetric_linear{ii} ./ den;
            x = linspace(0,1,length(y));
            plot(x,y,lineMarkers{ii},'LineWidth',3);
            hold on; grid on;
            legendCell1{legendCell1Idx} = sprintf('NA - \\Delta i=%0.03f', scanincrsToPlot(ii));
            legendCell1Idx = legendCell1Idx + 1;
        end
        title(sprintf('M=%d \\sigma^2=%d',M,noiseLevel));
        if(MIdx==1 && noiseLevelIdx==1)
            legend(legendCell1,'Location','SouthWest');
        end
        
        set(0,'CurrentFigure',figureVec(2));
        subplot(B{:});
        % compute and plot the z-score from the metric & numPts
        for ii=1:length(scanincrsToPlot)
            scanincr = scanincrsToPlot(ii);
            n = avgNumPts_parabola{ii};
            den = (2*(2*n+5))./(9*n.*(n-1));
            y = avgMetric_parabola{ii} ./ den;
            x = linspace(0,1,length(y));
            plot(x,y,lineMarkers{ii},'LineWidth',3);
            hold on; grid on;
        end
        title(sprintf('M=%d \\sigma^2=%d',M,noiseLevel));
        if(MIdx==1 && noiseLevelIdx==1)
            legend(legendCell1,'Location','SouthWest');
        end
        
        set(0,'CurrentFigure',figureVec(3));
        subplot(B{:});
        for ii=1:length(scanincrsToPlot)
            scanincr = scanincrsToPlot(ii);
            n = avgNumPts_cubic{ii};
            den = (2*(2*n+5))./(9*n.*(n-1));
            y = avgMetric_cubic{ii} ./ den;
            x = linspace(0,1,length(y));
            plot(x,y,lineMarkers{ii},'LineWidth',3);
            hold on; grid on;
        end
        title(sprintf('M=%d \\sigma^2=%d',M,noiseLevel));
        if(MIdx==1 && noiseLevelIdx==1)
            legend(legendCell1,'Location','SouthWest');
        end
        
        set(0,'CurrentFigure',figureVec(4));
        subplot(B{:});
        for ii=1:length(scanincrsToPlot)
            scanincr = scanincrsToPlot(ii);
            n = avgNumPts_sinu{ii};
            den = (2*(2*n+5))./(9*n.*(n-1));
            y = avgMetric_sinu{ii} ./ den;
            x = linspace(0,1,length(y));
            plot(x,y,lineMarkers{ii},'LineWidth',3);
            hold on; grid on;
        end
        title(sprintf('M=%d \\sigma^2=%d',M,noiseLevel));
        if(MIdx==1 && noiseLevelIdx==1)
            legend(legendCell1,'Location','SouthWest');
        end
        
        set(0,'CurrentFigure',figureVec(5));
        subplot(B{:});
        for ii=1:length(scanincrsToPlot)
            scanincr = scanincrsToPlot(ii);
            n = avgNumPts_hfsinu{ii};
            den = (2*(2*n+5))./(9*n.*(n-1));
            y = avgMetric_hfsinu{ii} ./ den;
            x = linspace(0,1,length(y));
            plot(x,y,lineMarkers{ii},'LineWidth',3);
            hold on; grid on;
        end
        title(sprintf('M=%d \\sigma^2=%d',M,noiseLevel));
        if(MIdx==1 && noiseLevelIdx==1)
            legend(legendCell1,'Location','SouthWest');
        end
        
        set(0,'CurrentFigure',figureVec(6));
        subplot(B{:});
        for ii=1:length(scanincrsToPlot)
            scanincr = scanincrsToPlot(ii);
            n = avgNumPts_fr{ii};
            den = (2*(2*n+5))./(9*n.*(n-1));
            y = avgMetric_fr{ii} ./ den;
            x = linspace(0,1,length(y));
            plot(x,y,lineMarkers{ii},'LineWidth',3);
            hold on; grid on;
        end
        title(sprintf('M=%d \\sigma^2=%d',M,noiseLevel));
        if(MIdx==1 && noiseLevelIdx==1)
            legend(legendCell1,'Location','SouthWest');
        end
        
        set(0,'CurrentFigure',figureVec(7));
        subplot(B{:});
        for ii=1:length(scanincrsToPlot)
            scanincr = scanincrsToPlot(ii);
            n = avgNumPts_step{ii};
            den = (2*(2*n+5))./(9*n.*(n-1));
            y = avgMetric_step{ii} ./ den;
            x = linspace(0,1,length(y));
            plot(x,y,lineMarkers{ii},'LineWidth',3);
            hold on; grid on;
        end
        title(sprintf('M=%d \\sigma^2=%d',M,noiseLevel));
        if(MIdx==1 && noiseLevelIdx==1)
            legend(legendCell1,'Location','SouthWest');
        end
        
        plotIdx = plotIdx + 1;
    end
    set(0,'CurrentFigure',figureVec(1)); figtitle('Linear');
    set(0,'CurrentFigure',figureVec(2)); figtitle('Parabolic');
    set(0,'CurrentFigure',figureVec(3)); figtitle('Cubic');
    set(0,'CurrentFigure',figureVec(4)); figtitle('Sinusoidal');
    set(0,'CurrentFigure',figureVec(5)); figtitle('HF-Sine');
    set(0,'CurrentFigure',figureVec(6)); figtitle('Fourth-Root');
    set(0,'CurrentFigure',figureVec(7)); figtitle('Step-Function');
end

%% compare the old CIM to the new CIM
clear; clc;

M = 500; xMin = 0; xMax = 1; num_noise = 30; noise = 3; l = 0;
x = rand(M,1)*(xMax-xMin)+xMin;
y1 = x + noise*(l/num_noise)*randn(M,1);
y2 = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
y3 = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
y4 = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
y5 = sin(16*pi*x) + noise*(l/num_noise)*randn(M,1);
y6 = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
y8 = (x > 0.5) + noise*5*l/num_noise*randn(M,1);
y9 =  rand(M,1)*(xMax-xMin)+xMin;

y_sig = y4;

cim(x,y_sig)
cim_v2(x,y_sig)

%% Generate statistical power curves for old CIM vs new CIM
% same methodology as Simon & Tibshirani:
% http://statweb.stanford.edu/~tibs/reshef/script.R

clear;
clc;

rng(1234);
dbstop if error;

nsim_null = 200;   % The number of null datasets we use to estimate our rejection reject regions for an alternative with level 0.05
nsim_alt  = 200;   % Number of alternative datasets we use to estimate our power

num_noise = 30;                    % The number of different noise levels used
noise = 3;                         % A constant to determine the amount of noise

M = 500;                % number of samples
numDepTests = 8;        % the number of different dependency tests we will conduct

minscanincrVal = 0.015625;
% Vectors holding the null "correlations" (for pearson, dcor and mic respectively) 
% for each of the nsim null datasets at a given noise level
cimNull = zeros(1,nsim_null);
cimv3Null = zeros(1,nsim_null);
cimv4MexNull = zeros(1,nsim_null);
cimv8aMexNull = zeros(1,nsim_null);
cimv8bMexNull = zeros(1,nsim_null);
cimv7Null = zeros(1,nsim_null);

cimAlt  = zeros(1,nsim_alt);
cimv3Alt = zeros(1,nsim_alt);
cimv4MexAlt = zeros(1,nsim_alt);
cimv8aMexAlt = zeros(1,nsim_alt);
cimv8bMexAlt = zeros(1,nsim_alt);
cimv7Alt = zeros(1,nsim_alt);

% Arrays holding the estimated power for each of the "correlation" types, 
% for each data type (linear, parabolic, etc...) with each noise level
cimPower = zeros(numDepTests, num_noise);
cimv3Power = zeros(numDepTests, num_noise);
cimv4MexPower = zeros(numDepTests, num_noise);
cimv8aMexPower = zeros(numDepTests,num_noise);
cimv8bMexPower = zeros(numDepTests,num_noise);
cimv7Power = zeros(numDepTests,num_noise);

% Simon & Tibshirani use xMin=0, xMax=1 for performing their analysis ...
xMin = 0;
xMax = 1;

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');
num_noise_test_min = 1;
num_noise_test_max = 30;
for l=num_noise_test_min:num_noise_test_max
    for typ=1:numDepTests
        dispstat(sprintf('Computing for noise level=%d Dependency Test=%d',l, typ),'keepthis', 'timestamp');
        % simulate data under the null w/ correct marginals
        parfor ii=1:nsim_null
            x = rand(M,1)*(xMax-xMin)+xMin;
            switch(typ)
                case 1
                    % linear
                    y = x + noise*(l/num_noise)*randn(M,1); 
                case 2
                    % parabolic
                    y = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
                case 3
                    % cubic
                    y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
                case 4
                    % low-freq sin
                    y = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
                case 5
                    % high-freq sin
                    y = sin(16*pi*x) + noise*(l/num_noise)*randn(M,1);
                case 6
                    % fourth root
                    y = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
                case 7
                    % circle
                    y=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
                case 8
                    % step function
                    y = (x > 0.5) + noise*5*l/num_noise*randn(M,1);
                otherwise
                    error('unknown dep type!');
            end
            % resimulate x so we have null scenario
            x = rand(M,1)*(xMax-xMin)+xMin;
            
            % calculate the metrics
            cimNull(ii)   = cim(x,y);
            cimv3Null(ii) = cim_v3(x,y,minscanincrVal);
            cimv4MexNull(ii) = cim_v4(x,y,minscanincrVal);
            cimv8aMexNull(ii) = cim_v5(x,y,minscanincrVal);
            cimv8bMexNull(ii) = cim_v6(x,y,minscanincrVal);
            cimv7Null(ii) = cim_v7(x,y,minscanincrVal);
        end
        
        % compute the rejection cutoffs
        cim_cut = quantile(cimNull, 0.95);
        cimv3_cut = quantile(cimv3Null, 0.95);
        cimv4Mex_cut = quantile(cimv4MexNull, 0.95);
        cimv8aMex_cut = quantile(cimv8aMexNull, 0.95);
        cimv8bMex_cut = quantile(cimv8bMexNull, 0.95);
        cimv7_cut = quantile(cimv7Null, 0.95);
        
        % resimulate the data under the alternative hypothesis
        parfor ii=1:nsim_alt
            x = rand(M,1)*(xMax-xMin)+xMin;
            switch(typ)
                case 1
                    % linear
                    y = x + noise*(l/num_noise)*randn(M,1); 
                case 2
                    % parabolic
                    y = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
                case 3
                    % cubic
                    y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3) + 10*noise*(l/num_noise)*randn(M,1);
                case 4
                    % low-freq sin
                    y = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
                case 5
                    % high-freq sin
                    y = sin(16*pi*x) + noise*(l/num_noise)*randn(M,1);
                case 6
                    % fourth root
                    y = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
                case 7
                    % circle
                    y=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
                case 8
                    % step function
                    y = (x > 0.5) + noise*5*l/num_noise*randn(M,1);
                otherwise
                    error('unknown dep type!');
            end
            
            % calculate the metrics
            cimAlt(ii)   = cim(x,y);
            cimv3Alt(ii) = cim_v3(x,y,minscanincrVal);
            cimv4MexAlt(ii) = cim_v4(x,y,minscanincrVal);
            cimv8aMexAlt(ii) = cim_v5(x,y,minscanincrVal);
            cimv8bMexAlt(ii) = cim_v6(x,y,minscanincrVal);
            cimv7Alt(ii) = cim_v7(x,y,minscanincrVal);
        end
        
        % compute the power
        cimPower(typ, l)   = sum(cimAlt > cim_cut)/nsim_alt;
        cimv3Power(typ, l)  = sum(cimv3Alt > cimv3_cut)/nsim_alt;
        cimv4MexPower(typ, l) = sum(cimv4MexAlt > cimv4Mex_cut)/nsim_alt;
        cimv8aMexPower(typ, l)  = sum(cimv8aMexAlt > cimv8aMex_cut)/nsim_alt;
        cimv8bMexPower(typ, l)  = sum(cimv8bMexAlt > cimv8bMex_cut)/nsim_alt;
        cimv7Power(typ, l)  = sum(cimv7Alt > cimv7_cut)/nsim_alt;
    end
end

% save the data
if(ispc)
    save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\cim_vstar_power_M_%d.mat', M));
elseif(ismac)
    save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/cim_vstar_power_M_%d.mat', M));
else
    save(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/cim_vstar_power_M_%d.mat', M));
end
%%
clear;
clc;
close all;
dbstop if error;

M = 500;  % which one do we want to plot?

if(ispc)
    load(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\cim_vstar_power_M_%d.mat', M));
elseif(ismac)
    load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/cim_vstar_power_M_%d.mat', M));
else
    load(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/cim_vstar_power_M_%d.mat', M));
end

powerMat = zeros(6,8,length(num_noise_test_min:num_noise_test_max));
powerMat(1,:,:) = cimPower;
powerMat(2,:,:) = cimv3Power;
powerMat(3,:,:) = cimv4MexPower;
powerMat(4,:,:) = cimv8aMexPower;
powerMat(5,:,:) = cimv8bMexPower;
powerMat(6,:,:) = cimv7Power;
noiseVec = (num_noise_test_min:num_noise_test_max)/10;

labels = {'CIM', 'CIMv3', 'CIMv4', 'CIMv5', 'CIMv6', 'CIMv7'};
plotStyle = 1;
plotPower(powerMat, M, labels, noiseVec, num_noise_test_min, num_noise_test_max, plotStyle)


%% Generate statistical power curves for CIMv4 vs CIMv8a/b
% same methodology as Simon & Tibshirani:
% http://statweb.stanford.edu/~tibs/reshef/script.R

clear;
clc;

rng(1234);
% rng(12345);
dbstop if error;

nsim_null = 250;   % The number of null datasets we use to estimate our rejection reject regions for an alternative with level 0.05
nsim_alt  = 250;   % Number of alternative datasets we use to estimate our power

num_noise_test_min = 0;
num_noise_test_max = 20;
noiseVec = num_noise_test_min:num_noise_test_max;
num_noise = length(noiseVec);                    
                                   % The number of different noise levels used
noise = 3;                         % A constant to determine the amount of noise

M = 500;                % number of samples
numDepTests = 8;        % the number of different dependency tests we will conduct

% Vectors holding the null "correlations" (for pearson, dcor and mic respectively) 
% for each of the nsim null datasets at a given noise level
cimv4Null     = zeros(1,nsim_null);
cimv4MexNull = zeros(1,nsim_null);
cimv8aMexNull = zeros(1,nsim_null);
cimv8bMexNull = zeros(1,nsim_null);

cimv4Alt     = zeros(1,nsim_alt);
cimv4MexAlt = zeros(1,nsim_alt);
cimv8aMexAlt = zeros(1,nsim_alt);
cimv8bMexAlt = zeros(1,nsim_alt);

% configuration parameter for cim algorithm
minScanIncr=0.015625;

% Arrays holding the estimated power for each of the "correlation" types, 
% for each data type (linear, parabolic, etc...) with each noise level
cimv4Power     = zeros(numDepTests,num_noise);
cimv4MexPower  = zeros(numDepTests,num_noise);
cimv8aMexPower = zeros(numDepTests,num_noise);
cimv8bMexPower = zeros(numDepTests,num_noise);

% Simon & Tibshirani use xMin=0, xMax=1 for performing their analysis ...
xMin = 0;
xMax = 1;

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');
for lIdx=1:num_noise
    l = noiseVec(lIdx);
    for typ=1:numDepTests
        dispstat(sprintf('Computing for noise level=%d Dependency Test=%d',l, typ),'keepthis', 'timestamp');
        % simulate data under the null w/ correct marginals
        parfor ii=1:nsim_null
            x = rand(M,1)*(xMax-xMin)+xMin;
            switch(typ)
                case 1
                    % linear
                    y = x + noise*(l/num_noise)*randn(M,1); 
                case 2
                    % parabolic
                    y = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
                case 3
                    % cubic
                    y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
                case 4
                    % low-freq sin
                    y = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
                case 5
                    % high-freq sin
                    y = sin(16*pi*x) + noise*(l/num_noise)*randn(M,1);
                case 6
                    % fourth root
                    y = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
                case 7
                    % circle
                    y=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
                case 8
                    % step function
                    y = (x > 0.5) + noise*5*l/num_noise*randn(M,1);
                otherwise
                    error('unknown dep type!');
            end
            % resimulate x so we have null scenario
            x = rand(M,1)*(xMax-xMin)+xMin;
            
            % calculate the metrics
            cimv4Null(ii)     = cim_v4(x,y,minScanIncr);
            cimv4MexNull(ii)  = cim_v4_cc_mex(x,y,minScanIncr);
            cimv8aMexNull(ii) = cim_v8a_cc_mex(x,y,minScanIncr);
            cimv8bMexNull(ii) = cim_v8b_cc_mex(x,y,minScanIncr);
        end
        
        % compute the rejection cutoffs
        cimv4_cut     = quantile(cimv4Null, 0.95);
        cimv4Mex_cut = quantile(cimv4MexNull, 0.95);
        cimv8aMex_cut = quantile(cimv8aMexNull, 0.95);
        cimv8bMex_cut = quantile(cimv8bMexNull, 0.95);
        
        % resimulate the data under the alternative hypothesis
        parfor ii=1:nsim_alt
            x = rand(M,1)*(xMax-xMin)+xMin;
            switch(typ)
                case 1
                    % linear
                    y = x + noise*(l/num_noise)*randn(M,1); 
                case 2
                    % parabolic
                    y = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
                case 3
                    % cubic
                    y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3) + 10*noise*(l/num_noise)*randn(M,1);
                case 4
                    % low-freq sin
                    y = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
                case 5
                    % high-freq sin
                    y = sin(16*pi*x) + noise*(l/num_noise)*randn(M,1);
                case 6
                    % fourth root
                    y = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
                case 7
                    % circle
                    y=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
                case 8
                    % step function
                    y = (x > 0.5) + noise*5*l/num_noise*randn(M,1);
                otherwise
                    error('unknown dep type!');
            end
            
            % calculate the metrics
            cimv4Alt(ii)     = cim_v4(x,y,minScanIncr);
            cimv4MexAlt(ii)  = cim_v4_cc_mex(x,y,minScanIncr);
            cimv8aMexAlt(ii) = cim_v8a_cc_mex(x,y,minScanIncr);
            cimv8bMexAlt(ii) = cim_v8b_cc_mex(x,y,minScanIncr);
        end
        
        % compute the power
        cimv4Power(typ,lIdx)       = sum(cimv4Alt > cimv4_cut)/nsim_alt;
        cimv4MexPower(typ, lIdx)   = sum(cimv4MexAlt > cimv4Mex_cut)/nsim_alt;
        cimv8aMexPower(typ, lIdx)  = sum(cimv8aMexAlt > cimv8aMex_cut)/nsim_alt;
        cimv8bMexPower(typ, lIdx)  = sum(cimv8bMexAlt > cimv8bMex_cut)/nsim_alt;
    end
end

% save the data
if(ispc)
    save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\cim_vstar4_8_power_M_%d.mat', M));
elseif(ismac)
    save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/cim_vstar4_8_power_M_%d.mat', M));
else
    save(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/cim_vstar4_8_power_M_%d.mat', M));
end

%%
clear;
clc;
close all;
dbstop if error;

M = 500;  % which one do we want to plot?

if(ispc)
    load(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\cim_vstar4_power_M_%d.mat', M));
elseif(ismac)
    load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/cim_vstar4_power_M_%d.mat', M));
else
    load(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/cim_vstar4_power_M_%d.mat', M));
end


powerMat = zeros(3,8,length(noiseVec));
powerMat(1,:,:) = cimv4MexPower(:,1:length(noiseVec));
powerMat(2,:,:) = cimv8aMexPower(:,1:length(noiseVec));

% load data for cc-mex after latest run
if(ispc)
    load(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\cim_vstar4_mexonly_power_M_%d.mat', M));
elseif(ismac)
    load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/cim_vstar4_mexonly_power_M_%d.mat', M));
else
    load(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/cim_vstar4_mexonly_power_M_%d.mat', M));
end


powerMat(3,:,:) = cimv8bMexPower(:,1:length(noiseVec));
% noiseVec = (num_noise_test_min:num_noise_test_max)/10;

% labels = {'CIMv4', 'CIMv8a', 'CIMv8b'};
labels = {'CIMv4', 'CIMv4cc', 'CIMv4ccMEX'};
plotStyle = 1;
plotPower(powerMat, M, labels, noiseVec, num_noise_test_min, num_noise_test_max, plotStyle)
