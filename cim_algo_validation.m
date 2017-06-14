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
        tauData = resultsMat(mIdx,lIdx,:);
        
        subplotVec = [p plotIdx];
        subplot(subplotVec);
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
        legendCell{legendCellIdx} = sprintf('\sigma^2=%d',noiseLevel);
        legendCellIdx = legendCellIdx + 1;
    end
    % overlay teh theoretical normal distribution of what we expect to
    % see
    y2 = normpdf(xi,0,(2*(2*M+5))/(9*M*(M-1)));
    plot(xi,y2,'-.');
    legendCell{legendCellIdx+1} = 'theoretical';

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

% WARNING: ENSURE THAT minepy/matlab/ is in the matlab path for MIC to
% work!

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
cimv4Null = zeros(1,nsim_null);
cim_smallmsi_Null = zeros(1,nsim_null);
cimv3_smallmsi_Null = zeros(1,nsim_null);
cimv4_smallmsi_Null = zeros(1,nsim_null);

cimAlt  = zeros(1,nsim_alt);
cimv3Alt = zeros(1,nsim_alt);
cimv4Alt = zeros(1,nsim_alt);
cim_smallmsi_Alt = zeros(1,nsim_alt);
cimv3_smallmsi_Alt = zeros(1,nsim_alt);
cimv4_smallmsi_Alt = zeros(1,nsim_alt);

% Arrays holding the estimated power for each of the "correlation" types, 
% for each data type (linear, parabolic, etc...) with each noise level
cimPower = zeros(numDepTests, num_noise);
cimv3Power = zeros(numDepTests, num_noise);
cimv4Power = zeros(numDepTests, num_noise);
cim_smallmsi_Power = zeros(numDepTests,num_noise);
cimv3_smallmsi_Power = zeros(numDepTests,num_noise);
cimv4_smallmsi_Power = zeros(numDepTests,num_noise);

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
%         for ii=1:nsim_null
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
            cimNull(ii) = cim(x, y);
            cimv3Null(ii) = cim_v3(x, y);
            cimv4Null(ii) = cim_v4(x,y);
            cim_smallmsi_Null(ii) = cim(x,y,minscanincrVal);
            cimv3_smallmsi_Null(ii) = cim_v3(x,y,minscanincrVal);
            cimv4_smallmsi_Null(ii) = cim_v4(x,y,minscanincrVal);
        end
        
        % compute the rejection cutoffs
        cim_cut = quantile(cimNull, 0.95);
        cimv3_cut = quantile(cimv3Null, 0.95);
        cimv4_cut = quantile(cimv4Null, 0.95);
        cim_smallmsi_cut = quantile(cim_smallmsi_Null, 0.95);
        cimv3_smallmsi_cut = quantile(cimv3_smallmsi_Null, 0.95);
        cimv4_smallmsi_cut = quantile(cimv4_smallmsi_Null, 0.95);
        
        % resimulate the data under the alternative hypothesis
        parfor ii=1:nsim_alt
%         for ii=1:nsim_alt
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
            cimAlt(ii) = cim(x, y);
            cimv3Alt(ii) = cim_v3(x, y);
            cimv4Alt(ii) = cim_v4(x,y);
            cim_smallmsi_Alt(ii) = cim(x,y,minscanincrVal);
            cimv3_smallmsi_Alt(ii) = cim_v3(x,y,minscanincrVal);
            cimv4_smallmsi_Alt(ii) = cim_v4(x,y,minscanincrVal);
        end
        
        % compute the power
        cimPower(typ, l)   = sum(cimAlt > cim_cut)/nsim_alt;
        cimv3Power(typ, l)  = sum(cimv3Alt > cimv3_cut)/nsim_alt;
        cimv4Power(typ, l) = sum(cimv4Alt > cimv4_cut)/nsim_alt;
        cim_smallmsi_Power(typ, l)  = sum(cim_smallmsi_Alt > cim_smallmsi_cut)/nsim_alt;
        cimv3_smallmsi_Power(typ, l)  = sum(cimv3_smallmsi_Alt > cimv3_smallmsi_cut)/nsim_alt;
        cimv4_smallmsi_Power(typ, l)  = sum(cimv4_smallmsi_Alt > cimv4_smallmsi_cut)/nsim_alt;
    end
end

% save the data
if(ispc)
    save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\cimv3_power_M_%d.mat', M));
elseif(ismac)
    save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/cimv3_power_M_%d.mat', M));
else
    save(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/cimv3_power_M_%d.mat', M));
end
%%
clear;
clc;
close all;

M = 500;  % which one do we want to plot?

if(ispc)
    load(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\cimv3_power_M_%d.mat', M));
elseif(ismac)
    load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/cimv3_power_M_%d.mat', M));
else
    load(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/cimv3_power_M_%d.mat', M));
end

% inlet plot configuration
M_inlet = 200;
if(M==500)
    inset_bufX = 0.0005; inset_bufY = 0.002;
else
    inset_bufX = 0.15; inset_bufY = 0.26;
end

inset_width = 0.1; inset_height = 0.08;

noiseVec = (num_noise_test_min:num_noise_test_max)/10;
figure;
h1 = subplot(2,2,1);
hh1 = plot(noiseVec, cimPower(1,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, cimv3Power(1,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, cim_smallmsi_Power(1,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, cimv3_smallmsi_Power(1,num_noise_test_min:num_noise_test_max), 'v-.', ...
     noiseVec, cimv4Power(1,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, cimv4_smallmsi_Power(1,num_noise_test_min:num_noise_test_max), 'p-.'); 
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel({'Noise Level','(a)'}, 'FontSize', 20); ylabel('Power', 'FontSize', 20); grid on;
% title('(a)', 'FontSize', 20);
h1.FontSize = 20; 
loc_inset = [h1.Position(1)+inset_bufX h1.Position(2)+inset_bufY inset_width inset_height];
ax1 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = tmp1;
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax1.Box = 'on'; ax1.XTick = []; ax1.YTick = [];
ax1.XLim = [min(tmp1) max(tmp1)];
ax1.YLim = [min(tmp2) max(tmp2)];
hh1(1).LineWidth = 1.5; 
hh1(2).LineWidth = 1.5; 
hh1(3).LineWidth = 1.5; 
hh1(4).LineWidth = 1.5; 
hh1(5).LineWidth = 1.5; 
hh1(6).LineWidth = 5; 

h2 = subplot(2,2,2);
hh2 = plot(noiseVec, cimPower(2,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, cimv3Power(2,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, cim_smallmsi_Power(2,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, cimv3_smallmsi_Power(2,num_noise_test_min:num_noise_test_max), 'v-.', ...
     noiseVec, cimv4Power(2,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, cimv4_smallmsi_Power(2,num_noise_test_min:num_noise_test_max), 'p-.');
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel({'Noise Level', '(b)'}, 'FontSize', 20); ylabel('Power', 'FontSize', 20); grid on;
% title('(b)', 'FontSize', 20);
h2.FontSize = 20; 
loc_inset = [h2.Position(1)+inset_bufX h2.Position(2)+inset_bufY inset_width inset_height];
ax2 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = 4*(tmp1-.5).^2;
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax2.Box = 'on'; ax2.XTick = []; ax2.YTick = [];
ax2.XLim = [min(tmp1) max(tmp1)];
ax2.YLim = [min(tmp2) max(tmp2)];
hh2(1).LineWidth = 1.5; 
hh2(2).LineWidth = 1.5; 
hh2(3).LineWidth = 1.5; 
hh2(4).LineWidth = 1.5; 
hh2(5).LineWidth = 1.5; 
hh2(6).LineWidth = 5.5; 

h3 = subplot(2,2,3); 
hh3 = plot(noiseVec, cimPower(3,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, cimv3Power(3,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, cim_smallmsi_Power(3,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, cimv3_smallmsi_Power(3,num_noise_test_min:num_noise_test_max), 'v-.', ...
     noiseVec, cimv4Power(3,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, cimv4_smallmsi_Power(3,num_noise_test_min:num_noise_test_max), 'p-.');
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel({'Noise Level', '(c)'}, 'FontSize', 20); ylabel('Power', 'FontSize', 20); grid on;
% title('(c)', 'FontSize', 20);
h3.FontSize = 20; 
loc_inset = [h3.Position(1)+inset_bufX h3.Position(2)+inset_bufY inset_width inset_height];
ax3 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = 128*(tmp1-1/3).^3-48*(tmp1-1/3).^3-12*(tmp1-1/3);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax3.Box = 'on'; ax3.XTick = []; ax3.YTick = [];
ax3.XLim = [min(tmp1) max(tmp1)];
ax3.YLim = [min(tmp2) max(tmp2)];
hh3(1).LineWidth = 1.5; 
hh3(2).LineWidth = 1.5; 
hh3(3).LineWidth = 1.5; 
hh3(4).LineWidth = 1.5; 
hh3(5).LineWidth = 1.5; 
hh3(6).LineWidth = 5; 

h4 = subplot(2,2,4); 
hh4 = plot(noiseVec, cimPower(4,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, cimv3Power(4,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, cim_smallmsi_Power(4,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, cimv3_smallmsi_Power(4,num_noise_test_min:num_noise_test_max), 'v-.', ...
     noiseVec, cimv4Power(4,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, cimv4_smallmsi_Power(4,num_noise_test_min:num_noise_test_max), 'p-.');
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel({'Noise Level', 'd'}, 'FontSize', 20); ylabel('Power', 'FontSize', 20); grid on;
% title('(d)', 'FontSize', 20);
h4.FontSize = 20; 
loc_inset = [h4.Position(1)+inset_bufX h4.Position(2)+inset_bufY inset_width inset_height];
ax4 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = sin(4*pi*tmp1);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax4.Box = 'on'; ax4.XTick = []; ax4.YTick = [];
ax4.XLim = [min(tmp1) max(tmp1)];
ax4.YLim = [min(tmp2) max(tmp2)];
hh4(1).LineWidth = 1.5; 
hh4(2).LineWidth = 1.5; 
hh4(3).LineWidth = 1.5; 
hh4(4).LineWidth = 1.5; 
hh4(5).LineWidth = 1.5; 
hh4(6).LineWidth = 5; 

figure;
h5 = subplot(2,2,1); 
hh5 = plot(noiseVec, cimPower(5,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, cimv3Power(5,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, cim_smallmsi_Power(5,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, cimv3_smallmsi_Power(5,num_noise_test_min:num_noise_test_max), 'v-.', ...
     noiseVec, cimv4Power(5,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, cimv4_smallmsi_Power(5,num_noise_test_min:num_noise_test_max), 'p-.');
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel({'Noise Level', '(e)'}, 'FontSize', 20); ylabel('Power', 'FontSize', 20); grid on;
% title('(e)', 'FontSize', 20);
h5.FontSize = 20; 
loc_inset = [h5.Position(1)+inset_bufX h5.Position(2)+inset_bufY inset_width inset_height];
ax5 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = sin(16*pi*tmp1);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax5.Box = 'on'; ax5.XTick = []; ax5.YTick = [];
ax5.XLim = [min(tmp1) max(tmp1)];
ax5.YLim = [min(tmp2) max(tmp2)];
hh5(1).LineWidth = 1.5; 
hh5(2).LineWidth = 1.5; 
hh5(3).LineWidth = 1.5; 
hh5(4).LineWidth = 1.5; 
hh5(5).LineWidth = 1.5; 
hh5(6).LineWidth = 5; 

h6 = subplot(2,2,2); 
hh6 = plot(noiseVec, cimPower(6,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, cimv3Power(6,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, cim_smallmsi_Power(6,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, cimv3_smallmsi_Power(6,num_noise_test_min:num_noise_test_max), 'v-.', ...
     noiseVec, cimv4Power(6,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, cimv4_smallmsi_Power(6,num_noise_test_min:num_noise_test_max), 'p-.');
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel({'Noise Level', '(f)'}, 'FontSize', 20); ylabel('Power', 'FontSize', 20); grid on;
% title('(f)', 'FontSize', 20);
h6.FontSize = 20; 
loc_inset = [h6.Position(1)+inset_bufX h6.Position(2)+inset_bufY inset_width inset_height];
ax6 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = tmp1.^(1/4);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax6.Box = 'on'; ax6.XTick = []; ax6.YTick = [];
ax6.XLim = [min(tmp1) max(tmp1)];
ax6.YLim = [min(tmp2) max(tmp2)];
hh6(1).LineWidth = 1.5; 
hh6(2).LineWidth = 1.5; 
hh6(3).LineWidth = 1.5; 
hh6(4).LineWidth = 1.5; 
hh6(5).LineWidth = 1.5; 
hh6(6).LineWidth = 5; 

h7 = subplot(2,2,3); 
hh7 = plot(noiseVec, cimPower(7,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, cimv3Power(7,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, cim_smallmsi_Power(7,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, cimv3_smallmsi_Power(7,num_noise_test_min:num_noise_test_max), 'v-.', ...
     noiseVec, cimv4Power(7,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, cimv4_smallmsi_Power(7,num_noise_test_min:num_noise_test_max), 'p-.');
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel({'Noise Level', '(g)'}, 'FontSize', 20); ylabel('Power', 'FontSize', 20); grid on;
% title('(g)', 'FontSize', 20);
h7.FontSize = 20; 
loc_inset = [h7.Position(1)+inset_bufX h7.Position(2)+inset_bufY inset_width inset_height];
ax7 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet/2);
tmp2 = (sqrt(1 - (2*tmp1 - 1).^2));
tmp3 = -(sqrt(1 - (2*tmp1 - 1).^2));
plot(tmp1,tmp2, 'k', 'LineWidth', 2); hold on;
plot(tmp1,tmp3, 'k', 'LineWidth', 2); 
ax7.Box = 'on'; ax7.XTick = []; ax7.YTick = [];
ax7.XLim = [min(tmp1) max(tmp1)];
ax7.YLim = [min(tmp3) max(tmp2)];
hh7(1).LineWidth = 1.5; 
hh7(2).LineWidth = 1.5; 
hh7(3).LineWidth = 1.5; 
hh7(4).LineWidth = 1.5; 
hh7(5).LineWidth = 1.5; 
hh7(6).LineWidth = 5; 

h8 = subplot(2,2,4); 
hh8 = plot(noiseVec, cimPower(8,num_noise_test_min:num_noise_test_max), 'o-.', ...
     noiseVec, cimv3Power(8,num_noise_test_min:num_noise_test_max), '+-.', ...
     noiseVec, cim_smallmsi_Power(8,num_noise_test_min:num_noise_test_max), 'd-.', ...
     noiseVec, cimv3_smallmsi_Power(8,num_noise_test_min:num_noise_test_max), 'v-.', ...
     noiseVec, cimv4Power(8,num_noise_test_min:num_noise_test_max), 's-.', ...
     noiseVec, cimv4_smallmsi_Power(8,num_noise_test_min:num_noise_test_max), 'p-.');
axis([min(noiseVec) max(noiseVec) 0 1]);
h8.FontSize = 20; 
legend('CIM', 'CIMv3', 'CIM(sMSI)', 'CIMv3(sMSI)', 'CIMv4', 'CIMv4(sMSI)');
xlabel({'Noise Level', '(h)'}, 'FontSize', 20); ylabel('Power', 'FontSize', 20); grid on;
% title('(h)', 'FontSize', 20);
loc_inset = [h8.Position(1)+inset_bufX h8.Position(2)+inset_bufY inset_width inset_height];
ax8 = axes('Position',loc_inset);
tmp1 = linspace(0,1,M_inlet);
tmp2 = (tmp1 > 0.5);
plot(tmp1,tmp2, 'k', 'LineWidth', 2);
ax8.Box = 'on'; ax8.XTick = []; ax8.YTick = [];
ax8.XLim = [min(tmp1) max(tmp1)];
hh8(1).LineWidth = 1.5; 
hh8(2).LineWidth = 1.5; 
hh8(3).LineWidth = 1.5; 
hh8(4).LineWidth = 1.5; 
hh8(5).LineWidth = 1.5; 
hh8(6).LineWidth = 5; 
