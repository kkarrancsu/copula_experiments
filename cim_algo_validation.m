%% Some experiments to validate the CIM algorithm

%% check kendall's tau's distribution for functionally dependent but P(conc)=P(disc) relationships
clear;
clc;

numMCSims = 1000;
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
        dispstat(sprintf('Computing for M=%d noise=%d', M, typ),'keepthis', 'timestamp');
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

numMCSims = 500;
MVec = [50, 100, 500, 5000];
xMin = 0;
xMax = 1;
num_noise = 30;                    % The number of different noise levels used
noise = 3; 
num_noise_test_min = 0;
num_noise_test_max = 30;
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
MVecToPlot = [50, 100, 500, 5000];
noiseLevelsToPlot = [0,10,20];
% scanincrsToPlot = [0.5, 0.25, 0.125, .0625, .03125];
scanincrsToPlot = [0.5, .25, 0.125];
subplotCfg = [length(MVecToPlot), length(noiseLevelsToPlot)];
numPlots = prod(subplotCfg);
lineMarkers = {'+-.','o-.','*-.','d-.','x-.','s.-'};

numDeps = 7;
figureVec = zeros(1,numDeps);
for ii=1:numDeps
    figureVec(ii) = figure(ii);
end

plotIdx = 1;
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
            getSignatureAvg(resultsMat_linear,mIdx,lIdx,numMCSims,scanincrsToPlot);
        [avgMetric_parabola, avgNumPts_parabola] = ...
            getSignatureAvg(resultsMat_parabola,mIdx,lIdx,numMCSims,scanincrsToPlot);
        [avgMetric_cubic, avgNumPts_cubic] = ...
            getSignatureAvg(resultsMat_cubic,mIdx,lIdx,numMCSims,scanincrsToPlot);
        [avgMetric_sinu, avgNumPts_sinu] = ...
            getSignatureAvg(resultsMat_sinu,mIdx,lIdx,numMCSims,scanincrsToPlot);
        [avgMetric_hfsinu, avgNumPts_hfsinu] = ...
            getSignatureAvg(resultsMat_hfsinu,mIdx,lIdx,numMCSims,scanincrsToPlot);
        [avgMetric_fr, avgNumPts_fr] = ...
            getSignatureAvg(resultsMat_fr,mIdx,lIdx,numMCSims,scanincrsToPlot);
        [avgMetric_step, avgNumPts_step] = ...
            getSignatureAvg(resultsMat_step,mIdx,lIdx,numMCSims,scanincrsToPlot);
        
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