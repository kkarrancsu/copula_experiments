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
            y1 = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
            % compute tau so we can look at the distribution of it, to see if
            % it matches the N(0,2(2n+5)/(9n(n+1))) as proved by Kendall
            tauVal = corr(x,y1,'type','kendall');
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
    y1 = normpdf(xi,0,(2*(2*M+5))/(9*M*(M-1)));
    plot(xi,y1,'-.');
    legendCell{legendCellIdx+1} = 'theoretical';

    title(sprintf('M=%d',M));
    grid on; xlabel('\tau');
    legend(legendCell);
    
end

%% test the CIM region finder
clear;
clc;

numMCSims = 100;
% MVec = 100:100:2000;
MVec = [100, 500, 2000];
xMin = 0;
xMax = 1;
num_noise = 30;                    % The number of different noise levels used
noise = 3; 
num_noise_test_min = 0;
num_noise_test_max = 30;
% noiseVec = num_noise_test_min:num_noise_test_max;
noiseVec = [0,10,20];
minscanincr = 0.025;

resultsMat_na_parabola = cell(length(MVec),length(noiseVec),numMCSims);
resultsMat_a_parabola  = cell(length(MVec),length(noiseVec),numMCSims);
resultsMat_na_sinu     = cell(length(MVec),length(noiseVec),numMCSims);
resultsMat_a_sinu      = cell(length(MVec),length(noiseVec),numMCSims);
resultsMat_na_cubic    = cell(length(MVec),length(noiseVec),numMCSims);
resultsMat_a_cubic     = cell(length(MVec),length(noiseVec),numMCSims);

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
            y1 = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
            y2 = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
            y3 = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
            
            signature_na_parabola = cim_region_finder(x,y1,minscanincr,0);
            signature_a_parabola = cim_region_finder(x,y1,minscanincr,1);
            signature_na_sinu = cim_region_finder(x,y2,minscanincr,0);
            signature_a_sinu = cim_region_finder(x,y2,minscanincr,1);
            signature_na_cubic = cim_region_finder(x,y3,minscanincr,0);
            signature_a_cubic = cim_region_finder(x,y3,minscanincr,1);
            
            resultsMat_na_parabola{mIdx,noiseIdx,ii} = signature_na_parabola;
            resultsMat_a_parabola{mIdx,noiseIdx,ii} = signature_a_parabola;
            resultsMat_na_sinu{mIdx,noiseIdx,ii} = signature_na_sinu;
            resultsMat_a_sinu{mIdx,noiseIdx,ii} = signature_a_sinu;
            resultsMat_na_cubic{mIdx,noiseIdx,ii} = signature_na_cubic;
            resultsMat_a_cubic{mIdx,noiseIdx,ii} = signature_a_cubic;
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
MVecToPlot = [100, 500, 2000];
noiseLevelsToPlot = [0,10,20];
scanincrsToPlot = [0.5, 0.25, 0.125,.0625];
subplotCfg = [length(MVecToPlot), length(noiseLevelsToPlot)];
numPlots = prod(subplotCfg);
lineMarkers = {'+-.','o-.','*-.','d-.','x-.','s.-'};

f1 = figure(1);
f2 = figure(2);
f3 = figure(3);

plotIdx = 1;
for MIdx=1:length(MVecToPlot)
    % find which index this goes into the resultsMat
    M = MVecToPlot(MIdx);
    mIdx = find(M==MVec);
    
    for noiseLevelIdx=1:length(noiseLevelsToPlot)
        legendCell1 = {}; legendCell1Idx = 1;
        legendCell2 = {}; legendCell2Idx = 1;
        legendCell3 = {}; legendCell3Idx = 1;

        noiseLevel = noiseLevelsToPlot(noiseLevelIdx);
        lIdx = find(noiseLevel==noiseVec);
        % get the data
        [avgMetric_na_parabola, avgNumPts_na_parabola] = ...
            getSignatureAvg(resultsMat_na_parabola,mIdx,lIdx,numMCSims,scanincrsToPlot);
        [avgMetric_a_parabola, avgNumPts_a_parabola] = ...
            getSignatureAvg(resultsMat_a_parabola,mIdx,lIdx,numMCSims,scanincrsToPlot);
        [avgMetric_na_sinu, avgNumPts_na_sinu] = ...
            getSignatureAvg(resultsMat_na_sinu,mIdx,lIdx,numMCSims,scanincrsToPlot);
        [avgMetric_a_sinu, avgNumPts_a_sinu] = ...
            getSignatureAvg(resultsMat_a_sinu,mIdx,lIdx,numMCSims,scanincrsToPlot);
        [avgMetric_na_cubic, avgNumPts_na_cubic] = ...
            getSignatureAvg(resultsMat_na_cubic,mIdx,lIdx,numMCSims,scanincrsToPlot);
        [avgMetric_a_cubic, avgNumPts_a_cubic] = ...
            getSignatureAvg(resultsMat_a_cubic,mIdx,lIdx,numMCSims,scanincrsToPlot);
        
        set(0,'CurrentFigure',f1);
        subplotCfgVec = [subplotCfg plotIdx]; 
        subplotCfgNum = subplotCfgVec(1)*100+subplotCfgVec(2)*10+subplotCfgVec(3);
        subplot(subplotCfgNum);
        % compute and plot the z-score from the metric & numPts
        for ii=1:length(scanincrsToPlot)
            scanincr = scanincrsToPlot(ii);
            n = avgNumPts_na_parabola{ii};
            den = (2*(2*n+5))./(9*n.*(n-1));
            x = scanincr:scanincr:1;
            y = avgMetric_na_parabola{ii} ./ den;
            plot(x,y,lineMarkers{ii},'LineWidth',3);
            hold on; grid on;
%             n = avgNumPts_a_parabola{ii};
%             den = (2*(2*n+5))./(9*n.*(n-1));
%             plot(avgMetric_a_parabola{ii} ./ den );
            legendCell1{legendCell1Idx} = sprintf('NA - \\Delta i=%0.03f', scanincrsToPlot(ii));
            legendCell1Idx = legendCell1Idx + 1;
%             legendCell1{legendCell1Idx} = sprintf('A - \\Delta i=%0.02f, \\sigma^2=%d',...
%                 scanincrsToPlot(ii),noiseLevel);
%             legendCell1Idx = legendCell1Idx + 1;
        end
        title(sprintf('M=%d \\sigma^2=%d',M,noiseLevel));
        if(MIdx==1 && noiseLevelIdx==1)
            legend(legendCell1,'Location','SouthWest');
        end
        
        set(0,'CurrentFigure',f2);
        subplot(subplotCfgNum);
        for ii=1:length(scanincrsToPlot)
            scanincr = scanincrsToPlot(ii);
            n = avgNumPts_na_sinu{ii};
            den = (2*(2*n+5))./(9*n.*(n-1));
            x = scanincr:scanincr:1;
            y = avgMetric_na_sinu{ii} ./ den;
            plot(x,y,lineMarkers{ii},'LineWidth',3);
            hold on; grid on;
%             n = avgNumPts_a_sinu{ii};
%             den = (2*(2*n+5))./(9*n.*(n-1));
%             plot(avgMetric_a_sinu{ii} ./ den );
        end
        title(sprintf('M=%d \\sigma^2=%d',M,noiseLevel));
        if(MIdx==1 && noiseLevelIdx==1)
            legend(legendCell1,'Location','SouthWest');
        end
        
        set(0,'CurrentFigure',f3);
        subplot(subplotCfgNum);
        for ii=1:length(scanincrsToPlot)
            scanincr = scanincrsToPlot(ii);
            n = avgNumPts_na_cubic{ii};
            den = (2*(2*n+5))./(9*n.*(n-1));
            x = scanincr:scanincr:1;
            y = avgMetric_na_cubic{ii} ./ den;
            plot(x,y,lineMarkers{ii},'LineWidth',3);
            hold on; grid on;
%             n = avgNumPts_a_cubic{ii};
%             den = (2*(2*n+5))./(9*n.*(n-1));
%             plot(avgMetric_a_cubic{ii} ./ den );
        end
        title(sprintf('M=%d \\sigma^2=%d',M,noiseLevel));
        if(MIdx==1 && noiseLevelIdx==1)
            legend(legendCell1,'Location','SouthWest');
        end
        
        plotIdx = plotIdx + 1;
    end
    
    set(0,'CurrentFigure',f1); figtitle('Parabolic');
    set(0,'CurrentFigure',f2); figtitle('Sinusoidal');
    set(0,'CurrentFigure',f3); figtitle('Cubic');
end