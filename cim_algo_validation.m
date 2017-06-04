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
            y = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
            % compute tau so we can look at the distribution of it, to see if
            % it matches the N(0,2(2n+5)/(9n(n+1))) as proved by Kendall
            tauVal = corr(x,y,'type','kendall');
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
    end
    % overlay teh theoretical normal distribution of what we expect to
    % see
    y = normpdf(xi,0,(2*(2*M+5))/(9*M*(M-1)));
    plot(xi,y,'-.');
    legendCell{legendCellIdx+1} = 'theoretical';

    title(sprintf('M=%d',M));
    grid on; xlabel('\tau');
    legend(legendCell);
    
end

%% test the CIM region finder