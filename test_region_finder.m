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

%% check kendall's tau's distribution for functionally dependent but P(conc)=P(disc) relationships
clear;
clc;
dbstop if error;

numMCSims = 100;
MVec = 50:50:500;
xMin = 0;
xMax = 1;
num_noise = 30;                    % The number of different noise levels used
noise = 3; 
num_noise_test_min = 0;
num_noise_test_max = 30;
noiseVec = num_noise_test_min:num_noise_test_max;

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
