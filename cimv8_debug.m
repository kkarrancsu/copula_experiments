%% Generate statistical power curves for CIMv4 vs CIMv8a/b
% same methodology as Simon & Tibshirani:
% http://statweb.stanford.edu/~tibs/reshef/script.R

clear;
clc;
dbstop if error;

% setup and start the parallel pool
myCluster = parcluster('local');
myCluster.NumWorkers = 6;
saveProfile(myCluster);
p = gcp

rng(1234);

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
% cimv4Null     = zeros(1,nsim_null);
cimv8aNull = zeros(1,nsim_null);
cimv8aRev3CCNull = zeros(1,nsim_null);
cimv8aRev3aCCNull = zeros(1,nsim_null);
cimv8aRev3bCCNull = zeros(1,nsim_null);

% cimv4Alt     = zeros(1,nsim_alt);
cimv8aAlt = zeros(1,nsim_alt);
cimv8aRev3CCAlt = zeros(1,nsim_alt);
cimv8aRev3aCCAlt = zeros(1,nsim_alt);
cimv8aRev3bCCAlt = zeros(1,nsim_alt);

% configuration parameter for cim algorithm
minScanIncr=0.015625;

% Arrays holding the estimated power for each of the "correlation" types, 
% for each data type (linear, parabolic, etc...) with each noise level
% cimv4Power     = zeros(numDepTests,num_noise);
cimv8aPower     = zeros(numDepTests,num_noise);
cimv8aRev3CCPower = zeros(numDepTests,num_noise);
cimv8aRev3aCCPower = zeros(numDepTests,num_noise);
cimv8aRev3bCCPower = zeros(numDepTests,num_noise);

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
%             cimv4Null(ii)        = cim_v4(x,y,minScanIncr);
            cimv8aNull(ii)        = cim_v8a(x,y,minScanIncr);
            cimv8aRev3CCNull(ii)  = cim_v8a_rev3cc(x,y,minScanIncr);
            cimv8aRev3aCCNull(ii) = cim_v8a_rev3acc(x,y,minScanIncr);
            cimv8aRev3bCCNull(ii) = cim_v8a_rev3bcc(x,y,minScanIncr);
        end
        
        % compute the rejection cutoffs
%         cimv4_cut        = quantile(cimv4Null, 0.95);
        cimv8a_cut       = quantile(cimv8aNull, 0.95);
        cimv8aRev3CC_cut = quantile(cimv8aRev3CCNull, 0.95);
        cimv8aRev3aCC_cut = quantile(cimv8aRev3aCCNull, 0.95);
        cimv8aRev3bCC_cut = quantile(cimv8aRev3bCCNull, 0.95);
        
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
%             cimv4Alt(ii)        = cim_v4(x,y,minScanIncr);
            cimv8aAlt(ii)       = cim_v8a(x,y,minScanIncr);
            cimv8aRev3CCAlt(ii) = cim_v8a_rev3cc(x,y,minScanIncr);
            cimv8aRev3aCCAlt(ii) = cim_v8a_rev3acc(x,y,minScanIncr);
            cimv8aRev3bCCAlt(ii) = cim_v8a_rev3bcc(x,y,minScanIncr);
        end
        
        % compute the power
%         cimv4Power(typ,lIdx)          = sum(cimv4Alt > cimv4_cut)/nsim_alt;
        cimv8aPower(typ, lIdx)        = sum(cimv8aAlt > cimv8a_cut)/nsim_alt;
        cimv8aRev3CCPower(typ, lIdx)  = sum(cimv8aRev3CCAlt > cimv8aRev3CC_cut)/nsim_alt;
        cimv8aRev3aCCPower(typ, lIdx)  = sum(cimv8aRev3aCCAlt > cimv8aRev3aCC_cut)/nsim_alt;
        cimv8aRev3bCCPower(typ, lIdx)  = sum(cimv8aRev3bCCAlt > cimv8aRev3bCC_cut)/nsim_alt;
    end
    
    % save the data in between data points so we can chart progress easily
    if(ispc)
        save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\cim_v8rev3_debug_power_M_%d.mat', M));
    elseif(ismac)
        save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/cim_v8rev3_debug_power_M_%d.mat', M));
    else
        save(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/cim_v8rev3_debug_power_M_%d.mat', M));
    end
end

% save the data
if(ispc)
    save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\cim_v8rev3_debug_power_M_%d.mat', M));
elseif(ismac)
    save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/cim_v8rev3_debug_power_M_%d.mat', M));
else
    save(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/cim_v8rev3_debug_power_M_%d.mat', M));
end
