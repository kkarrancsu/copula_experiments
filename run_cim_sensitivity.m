%% run script for testing the sensitivity of the CIM algorithm

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

M = 200;
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
num_noise_test_max = 30;
noiseVec = num_noise_test_min:num_noise_test_max;
plotPower(powerCurve, M, labels, noiseVec, num_noise_test_min, num_noise_test_max, plotStyle)

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

M = 200;
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

num_noise_test_min = 0;
num_noise_test_max = 30;
noiseVec = num_noise_test_min:num_noise_test_max;


plotAlgoSensitivity(algoSensitivityData, scanincrsToTest, noiseVec)