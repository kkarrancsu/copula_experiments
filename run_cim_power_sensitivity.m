%% run script for testing the sensitivity of the CIM algorithm

clear;
clc;
close all;

cimVersion = 6;
scanincrsToTest = [0.25, 0.125, .0625, .03125, .015625];
M = 500;
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

powerCurve = cim_power_sensitivity(cimfunc,M,scanincrsToTest);

% save the results
if(ispc)
    save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\%s_powerSensitivity_M_%d.mat', fnameStr, M));
elseif(ismac)
    save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/%s_powerSensitivity_M_%d.mat', fnameStr, M));
else
    save(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/%s_powerSensitivity_M_%d.mat', fnameStr, M));
end

%% plot the power sensitivity
clear;
clc;
close all;

M = 500;
cimVersion = 6;

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

