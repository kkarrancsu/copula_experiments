%% CIM power debug
clear;
clc;
close all;
dbstop if error;

% setup and start the parallel pool
myCluster = parcluster('local');
myCluster.NumWorkers = 6;
saveProfile(myCluster);
p = gcp

rng(1230);

M = 500;
minScanIncr = 0.015625;

nameIdxCorrelationCell = {'CIMv4', 'CIMv4_cc','CIMv8a_cc','CIMv8b_cc'};

functionHandlesCell = {@cim_v4;
                       @cim_v4_cc_mex;
                       @cim_v8a_cc_mex;
                       @cim_v8b_cc_mex;};
functionArgsCell    = {{minScanIncr};
                       {minScanIncr};
                       {minScanIncr};
                       {minScanIncr};};
[powerCurve] = compute_power_curves(M,functionHandlesCell, functionArgsCell, 200, 200);

% save the data
if(ispc)
    save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\cim_vstar_power_M_%d.mat',M));
elseif(ismac)
    save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/cim_vstar_power_M_%d.mat',M));
else
    save(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/cim_vstar_power_M_%d.mat',M));
end
