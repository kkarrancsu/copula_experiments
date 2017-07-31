%% Characterize the region detection of CIM_V8

clear;
clc;
close all;
dbstop if error;

numMCSims = 500;
MVec = [50 100 200:100:1000];
noiseVec = 0:30;
% simulate different change rates
cornerPts = 0.1:0.1:0.9;
minScanIncr = 0.015625;

rdData_up = zeros(length(MVec),length(noiseVec),length(cornerPts),numMCSims);
rdData_down = zeros(length(MVec),length(noiseVec),length(cornerPts),numMCSims);

num_noise = 30;                    % The number of different noise levels used
noise = 3;                         % A constant to determine the amount of noise

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');
for mIdx=1:length(MVec)
    M = MVec(mIdx);
    for noiseIdx=1:length(noiseVec)
        l = noiseVec(noiseIdx);
        dispstat(sprintf('Computing for M=%d noise=%d', M, l),'keepthis', 'timestamp');
        for cornerPtIdx=1:length(cornerPts)
            cornerPt = cornerPts(cornerPtIdx);
            parfor ii=1:numMCSims
                x = rand(M,1);
                y1 = 4*(x-cornerPt).^2 + noise*(l/num_noise)*randn(M,1);
                y2 = -4*(x-cornerPt).^2 + noise*(l/num_noise)*randn(M,1);
                
                [~,rect1] = cim_v8_cc_mex(x,y1,minScanIncr);
                [~,rect2] = cim_v8_cc_mex(x,y1,minScanIncr);
                
                % parse where the region detection found the division
                % marker for the down parabola
                
                if(rect1(4,1)==1)
                    % assume we are u/v orientation
                    r1 = rect1(2,1);
                else
                    % assume we are v/u orientation
                    r1 = rect1(4,1);
                end
                
                % parse where the region detection found the division
                % marker for the up parabola
                if(rect2(4,1)==1)
                    % assume we are u/v orientation
                    r2 = rect2(2,1);
                else
                    % assume we are v/u orientation
                    r2 = rect2(4,1);
                end
                fprintf('r1=%0.02f r2=%0.02f\n', r1, r2);
                rdData_up(mIdx,noiseIdx,cornerPtIdx,ii) = r1;
                rdData_down(mIdx,noiseIdx,cornerPtIdx,ii) = r2;
            end
        end
    end
    % store the results
    if(ispc)
        save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\cim_v8_regionDetection_%d.mat',M));
    elseif(ismac)
        save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/cim_v8_regionDetection_%d.mat',M));
    else
        save(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/cim_v8_regionDetection_%d.mat',M));
    end
end

