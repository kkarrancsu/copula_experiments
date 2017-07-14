%% debug cim_v8a_rev3cc

clear;
clc;
close all;

rng(1234);

M = 500;
minscanincr = 0.015625;
numMCSim = 1000;

noise = 3; l = 10; num_noise = 30;

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');
for ii=1:numMCSim
    dispstat(sprintf('%d/%d',ii, numMCSim),'timestamp');
    
    x = rand(M,1);
    y = x + noise*(l/num_noise)*randn(M,1);
    m1 = cim_v8a(x,y,minscanincr);
    m2 = cim_v8a_rev4cc_mex(x,y,minscanincr);
%     m3 = cim_v4(x,y,minscanincr);
    
%     if(abs(m3-m1)>1*eps)
%         fprintf('Mismatch with v4 and v8\n');
%         pause;
%     end
    
    if(abs(m2-m1)>1*eps)
        fprintf('Mismatch with v8 and v8_rev4\n');
        pause;
    end
end

%% Dig into the mismatched data to see what is up :D
clear;
clc;
close all;

if(ispc)
    load('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\cim_mismatch\\mismatch_v4_v8.mat');
elseif(ismac)
    load('/Users/Kiran/ownCloud/PhD/sim_results/cim_mismatch/mismatch_v4_v8.mat');
else
    load('/home/kiran/ownCloud/PhD/sim_results/cim_mismatch/mismatch_v4_v8.mat');
end

minscanincr = 0.25;

m1 = cim_v4(x,y,minscanincr)
m2 = cim_v8a(x,y,minscanincr)