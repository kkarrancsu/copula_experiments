%% Test the MEX implementation of CIM_V4
clear;
clc;

numMCSim = 25;

num_noise = 30;                    % The number of different noise levels used
noise = 3;                         % A constant to determine the amount of noise
l = 10;

sseVec = zeros(8,numMCSim);

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');
M = 500;

for mcSimNum=1:numMCSim
    dispstat(sprintf('%d/%d',mcSimNum, numMCSim),'keepthis', 'timestamp');
    x = rand(M,1); 
    
    y = x+ noise*(l/num_noise)*randn(M,1); 
    c1 = cim_v4(x,y); 
    c2 = cim_v4_cc_mex(x,y);
    diffVal = (c1-c2)^2;
    sseVec(1,mcSimNum) = diffVal;
    
    y = 4*(x-0.5).^2+ noise*(l/num_noise)*randn(M,1);
    c1 = cim_v4(x,y); 
    c2 = cim_v4_cc_mex(x,y);
    diffVal = (c1-c2)^2;
    sseVec(2,mcSimNum) = diffVal;
    
    y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
    c1 = cim_v4(x,y); 
    c2 = cim_v4_cc_mex(x,y);
    diffVal = (c1-c2)^2;
    sseVec(3,mcSimNum) = diffVal;
    
    y = sin(4*pi*x)+ 2*noise*(l/num_noise)*randn(M,1);
    c1 = cim_v4(x,y); 
    c2 = cim_v4_cc_mex(x,y);
    diffVal = (c1-c2)^2;
    sseVec(4,mcSimNum) = diffVal;
    
    y = sin(16*pi*x)+ noise*(l/num_noise)*randn(M,1);
    c1 = cim_v4(x,y); 
    c2 = cim_v4_cc_mex(x,y);
    diffVal = (c1-c2)^2;
    sseVec(5,mcSimNum) = diffVal;
    
    y = x.^(1/4)+ noise*(l/num_noise)*randn(M,1);
    c1 = cim_v4(x,y); 
    c2 = cim_v4_cc_mex(x,y);
    diffVal = (c1-c2)^2;
    sseVec(6,mcSimNum) = diffVal;
    
    y=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2))+ noise/4*l/num_noise*randn(M,1);
    c1 = cim_v4(x,y); 
    c2 = cim_v4_cc_mex(x,y);
    diffVal = (c1-c2)^2;
    sseVec(7,mcSimNum) = diffVal;
    
    y = double((x > 0.5))+ noise*5*l/num_noise*randn(M,1);
    c1 = cim_v4(x,y); 
    c2 = cim_v4_cc_mex(x,y);
    diffVal = (c1-c2)^2;
    sseVec(8,mcSimNum) = diffVal;
    
end

sum(sseVec,2)