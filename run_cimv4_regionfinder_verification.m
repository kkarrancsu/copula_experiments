%% verify the CIM-V4 region finder
clear;
clc;
dbstop if error;

M = 500;
x = rand(M,1);

noise = 3; l = 1; num_noise = 30;

y1 = x + noise*(l/num_noise)*randn(M,1); 
y2 = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
y3 = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
y4 = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
y5 = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
y6 = (2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
% [z1,z2] = cim_v4(x,y1)

u = copularnd('gaussian',0.9,500);
[z1,z2] = cim_v4(u(:,1),u(:,2))
