% Test CoS-2D

clear;
clc;
dbstop if error;

M = 500;

x = rand(M,1);
num_noise = 3;
noise = 3;                         % A constant to determine the amount of noise
l = 0.50;

y1 = x + noise*(l/num_noise)*randn(M,1); 
y2 = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
y3 = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
y4 = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
y5 = sin(16*pi*x) + noise*(l/num_noise)*randn(M,1);
y6 = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
y7 =(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
y8 = (x > 0.5) + noise*5*l/num_noise*randn(M,1);

y = y3;

% rho = .7;
% Z = mvnrnd([0 0], [1 rho; rho 1], n);
% U = normcdf(Z);
% X = [gaminv(U(:,1),2,1) tinv(U(:,2),5)];
% x = X(:,1); y = X(:,2);

u = pobs(x);
v = pobs(y);

subplot(2,2,2);
scatter(x,y);
grid on;
xlabel('x'); ylabel('y');

subplot(2,2,1);
scatter(u,v)
grid on;
xlabel('u'); ylabel('v');

[cosValue,RR] = cos2d_v3(x,y);

W = max(0,1-2+RR(:,2)+RR(:,3));
M = min(RR(:,2),RR(:,3));

% fprintf('CoS(x,y1)=%0.2f CoS(x,y2)=%0.02f\n',cosValue1,cosValue2);

% plot the empirical copula -- color coded by the domain
colors = {'r', 'b'};
numDomains = unique(RR(:,end));
numDomains = numDomains';

subplot(2,2,3);
hold on;
for domain=numDomains
    if(mod(domain,2)==1)
        c = 'r+-.';
    else
        c = 'b+-.';
    end
    II = find(RR(:,end)==domain);
    uu = RR(II,2);
    vals = RR(II,end-2);
    plot(uu,vals,c);
end
% plot(RR(:,2),RR(:,4),'o-.'); hold on;
% plot(RR(:,2),M,'+-.'); 
% plot(RR(:,2),W,'d-.');
plot(RR(:,2),RR(:,5),'gv-.'); hold on;
xlabel('u'); ylabel('C(u,v)');
title(sprintf('CoS=%0.02f',cosValue));
% legend('C','M','W','Pi','location','northwest')
% legend('C','Pi','location','northwest')
grid on;

subplot(2,2,4);
K = 25;
U = [u v]; 
C_est = empcopulacdf(U, K, 'deheuvels');
u = linspace(0.01,0.99,K);
[U1,U2] = ndgrid(u);
mesh(U1,U2,C_est);
xlabel('u')
ylabel('v')
title('C')
rotate3d on

% figure;
% subplot(2,2,1);
% plot(x);
% xlabel('t'); ylabel('x');
% subplot(2,2,2);
% plot(y);
% xlabel('t'); ylabel('y');
% subplot(2,2,3);
% pwelch(x);
% 
% subplot(2,2,4);
% pwelch(y);

%%

clear;
clc;

numMCSim = 500;

M = 2000;

R = [1 0.1; 0.1 1];
cosValVec = zeros(1,numMCSim);
for mcSim=1:numMCSim
    X = mvnrnd([0 0],R,M);
    cosVal = cos2d_v3(X(:,1),X(:,2));
    cosValVec(mcSim) = cosVal;    
end

fprintf('rho=0.1 E[CoS]=%0.02f sigma[CoS]=%0.02f\n', mean(cosValVec), std(cosValVec));

R = [1 0.3; 0.3 1];
cosValVec = zeros(1,numMCSim);
for mcSim=1:numMCSim
    X = mvnrnd([0 0],R,M);
    cosVal = cos2d_v3(X(:,1),X(:,2));
    cosValVec(mcSim) = cosVal;
end

fprintf('rho=0.3 E[CoS]=%0.02f sigma[CoS]=%0.02f\n', mean(cosValVec), std(cosValVec));

R = [1 0.5; 0.5 1];
cosValVec = zeros(1,numMCSim);
for mcSim=1:numMCSim
    X = mvnrnd([0 0],R,M);
    cosVal = cos2d_v3(X(:,1),X(:,2));
    cosValVec(mcSim) = cosVal;
end
fprintf('rho=0.5 E[CoS]=%0.02f sigma[CoS]=%0.02f\n', mean(cosValVec), std(cosValVec));

R = [1 0.7; 0.7 1];
cosValVec = zeros(1,numMCSim);
for mcSim=1:numMCSim
    X = mvnrnd([0 0],R,M);
    cosVal = cos2d_v3(X(:,1),X(:,2));
    cosValVec(mcSim) = cosVal;
end

fprintf('rho=0.7 E[CoS]=%0.02f sigma[CoS]=%0.02f\n', mean(cosValVec), std(cosValVec));


R = [1 0.9; 0.9 1];
cosValVec = zeros(1,numMCSim);
for mcSim=1:numMCSim
    X = mvnrnd([0 0],R,M);
    cosVal = cos2d_v3(X(:,1),X(:,2));
    cosValVec(mcSim) = cosVal;
end

fprintf('rho=0.9 E[CoS]=%0.02f sigma[CoS]=%0.02f\n', mean(cosValVec), std(cosValVec));
