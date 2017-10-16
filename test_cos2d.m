% Test CoS-2D

clear;
clc;
dbstop if error;

n = 500;

x = rand(n,1);
y1 = (x-0.5).^2;
y2 = sin(4*pi*x);
y3 = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3);

nVar = 0.1;

y = y3 + randn(n,1)*nVar;

% rho = .7;
% Z = mvnrnd([0 0], [1 rho; rho 1], n);
% U = normcdf(Z);
% X = [gaminv(U(:,1),2,1) tinv(U(:,2),5)];
% x = X(:,1); y = X(:,2);

u = pobs(x);
v = pobs(y);
W = min([u v],[],2);
M = max(u+v-1,0);

subplot(2,2,1);
scatter(x,y);
grid on;
xlabel('x'); ylabel('y');

subplot(2,2,2);
scatter(u,v)
grid on;
xlabel('u'); ylabel('v');

% [cosValue1,RR] = cos2d_v2(x,y1);
[cosValue2,RR] = cos2d_v2(x,y);

% fprintf('CoS(x,y1)=%0.2f CoS(x,y2)=%0.02f\n',cosValue1,cosValue2);

% plot the empirical copula -- color coded by the domain
colors = {'r', 'b'};
numDomains = unique(RR(:,end));
numDomains = numDomains';

subplot(2,2,3);
hold on;
for domain=numDomains
    if(mod(domain,2)==1)
        c = 'r';
    else
        c = 'b';
    end
    II = find(RR(:,end)==domain);
    uu = RR(II,2);
    vals = RR(II,end-1);
    plot(uu,vals,c);
end
grid on;
%scatter(u,M); 
%scatter(u,W);
xlabel('u'); ylabel('C(u,v)');
title(sprintf('CoS=%0.02f',cosValue2));
%legend('C','M','W')

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

%%

clear;
clc;

numMCSim = 500;

M = 2000;

R = [1 0.1; 0.1 1];
cosValVec = zeros(1,numMCSim);
for mcSim=1:numMCSim
    X = mvnrnd([0 0],R,M);
    cosVal = cos2d_v2(X(:,1),X(:,2));
    cosValVec(mcSim) = cosVal;    
end

fprintf('rho=0.1 E[CoS]=%0.02f sigma[CoS]=%0.02f\n', mean(cosValVec), std(cosValVec));

R = [1 0.3; 0.3 1];
cosValVec = zeros(1,numMCSim);
for mcSim=1:numMCSim
    X = mvnrnd([0 0],R,M);
    cosVal = cos2d_v2(X(:,1),X(:,2));
    cosValVec(mcSim) = cosVal;
end

fprintf('rho=0.3 E[CoS]=%0.02f sigma[CoS]=%0.02f\n', mean(cosValVec), std(cosValVec));

R = [1 0.5; 0.5 1];
cosValVec = zeros(1,numMCSim);
for mcSim=1:numMCSim
    X = mvnrnd([0 0],R,M);
    cosVal = cos2d_v2(X(:,1),X(:,2));
    cosValVec(mcSim) = cosVal;
end
fprintf('rho=0.5 E[CoS]=%0.02f sigma[CoS]=%0.02f\n', mean(cosValVec), std(cosValVec));

R = [1 0.7; 0.7 1];
cosValVec = zeros(1,numMCSim);
for mcSim=1:numMCSim
    X = mvnrnd([0 0],R,M);
    cosVal = cos2d_v2(X(:,1),X(:,2));
    cosValVec(mcSim) = cosVal;
end

fprintf('rho=0.7 E[CoS]=%0.02f sigma[CoS]=%0.02f\n', mean(cosValVec), std(cosValVec));


R = [1 0.9; 0.9 1];
cosValVec = zeros(1,numMCSim);
for mcSim=1:numMCSim
    X = mvnrnd([0 0],R,M);
    cosVal = cos2d_v2(X(:,1),X(:,2));
    cosValVec(mcSim) = cosVal;
end

fprintf('rho=0.9 E[CoS]=%0.02f sigma[CoS]=%0.02f\n', mean(cosValVec), std(cosValVec));
