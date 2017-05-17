%% configuration
clear;
clc;

xMin = 0;
xMax = 1;
num_noise = 30;                    % The number of different noise levels used
noise = 3;                         % A constant to determine the amount of noise

MVec = 100:100:1000;

options=optimset('Display', 'off');

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');

K = 50;
weightVec = ones(K,1);
numMCSim = 100;

for M=MVec
    % Setup M & W copulas
    x = rand(M,1)*(xMax-xMin)+xMin;
    y1 = x;
    y2 = -x;
    
    M_copula = empcopulacdf(pobs([x,y1]),K,'deheuvels');
    W_copula = empcopulacdf(pobs([x,y2]),K,'deheuvels');

    % Test Empirical copula distance as a distance measure under various
    % amounts of noise
    num_noise_test_min = 1;
    num_noise_test_max = 30;

    noiseVec = num_noise_test_min:num_noise_test_max;
    distancesTensor = zeros(5,length(noiseVec));
    distancesVec_parabola = zeros(5,length(noiseVec));
    distancesVec_cubic = zeros(5,length(noiseVec));
    distancesVec_sinusoid = zeros(5,length(noiseVec));
                    % 1 - SSE to M copula
                    % 2 - SSE to W copula
                    % 3 - eigenvalue metric
                    % 4 - EMD to M copula
                    % 5 - EMD to W copula
    for l=noiseVec
        dispstat(sprintf('Computing for noise level=%d >> M=%d',l,M),'keepthis', 'timestamp');
        M_dist_linear = 0; M_dist_quadratic = 0; M_dist_cubic = 0; M_dist_sinusoid = 0;
        W_dist_linear = 0; W_dist_quadratic = 0; W_dist_cubic = 0; W_dist_sinusoid = 0;
        eig_linear = 0; eig_quadratic = 0; eig_cubic = 0; eig_sinusoid = 0;
        emd_M_linear = 0; emd_M_quadratic = 0; emd_M_cubic = 0; emd_M_sinusoid = 0;
        emd_W_linear = 0; emd_W_quadratic = 0; emd_W_cubic = 0; emd_W_sinusoid = 0;
        parfor mcSimNum=1:numMCSim
            x = rand(M,1)*(xMax-xMin)+xMin;
            y1 = x + noise*(l/num_noise)*randn(M,1);
            y2 = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
            y3 = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
            y4 = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);

            U1 = pobs([x y1]);
            U2 = pobs([x y2]);
            U3 = pobs([x y3]);
            U4 = pobs([x y4]);

            C1 = empcopulacdf(U1,K,'deheuvels'); 
            C2 = empcopulacdf(U2,K,'deheuvels');
            C3 = empcopulacdf(U3,K,'deheuvels');
            C4 = empcopulacdf(U4,K,'deheuvels');
            M_dist_linear = M_dist_linear + empCopulaDistance(C1,M_copula,'sse');
            W_dist_linear = W_dist_linear + empCopulaDistance(C1,W_copula,'sse');
            eig_linear = eig_linear + eigValueMetric(U1);
            [~,emdVal] = emd(C1,M_copula,weightVec,weightVec,@gdf);
            emd_M_linear = emd_M_linear + emdVal;
            [~,emdVal] = emd(C1,W_copula,weightVec,weightVec,@gdf);
            emd_W_linear = emd_W_linear + emdVal;
            
            M_dist_quadratic = M_dist_quadratic + empCopulaDistance(C2,M_copula,'sse');
            W_dist_quadratic = W_dist_quadratic + empCopulaDistance(C2,W_copula,'sse');
            eig_quadratic = eig_quadratic + eigValueMetric(U2);
            [~,emdVal] = emd(C2,M_copula,weightVec,weightVec,@gdf);
            emd_M_quadratic = emd_M_quadratic + emdVal;
            [~,emdVal] = emd(C2,W_copula,weightVec,weightVec,@gdf);
            emd_W_quadratic = emd_W_quadratic + emdVal;

            M_dist_cubic = M_dist_cubic + empCopulaDistance(C3,M_copula,'sse');
            W_dist_cubic = W_dist_cubic + empCopulaDistance(C3,W_copula,'sse');
            eig_cubic = eig_cubic + eigValueMetric(U3);
            [~,emdVal] = emd(C3,M_copula,weightVec,weightVec,@gdf);
            emd_M_cubic = emd_M_cubic + emdVal;
            [~,emdVal] = emd(C3,W_copula,weightVec,weightVec,@gdf);
            emd_W_cubic = emd_W_cubic + emdVal;

            M_dist_sinusoid = M_dist_sinusoid + empCopulaDistance(C4,M_copula,'sse');
            W_dist_sinusoid = W_dist_sinusoid + empCopulaDistance(C4,W_copula,'sse');
            eig_sinusoid = eig_sinusoid + eigValueMetric(U4);
            [~,emdVal] = emd(C4,M_copula,weightVec,weightVec,@gdf);
            emd_M_sinusoid = emd_M_sinusoid + emdVal;
            [~,emdVal] = emd(C4,W_copula,weightVec,weightVec,@gdf);
            emd_W_sinusoid = emd_W_sinusoid + emdVal;
        end

        distancesTensor(1,l) = M_dist_linear/numMCSim;
        distancesTensor(2,l) = W_dist_linear/numMCSim;
        distancesTensor(3,l) = eig_linear/numMCSim;
        distancesTensor(4,l) = emd_M_linear/numMCSim;
        distancesTensor(5,l) = emd_W_linear/numMCSim;

        distancesVec_parabola(1,l) = M_dist_quadratic/numMCSim;
        distancesVec_parabola(2,l) = W_dist_quadratic/numMCSim;
        distancesVec_parabola(3,l) = eig_quadratic/numMCSim;
        distancesVec_parabola(4,l) = emd_M_quadratic/numMCSim;
        distancesVec_parabola(5,l) = emd_W_quadratic/numMCSim;

        distancesVec_cubic(1,l) = M_dist_sinusoid/numMCSim;
        distancesVec_cubic(2,l) = W_dist_cubic/numMCSim;
        distancesVec_cubic(3,l) = eig_cubic/numMCSim;
        distancesVec_cubic(4,l) = emd_M_cubic/numMCSim;
        distancesVec_cubic(5,l) = emd_W_cubic/numMCSim;

        distancesVec_sinusoid(1,l) = M_dist_sinusoid/numMCSim;
        distancesVec_sinusoid(2,l) = W_dist_sinusoid/numMCSim;
        distancesVec_sinusoid(3,l) = eig_sinusoid/numMCSim;
        distancesVec_sinusoid(4,l) = emd_M_sinusoid/numMCSim;
        distancesVec_sinusoid(5,l) = emd_W_sinusoid/numMCSim;
    end

    % save the data
    if(ispc)
        save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\distVec_M_%d.mat', M));
    elseif(ismac)
        save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/clustering/distVec_M_%d.mat', M));
    else
        save(sprintf('/home/kiran/ownCloud/PhD/sim_results/clustering/distVec_M_%d.mat', M));
    end
end

%%

M = 100;
if(ispc)
    load(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\distVec_M_%d.mat', M));
elseif(ismac)
    load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/clustering/distVec_M_%d.mat', M));
else
    load(sprintf('/home/kiran/ownCloud/PhD/sim_results/clustering/distVec_M_%d.mat', M));
end

figure;
subplot(2,2,1);
plot(noiseVec,distancesTensor(1,:),'b.-', ...
     noiseVec,distancesTensor(2,:),'r.-', ...
     noiseVec,distancesTensor(3,:),'k.-');
grid on;
xlabel('noise');
ylabel('Metric');
legend('SSE[M(u,v)]','SSE[W(u,v)]','EIG');
title(sprintf('Linear Dependency M=%d',M));

subplot(2,2,2);
plot(noiseVec,distancesVec_parabola(1,:),'o-.', ...
     noiseVec,distancesVec_parabola(2,:),'+-.', ...
     noiseVec,distancesVec_parabola(3,:),'d-.');
grid on;
xlabel('noise');
ylabel('Metric');
legend('SSE[M(u,v)]','SSE[W(u,v)]','EIG');
title(sprintf('Quadratic Dependency M=%d',M));

subplot(2,2,3);
plot(noiseVec,distancesVec_cubic(1,:),'o-.', ...
     noiseVec,distancesVec_cubic(2,:),'+-.', ...
     noiseVec,distancesVec_cubic(3,:),'d-.');
grid on;
xlabel('noise');
ylabel('Metric');
legend('SSE[M(u,v)]','SSE[W(u,v)]','EIG');
title(sprintf('Cubic Dependency M=%d',M));

subplot(2,2,4);
plot(noiseVec,distancesVec_sinusoid(1,:),'o-.', ...
     noiseVec,distancesVec_sinusoid(2,:),'+-.', ...
     noiseVec,distancesVec_sinusoid(3,:),'d-.');
grid on;
xlabel('noise');
ylabel('Metric');
legend('SSE[M(u,v)]','SSE[W(u,v)]','EIG');
title(sprintf('Sinusoidal Dependency M=%d',M));

figure;
subplot(1,3,1);
plot(noiseVec,distancesTensor(1,:),'o-.', ...
     noiseVec,distancesVec_parabola(1,:),'+-.', ...
     noiseVec,distancesVec_cubic(1,:),'d-.', ...
     noiseVec,distancesVec_sinusoid(1,:),'v-.');
grid on;
xlabel('noise');
ylabel('Metric');
legend('Linear','Parabola','Cubic','Sinusoid');
title(sprintf('SSE[M(u,v)] M=%d',M));

subplot(1,3,2);
plot(noiseVec,distancesTensor(2,:),'o-.', ...
     noiseVec,distancesVec_parabola(2,:),'+-.', ...
     noiseVec,distancesVec_cubic(2,:),'d-.', ...
     noiseVec,distancesVec_sinusoid(2,:),'v-.');
grid on;
xlabel('noise');
ylabel('Metric');
legend('Linear','Parabola','Cubic','Sinusoid');
title(sprintf('SSE[W(u,v)] M=%d',M));

subplot(1,3,3);
plot(noiseVec,distancesTensor(3,:),'o-.', ...
     noiseVec,distancesVec_parabola(3,:),'+-.', ...
     noiseVec,distancesVec_cubic(3,:),'d-.', ...
     noiseVec,distancesVec_sinusoid(3,:),'v-.');
grid on;
xlabel('noise');
ylabel('Metric');
legend('Linear','Parabola','Cubic','Sinusoid');
title(sprintf('EIG M=%d',M));