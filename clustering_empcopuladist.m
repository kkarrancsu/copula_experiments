%% configuration
clear;
clc;

xMin = 0;
xMax = 1;
num_noise = 30;                    % The number of different noise levels used
noise = 3;                         % A constant to determine the amount of noise

MVecResults = 100:100:1000;

options=optimset('Display', 'off');

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');

K = 50;
weightVec = ones(K,1);
numMCSim = 100;

for M=MVecResults
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
    distancesTensor = zeros(3,length(noiseVec));
    distancesVec_parabola = zeros(3,length(noiseVec));
    distancesVec_cubic = zeros(3,length(noiseVec));
    distancesVec_sinusoid = zeros(3,length(noiseVec));
                    % 1 - SSE to M copula
                    % 2 - SSE to W copula
                    % 3 - eigenvalue metric
    for l=noiseVec
        dispstat(sprintf('Computing for noise level=%d >> M=%d',l,M),'keepthis', 'timestamp');
        M_dist_linear = 0; M_dist_quadratic = 0; M_dist_cubic = 0; M_dist_sinusoid = 0;
        W_dist_linear = 0; W_dist_quadratic = 0; W_dist_cubic = 0; W_dist_sinusoid = 0;
        eig_linear = 0; eig_quadratic = 0; eig_cubic = 0; eig_sinusoid = 0;
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
            
            M_dist_quadratic = M_dist_quadratic + empCopulaDistance(C2,M_copula,'sse');
            W_dist_quadratic = W_dist_quadratic + empCopulaDistance(C2,W_copula,'sse');
            eig_quadratic = eig_quadratic + eigValueMetric(U2);

            M_dist_cubic = M_dist_cubic + empCopulaDistance(C3,M_copula,'sse');
            W_dist_cubic = W_dist_cubic + empCopulaDistance(C3,W_copula,'sse');
            eig_cubic = eig_cubic + eigValueMetric(U3);

            M_dist_sinusoid = M_dist_sinusoid + empCopulaDistance(C4,M_copula,'sse');
            W_dist_sinusoid = W_dist_sinusoid + empCopulaDistance(C4,W_copula,'sse');
            eig_sinusoid = eig_sinusoid + eigValueMetric(U4);
        end

        distancesTensor(1,l) = M_dist_linear/numMCSim;
        distancesTensor(2,l) = W_dist_linear/numMCSim;
        distancesTensor(3,l) = eig_linear/numMCSim;

        distancesVec_parabola(1,l) = M_dist_quadratic/numMCSim;
        distancesVec_parabola(2,l) = W_dist_quadratic/numMCSim;
        distancesVec_parabola(3,l) = eig_quadratic/numMCSim;

        distancesVec_cubic(1,l) = M_dist_sinusoid/numMCSim;
        distancesVec_cubic(2,l) = W_dist_cubic/numMCSim;
        distancesVec_cubic(3,l) = eig_cubic/numMCSim;

        distancesVec_sinusoid(1,l) = M_dist_sinusoid/numMCSim;
        distancesVec_sinusoid(2,l) = W_dist_sinusoid/numMCSim;
        distancesVec_sinusoid(3,l) = eig_sinusoid/numMCSim;
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

M = 500;
if(ispc)
    load(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\distVec_M_%d.mat', M));
elseif(ismac)
    load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/clustering/distVec_M_%d.mat', M));
else
    load(sprintf('/home/kiran/ownCloud/PhD/sim_results/clustering/distVec_M_%d.mat', M));
end

figure;
subplot(2,2,1);
hh1 = plot(noiseVec,distancesTensor(1,:),'o-.', ...
     noiseVec,distancesTensor(2,:),'+-.', ...
     noiseVec,distancesTensor(3,:),'d-.');
grid on;
xlabel('noise');
ylabel('Metric');
legend('SSE[M(u,v)]','SSE[W(u,v)]','EIG');
title(sprintf('Linear Dependency M=%d',M));
hh1(1).LineWidth = 1.5; 
hh1(2).LineWidth = 1.5; 
hh1(3).LineWidth = 1.5; 

subplot(2,2,2);
hh1 = plot(noiseVec,distancesVec_parabola(1,:),'o-.', ...
     noiseVec,distancesVec_parabola(2,:),'+-.', ...
     noiseVec,distancesVec_parabola(3,:),'d-.');
grid on;
xlabel('noise');
ylabel('Metric');
legend('SSE[M(u,v)]','SSE[W(u,v)]','EIG');
title(sprintf('Quadratic Dependency M=%d',M));
hh1(1).LineWidth = 1.5; 
hh1(2).LineWidth = 1.5; 
hh1(3).LineWidth = 1.5; 

subplot(2,2,3);
hh1 = plot(noiseVec,distancesVec_cubic(1,:),'o-.', ...
     noiseVec,distancesVec_cubic(2,:),'+-.', ...
     noiseVec,distancesVec_cubic(3,:),'d-.');
grid on;
xlabel('noise');
ylabel('Metric');
legend('SSE[M(u,v)]','SSE[W(u,v)]','EIG');
title(sprintf('Cubic Dependency M=%d',M));
hh1(1).LineWidth = 1.5; 
hh1(2).LineWidth = 1.5; 
hh1(3).LineWidth = 1.5; 

subplot(2,2,4);
hh1 = plot(noiseVec,distancesVec_sinusoid(1,:),'o-.', ...
     noiseVec,distancesVec_sinusoid(2,:),'+-.', ...
     noiseVec,distancesVec_sinusoid(3,:),'d-.');
grid on;
xlabel('noise');
ylabel('Metric');
legend('SSE[M(u,v)]','SSE[W(u,v)]','EIG');
title(sprintf('Sinusoidal Dependency M=%d',M));
hh1(1).LineWidth = 1.5; 
hh1(2).LineWidth = 1.5; 
hh1(3).LineWidth = 1.5; 

figure;
subplot(1,3,1);
hh1 = plot(noiseVec,distancesTensor(1,:),'o-.', ...
     noiseVec,distancesVec_parabola(1,:),'+-.', ...
     noiseVec,distancesVec_cubic(1,:),'d-.', ...
     noiseVec,distancesVec_sinusoid(1,:),'v-.');
grid on;
xlabel('noise');
ylabel('Metric');
legend('Linear','Parabola','Cubic','Sinusoid');
title(sprintf('SSE[M(u,v)] M=%d',M));
hh1(1).LineWidth = 1.5; 
hh1(2).LineWidth = 1.5; 
hh1(3).LineWidth = 1.5; 
hh1(4).LineWidth = 1.5; 

subplot(1,3,2);
hh1 = plot(noiseVec,distancesTensor(2,:),'o-.', ...
     noiseVec,distancesVec_parabola(2,:),'+-.', ...
     noiseVec,distancesVec_cubic(2,:),'d-.', ...
     noiseVec,distancesVec_sinusoid(2,:),'v-.');
grid on;
xlabel('noise');
ylabel('Metric');
legend('Linear','Parabola','Cubic','Sinusoid');
title(sprintf('SSE[W(u,v)] M=%d',M));
hh1(1).LineWidth = 1.5; 
hh1(2).LineWidth = 1.5; 
hh1(3).LineWidth = 1.5; 
hh1(4).LineWidth = 1.5; 

subplot(1,3,3);
hh1 = plot(noiseVec,distancesTensor(3,:),'o-.', ...
     noiseVec,distancesVec_parabola(3,:),'+-.', ...
     noiseVec,distancesVec_cubic(3,:),'d-.', ...
     noiseVec,distancesVec_sinusoid(3,:),'v-.');
grid on;
xlabel('noise');
ylabel('Metric');
legend('Linear','Parabola','Cubic','Sinusoid');
title(sprintf('EIG M=%d',M));
hh1(1).LineWidth = 1.5; 
hh1(2).LineWidth = 1.5; 
hh1(3).LineWidth = 1.5; 
hh1(4).LineWidth = 1.5; 

%% Compute the first difference in the V dimension, after converting to pobs
% to see if we can use that as a measure of whether it is "function" or an
% "association"
clear;
clc;

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');

xMin = 0;
xMax = 1;
num_noise = 30;                    % The number of different noise levels used
noise = 3;                         % A constant to determine the amount of noise
numNoisePtsToCompute = 10;

numMCSims = 1000;
MVecResults = 100:100:1000;
num_noise_test_min = 1;
num_noise_test_max = 30;

M = MVecResults(1);
noiseVec = linspace(num_noise_test_min,num_noise_test_max,numNoisePtsToCompute);

noiseIdx = 1;
v1_sig = zeros(numNoisePtsToCompute,M-1); v2_sig = zeros(numNoisePtsToCompute,M-1); 
v3_sig = zeros(numNoisePtsToCompute,M-1); v4_sig = zeros(numNoisePtsToCompute,M-1);
v5_sig = zeros(numNoisePtsToCompute,M-1); v6_sig = zeros(numNoisePtsToCompute,M-1); 
v7_sig = zeros(numNoisePtsToCompute,M-1); v8_sig = zeros(numNoisePtsToCompute,M-1);

for l=noiseVec
    dispstat(sprintf('Computing for noise level=%0.02f >> M=%d',l,M),'keepthis', 'timestamp');
    v1_tmp_sig = zeros(1,M-1); v2_tmp_sig = zeros(1,M-1);
    v3_tmp_sig = zeros(1,M-1); v4_tmp_sig = zeros(1,M-1);
    v5_tmp_sig = zeros(1,M-1); v6_tmp_sig = zeros(1,M-1);
    v7_tmp_sig = zeros(1,M-1); v8_tmp_sig = zeros(1,M-1);
    parfor mcSimNum=1:numMCSims
        x = rand(M,1)*(xMax-xMin)+xMin;
        y1 = x + noise*(l/num_noise)*randn(M,1);
        y2 = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
        y3 = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
        y4 = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
        y5 = sin(16*pi*x) + noise*(l/num_noise)*randn(M,1);
        y6 = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
        y7=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
        y8 = (x > 0.5) + noise*5*l/num_noise*randn(M,1);

        U1 = pobs([x y1]); [~,I1] = sort(U1(:,1)); zz1 = U1(:,2); vv1 = zz1(I1);
        U2 = pobs([x y2]); [~,I2] = sort(U2(:,1)); zz2 = U2(:,2); vv2 = zz2(I2);
        U3 = pobs([x y3]); [~,I3] = sort(U3(:,1)); zz3 = U3(:,2); vv3 = zz3(I3);
        U4 = pobs([x y4]); [~,I4] = sort(U4(:,1)); zz4 = U4(:,2); vv4 = zz4(I4);
        U5 = pobs([x y5]); [~,I5] = sort(U5(:,1)); zz5 = U5(:,2); vv5 = zz5(I5);
        U6 = pobs([x y6]); [~,I6] = sort(U6(:,1)); zz6 = U6(:,2); vv6 = zz6(I6);
        U7 = pobs([x y7]); [~,I7] = sort(U7(:,1)); zz7 = U7(:,2); vv7 = zz7(I7);
        U8 = pobs([x y8]); [~,I8] = sort(U8(:,1)); zz8 = U8(:,2); vv8 = zz8(I8);
        
        v1_tmp_sig = v1_tmp_sig + diff(vv1)';
        v2_tmp_sig = v2_tmp_sig + diff(vv2)';
        v3_tmp_sig = v3_tmp_sig + diff(vv3)';
        v4_tmp_sig = v4_tmp_sig + diff(vv4)';
        v5_tmp_sig = v5_tmp_sig + diff(vv5)';
        v6_tmp_sig = v6_tmp_sig + diff(vv6)';
        v7_tmp_sig = v7_tmp_sig + diff(vv7)';
        v8_tmp_sig = v8_tmp_sig + diff(vv8)';
    end
    
    v1_sig(noiseIdx,:) = v1_tmp_sig/numMCSims; v2_sig(noiseIdx,:) = v2_tmp_sig/numMCSims;
    v3_sig(noiseIdx,:) = v3_tmp_sig/numMCSims; v4_sig(noiseIdx,:) = v4_tmp_sig/numMCSims;
    v5_sig(noiseIdx,:) = v5_tmp_sig/numMCSims; v6_sig(noiseIdx,:) = v6_tmp_sig/numMCSims;
    v7_sig(noiseIdx,:) = v7_tmp_sig/numMCSims; v8_sig(noiseIdx,:) = v8_tmp_sig/numMCSims;
    
    noiseIdx = noiseIdx + 1;
end

idxToPlot = floor(numNoisePtsToCompute/2);

subplot(2,4,1);
hh1 = plot(1:M-1,v1_sig(idxToPlot,:),'o-.');
hh1(1).LineWidth = 1.5; 
grid on; title('Linear');

subplot(2,4,2);
hh1 = plot(1:M-1,v2_sig(idxToPlot,:),'o-.');
hh1(1).LineWidth = 1.5; 
grid on; title('Quadratic');

subplot(2,4,3);
hh1 = plot(1:M-1,v3_sig(idxToPlot,:),'o-.');
hh1(1).LineWidth = 1.5; 
grid on; title('Cubic');

subplot(2,4,4);
hh1 = plot(1:M-1,v4_sig(idxToPlot,:),'o-.');
hh1(1).LineWidth = 1.5; 
grid on; title('Low-Freq Sin');

subplot(2,4,5);
hh1 = plot(1:M-1,v5_sig(idxToPlot,:),'o-.');
hh1(1).LineWidth = 1.5; 
grid on; title('High-Freq Sin');

subplot(2,4,6);
hh1 = plot(1:M-1,v6_sig(idxToPlot,:),'o-.');
hh1(1).LineWidth = 1.5; 
grid on; title('Fourth-Root');

subplot(2,4,7);
hh1 = plot(1:M-1,v7_sig(idxToPlot,:),'o-.');
hh1(1).LineWidth = 1.5; 
grid on; title('Circular');

subplot(2,4,8);
hh1 = plot(1:M-1,v8_sig(idxToPlot,:),'o-.');
hh1(1).LineWidth = 1.5; 
grid on; title('Step');

%% Simply plot what the data would look like

clear;
clc;

M = 500;

xMin = 0;
xMax = 1;
num_noise = 30;                    % The number of different noise levels used
noise = 3;                         % A constant to determine the amount of noise

l = 1;

x = rand(M,1)*(xMax-xMin)+xMin;
y1 = x + noise*(l/num_noise)*randn(M,1);
y2 = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
y3 = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
y4 = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
y5 = sin(16*pi*x) + noise*(l/num_noise)*randn(M,1);
y6 = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
y7=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
y8 = (x > 0.5) + noise*5*l/num_noise*randn(M,1);
y9 = rand(M,1)*(xMax-xMin)+xMin;

U1 = pobs([x y1]);
U2 = pobs([x y2]);
U3 = pobs([x y3]);
U4 = pobs([x y4]);
U5 = pobs([x y5]);
U6 = pobs([x y6]);
U7 = pobs([x y7]);
U8 = pobs([x y8]);
U9 = pobs([x y9]);

figure;
subplot(3,3,1); scatter(U1(:,1),U1(:,2)); grid on; title('Linear');
subplot(3,3,2); scatter(U2(:,1),U2(:,2)); grid on; title('Quadratic');
subplot(3,3,3); scatter(U3(:,1),U3(:,2)); grid on; title('Cubic');
subplot(3,3,4); scatter(U4(:,1),U4(:,2)); grid on; title('Low-Freq Sin');
subplot(3,3,5); scatter(U5(:,1),U5(:,2)); grid on; title('High-Freq Sin');
subplot(3,3,6); scatter(U6(:,1),U6(:,2)); grid on; title('Fourth-Root');
subplot(3,3,7); scatter(U7(:,1),U7(:,2)); grid on; title('Circular');
subplot(3,3,8); scatter(U8(:,1),U8(:,2)); grid on; title('Step');
subplot(3,3,9); scatter(U9(:,1),U9(:,2)); grid on; title('Independence');
figtitle(sprintf('Noise=%0.02f\n', l));

%% E-Copula experiments, why does distance never goto 0?

clear;
clc;

M = 200;
x = rand(M,1);
y = x;
[U1,V1] = pobs_sorted(x,y); UU1 = [U1 V1];
C_uv = ecopula(U1,V1);
M_uv = min(U1,V1)';
W_uv = max(U1+V1-1,0)';

d1 = empCopulaDistance(C_uv,M_uv,'sse')
d2 = empCopulaDistance(C_uv,W_uv,'sse')

plot(1:M,C_uv,1:M,M_uv,1:M,W_uv); grid on;
legend('C','M','W')


%% test how adding new points to the empirical copula affects the SSE distance metric

clear;
clc;

xMin = 0;
xMax = 1;
num_noise = 30;                    % The number of different noise levels used
noise = 3;                         % A constant to determine the amount of noise

MVecResults = 100:100:1000;
num_noise_test_min = 0;
num_noise_test_max = 30;
noiseVec = num_noise_test_min:num_noise_test_max;
numMCSim = 100;
numDepTypes = 9;

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');

distCellM = cell(length(noiseVec),numDepTypes);
distCellW = cell(length(noiseVec),numDepTypes);
for M=MVecResults
    lIdx = 1;
    for l=noiseVec
        distVecM_dep1_mc = zeros(1,M-1); distVecW_dep1_mc = zeros(1,M-1);
        distVecM_dep2_mc = zeros(1,M-1); distVecW_dep2_mc = zeros(1,M-1);
        distVecM_dep3_mc = zeros(1,M-1); distVecW_dep3_mc = zeros(1,M-1);
        distVecM_dep4_mc = zeros(1,M-1); distVecW_dep4_mc = zeros(1,M-1);
        distVecM_dep5_mc = zeros(1,M-1); distVecW_dep5_mc = zeros(1,M-1);
        distVecM_dep6_mc = zeros(1,M-1); distVecW_dep6_mc = zeros(1,M-1);
        distVecM_dep7_mc = zeros(1,M-1); distVecW_dep7_mc = zeros(1,M-1);
        distVecM_dep8_mc = zeros(1,M-1); distVecW_dep8_mc = zeros(1,M-1);
        distVecM_dep9_mc = zeros(1,M-1); distVecW_dep9_mc = zeros(1,M-1);

        dispstat(sprintf('Computing for noise level=%d >> M=%d',l,M),'keepthis', 'timestamp');
        
        for mcSimNum=1:numMCSim
            x = rand(M,1)*(xMax-xMin)+xMin;
            y1 = x + noise*(l/num_noise)*randn(M,1);
            y2 = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
            y3 = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
            y4 = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
            y5 = sin(16*pi*x) + noise*(l/num_noise)*randn(M,1);
            y6 = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
            y7=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
            y8 = (x > 0.5) + noise*5*l/num_noise*randn(M,1);
            y9 = rand(M,1)*(xMax-xMin)+xMin;

            [U1,V1] = pobs_sorted(x,y1); UU1 = [U1 V1];
            [U2,V2] = pobs_sorted(x,y2); UU2 = [U2 V2];
            [U3,V3] = pobs_sorted(x,y3); UU3 = [U3 V3];
            [U4,V4] = pobs_sorted(x,y4); UU4 = [U4 V4];
            [U5,V5] = pobs_sorted(x,y5); UU5 = [U5 V5];
            [U6,V6] = pobs_sorted(x,y6); UU6 = [U6 V6];
            [U7,V7] = pobs_sorted(x,y7); UU7 = [U7 V7];
            [U8,V8] = pobs_sorted(x,y8); UU8 = [U8 V8];
            [U9,V9] = pobs_sorted(x,y9); UU9 = [U9 V9];
            
            distVecM_dep1 = zeros(1,M); distVecW_dep1 = zeros(1,M);
            distVecM_dep2 = zeros(1,M); distVecW_dep2 = zeros(1,M);
            distVecM_dep3 = zeros(1,M); distVecW_dep3 = zeros(1,M);
            distVecM_dep4 = zeros(1,M); distVecW_dep4 = zeros(1,M);
            distVecM_dep5 = zeros(1,M); distVecW_dep5 = zeros(1,M);
            distVecM_dep6 = zeros(1,M); distVecW_dep6 = zeros(1,M);
            distVecM_dep7 = zeros(1,M); distVecW_dep7 = zeros(1,M);
            distVecM_dep8 = zeros(1,M); distVecW_dep8 = zeros(1,M);
            distVecM_dep9 = zeros(1,M); distVecW_dep9 = zeros(1,M);
            
            parfor ii=2:M
                pts = UU1(1:ii,:); ptsU = pts(:,1); ptsV = pts(:,2);
                C_uv = ecopula(ptsU,ptsV);
                M_uv = min(ptsU,ptsV)';
                W_uv = max(ptsU+ptsV-1,0)';
                distVecM_dep1(ii) = empCopulaDistance(C_uv,M_uv,'sse');
                distVecW_dep1(ii) = empCopulaDistance(C_uv,W_uv,'sse');
                
                pts = UU2(1:ii,:); ptsU = pts(:,1); ptsV = pts(:,2);
                C_uv = ecopula(ptsU,ptsV);
                M_uv = min(ptsU,ptsV)';
                W_uv = max(ptsU+ptsV-1,0)';
                distVecM_dep2(ii) = empCopulaDistance(C_uv,M_uv,'sse');
                distVecW_dep2(ii) = empCopulaDistance(C_uv,W_uv,'sse');
                
                pts = UU3(1:ii,:); ptsU = pts(:,1); ptsV = pts(:,2);
                C_uv = ecopula(ptsU,ptsV);
                M_uv = min(ptsU,ptsV)';
                W_uv = max(ptsU+ptsV-1,0)';
                distVecM_dep3(ii) = empCopulaDistance(C_uv,M_uv,'sse');
                distVecW_dep3(ii) = empCopulaDistance(C_uv,W_uv,'sse');
                
                pts = UU4(1:ii,:); ptsU = pts(:,1); ptsV = pts(:,2);
                C_uv = ecopula(ptsU,ptsV);
                M_uv = min(ptsU,ptsV)';
                W_uv = max(ptsU+ptsV-1,0)';
                distVecM_dep4(ii) = empCopulaDistance(C_uv,M_uv,'sse');
                distVecW_dep4(ii) = empCopulaDistance(C_uv,W_uv,'sse');
                
                pts = UU5(1:ii,:); ptsU = pts(:,1); ptsV = pts(:,2);
                C_uv = ecopula(ptsU,ptsV);
                M_uv = min(ptsU,ptsV)';
                W_uv = max(ptsU+ptsV-1,0)';
                distVecM_dep5(ii) = empCopulaDistance(C_uv,M_uv,'sse');
                distVecW_dep5(ii) = empCopulaDistance(C_uv,W_uv,'sse');
                
                pts = UU6(1:ii,:); ptsU = pts(:,1); ptsV = pts(:,2);
                C_uv = ecopula(ptsU,ptsV);
                M_uv = min(ptsU,ptsV)';
                W_uv = max(ptsU+ptsV-1,0)';
                distVecM_dep6(ii) = empCopulaDistance(C_uv,M_uv,'sse');
                distVecW_dep6(ii) = empCopulaDistance(C_uv,W_uv,'sse');
                
                pts = UU7(1:ii,:); ptsU = pts(:,1); ptsV = pts(:,2);
                C_uv = ecopula(ptsU,ptsV);
                M_uv = min(ptsU,ptsV)';
                W_uv = max(ptsU+ptsV-1,0)';
                distVecM_dep7(ii) = empCopulaDistance(C_uv,M_uv,'sse');
                distVecW_dep7(ii) = empCopulaDistance(C_uv,W_uv,'sse');
                
                pts = UU8(1:ii,:); ptsU = pts(:,1); ptsV = pts(:,2);
                C_uv = ecopula(ptsU,ptsV);
                M_uv = min(ptsU,ptsV)';
                W_uv = max(ptsU+ptsV-1,0)';
                distVecM_dep8(ii) = empCopulaDistance(C_uv,M_uv,'sse');
                distVecW_dep8(ii) = empCopulaDistance(C_uv,W_uv,'sse');
                
                pts = UU9(1:ii,:); ptsU = pts(:,1); ptsV = pts(:,2);
                C_uv = ecopula(ptsU,ptsV);
                M_uv = min(ptsU,ptsV)';
                W_uv = max(ptsU+ptsV-1,0)';
                distVecM_dep9(ii) = empCopulaDistance(C_uv,M_uv,'sse');
                distVecW_dep9(ii) = empCopulaDistance(C_uv,W_uv,'sse');
            end
            distVecM_dep1 = distVecM_dep1(2:end); distVecW_dep1 = distVecW_dep1(2:end);
            distVecM_dep2 = distVecM_dep2(2:end); distVecW_dep2 = distVecW_dep2(2:end);
            distVecM_dep3 = distVecM_dep3(2:end); distVecW_dep3 = distVecW_dep3(2:end);
            distVecM_dep4 = distVecM_dep4(2:end); distVecW_dep4 = distVecW_dep4(2:end);
            distVecM_dep5 = distVecM_dep5(2:end); distVecW_dep5 = distVecW_dep5(2:end);
            distVecM_dep6 = distVecM_dep6(2:end); distVecW_dep6 = distVecW_dep6(2:end);
            distVecM_dep7 = distVecM_dep7(2:end); distVecW_dep7 = distVecW_dep7(2:end);
            distVecM_dep8 = distVecM_dep8(2:end); distVecW_dep8 = distVecW_dep8(2:end);
            distVecM_dep9 = distVecM_dep9(2:end); distVecW_dep9 = distVecW_dep9(2:end);
            
            distVecM_dep1_mc = distVecM_dep1_mc + distVecM_dep1; distVecW_dep1_mc = distVecW_dep1_mc + distVecW_dep1;
            distVecM_dep2_mc = distVecM_dep2_mc + distVecM_dep2; distVecW_dep2_mc = distVecW_dep2_mc + distVecW_dep2;
            distVecM_dep3_mc = distVecM_dep3_mc + distVecM_dep3; distVecW_dep3_mc = distVecW_dep3_mc + distVecW_dep3;
            distVecM_dep4_mc = distVecM_dep4_mc + distVecM_dep4; distVecW_dep4_mc = distVecW_dep4_mc + distVecW_dep4;
            distVecM_dep5_mc = distVecM_dep5_mc + distVecM_dep5; distVecW_dep5_mc = distVecW_dep5_mc + distVecW_dep5;
            distVecM_dep6_mc = distVecM_dep6_mc + distVecM_dep6; distVecW_dep6_mc = distVecW_dep6_mc + distVecW_dep6;
            distVecM_dep7_mc = distVecM_dep7_mc + distVecM_dep7; distVecW_dep7_mc = distVecW_dep7_mc + distVecW_dep7;
            distVecM_dep8_mc = distVecM_dep8_mc + distVecM_dep8; distVecW_dep8_mc = distVecW_dep8_mc + distVecW_dep8;
            distVecM_dep9_mc = distVecM_dep9_mc + distVecM_dep9; distVecW_dep9_mc = distVecW_dep9_mc + distVecW_dep9;
            
        end
        
        distVecM_dep1_mc = distVecM_dep1_mc/numMCSim; distVecW_dep1_mc = distVecW_dep1_mc/numMCSim;
        distVecM_dep2_mc = distVecM_dep2_mc/numMCSim; distVecW_dep2_mc = distVecW_dep2_mc/numMCSim;
        distVecM_dep3_mc = distVecM_dep3_mc/numMCSim; distVecW_dep3_mc = distVecW_dep3_mc/numMCSim;
        distVecM_dep4_mc = distVecM_dep4_mc/numMCSim; distVecW_dep4_mc = distVecW_dep4_mc/numMCSim;
        distVecM_dep5_mc = distVecM_dep5_mc/numMCSim; distVecW_dep5_mc = distVecW_dep5_mc/numMCSim;
        distVecM_dep6_mc = distVecM_dep6_mc/numMCSim; distVecW_dep6_mc = distVecW_dep6_mc/numMCSim;
        distVecM_dep7_mc = distVecM_dep7_mc/numMCSim; distVecW_dep7_mc = distVecW_dep7_mc/numMCSim;
        distVecM_dep8_mc = distVecM_dep8_mc/numMCSim; distVecW_dep8_mc = distVecW_dep8_mc/numMCSim;
        distVecM_dep9_mc = distVecM_dep9_mc/numMCSim; distVecW_dep9_mc = distVecW_dep9_mc/numMCSim;
        
        distCellM{lIdx,1} = distVecM_dep1_mc; distCellW{lIdx,1} = distVecW_dep1_mc;
        distCellM{lIdx,2} = distVecM_dep2_mc; distCellW{lIdx,2} = distVecW_dep2_mc;
        distCellM{lIdx,3} = distVecM_dep3_mc; distCellW{lIdx,3} = distVecW_dep3_mc;
        distCellM{lIdx,4} = distVecM_dep4_mc; distCellW{lIdx,4} = distVecW_dep4_mc;
        distCellM{lIdx,5} = distVecM_dep5_mc; distCellW{lIdx,5} = distVecW_dep5_mc;
        distCellM{lIdx,6} = distVecM_dep6_mc; distCellW{lIdx,6} = distVecW_dep6_mc;
        distCellM{lIdx,7} = distVecM_dep7_mc; distCellW{lIdx,7} = distVecW_dep7_mc;
        distCellM{lIdx,8} = distVecM_dep8_mc; distCellW{lIdx,8} = distVecW_dep8_mc;
        distCellM{lIdx,9} = distVecM_dep9_mc; distCellW{lIdx,9} = distVecW_dep9_mc;
        
        lIdx = lIdx + 1;
    end
    % save the data for post-analysis
    if(ispc)
        save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\pointAddition_M_%d.mat', M));
    elseif(ismac)
        save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/clustering/pointAddition_M_%d.mat', M));
    else
        save(sprintf('/home/kiran/ownCloud/PhD/sim_results/clustering/pointAddition_M_%d.mat', M));
    end
end

%% plot the results for the point-additionp

M = 300;
if(ispc)
    load(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\pointAddition_M_%d.mat', M));
elseif(ismac)
    load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/clustering/pointAddition_M_%d.mat', M));
else
    load(sprintf('/home/kiran/ownCloud/PhD/sim_results/clustering/pointAddition_M_%d.mat', M));
end

noiseLevelsToPlot = [0,15];

depIdx = 1;
subplot(3,3,depIdx);
for noiseIdx=noiseLevelsToPlot
    y = distCellM{noiseIdx+1,depIdx};
    x = 1:length(y);
    plot(x,y,'-*'); hold on;
end
for noiseIdx=noiseLevelsToPlot
    y = distCellW{noiseIdx+1,depIdx};
    x = 1:length(y);
    plot(x,y,'-.v'); hold on;
end
grid on;
xlabel('ii'); ylabel('Distance');
legendCell = cell(1,length(noiseLevelsToPlot)*2);
legendCellIdx = 1;
for jj=1:length(noiseLevelsToPlot)
    for ii=1:length(noiseLevelsToPlot)
        if(jj==1)
            legendCell{legendCellIdx} = sprintf('M(%d)',noiseLevelsToPlot(ii));
        else
            legendCell{legendCellIdx} = sprintf('W(%d)',noiseLevelsToPlot(ii));
        end
        legendCellIdx = legendCellIdx + 1;
    end
end
legend(legendCell);
title('Linear');

depIdx = depIdx + 1;
subplot(3,3,depIdx);
for noiseIdx=noiseLevelsToPlot
    y = distCellM{noiseIdx+1,depIdx};
    x = 1:length(y);
    plot(x,y,'-*'); hold on;
end
for noiseIdx=noiseLevelsToPlot
    y = distCellW{noiseIdx+1,depIdx};
    x = 1:length(y);
    plot(x,y,'-.v'); hold on;
end
grid on;
xlabel('ii'); ylabel('Distance');
legendCell = cell(1,length(noiseLevelsToPlot)*2);
legendCellIdx = 1;
for jj=1:length(noiseLevelsToPlot)
    for ii=1:length(noiseLevelsToPlot)
        if(jj==1)
            legendCell{legendCellIdx} = sprintf('M(%d)',noiseLevelsToPlot(ii));
        else
            legendCell{legendCellIdx} = sprintf('W(%d)',noiseLevelsToPlot(ii));
        end
        legendCellIdx = legendCellIdx + 1;
    end
end
legend(legendCell);
title('Quadratic');

depIdx = depIdx + 1;
subplot(3,3,depIdx);
for noiseIdx=noiseLevelsToPlot
    y = distCellM{noiseIdx+1,depIdx};
    x = 1:length(y);
    plot(x,y,'-*'); hold on;
end
for noiseIdx=noiseLevelsToPlot
    y = distCellW{noiseIdx+1,depIdx};
    x = 1:length(y);
    plot(x,y,'-.v'); hold on;
end
grid on;
xlabel('ii'); ylabel('Distance');
legendCell = cell(1,length(noiseLevelsToPlot)*2);
legendCellIdx = 1;
for jj=1:length(noiseLevelsToPlot)
    for ii=1:length(noiseLevelsToPlot)
        if(jj==1)
            legendCell{legendCellIdx} = sprintf('M(%d)',noiseLevelsToPlot(ii));
        else
            legendCell{legendCellIdx} = sprintf('W(%d)',noiseLevelsToPlot(ii));
        end
        legendCellIdx = legendCellIdx + 1;
    end
end
legend(legendCell);
title('Cubic');

depIdx = depIdx + 1;
subplot(3,3,depIdx);
for noiseIdx=noiseLevelsToPlot
    y = distCellM{noiseIdx+1,depIdx};
    x = 1:length(y);
    plot(x,y,'-*'); hold on;
end
for noiseIdx=noiseLevelsToPlot
    y = distCellW{noiseIdx+1,depIdx};
    x = 1:length(y);
    plot(x,y,'-.v'); hold on;
end
grid on;
xlabel('ii'); ylabel('Distance');
legendCell = cell(1,length(noiseLevelsToPlot)*2);
legendCellIdx = 1;
for jj=1:length(noiseLevelsToPlot)
    for ii=1:length(noiseLevelsToPlot)
        if(jj==1)
            legendCell{legendCellIdx} = sprintf('M(%d)',noiseLevelsToPlot(ii));
        else
            legendCell{legendCellIdx} = sprintf('W(%d)',noiseLevelsToPlot(ii));
        end
        legendCellIdx = legendCellIdx + 1;
    end
end
legend(legendCell);
title('LF-Sin');

depIdx = depIdx + 1;
subplot(3,3,depIdx);
for noiseIdx=noiseLevelsToPlot
    y = distCellM{noiseIdx+1,depIdx};
    x = 1:length(y);
    plot(x,y,'-*'); hold on;
end
for noiseIdx=noiseLevelsToPlot
    y = distCellW{noiseIdx+1,depIdx};
    x = 1:length(y);
    plot(x,y,'-.v'); hold on;
end
grid on;
xlabel('ii'); ylabel('Distance');
legendCell = cell(1,length(noiseLevelsToPlot)*2);
legendCellIdx = 1;
for jj=1:length(noiseLevelsToPlot)
    for ii=1:length(noiseLevelsToPlot)
        if(jj==1)
            legendCell{legendCellIdx} = sprintf('M(%d)',noiseLevelsToPlot(ii));
        else
            legendCell{legendCellIdx} = sprintf('W(%d)',noiseLevelsToPlot(ii));
        end
        legendCellIdx = legendCellIdx + 1;
    end
end
legend(legendCell);
title('HF-Sin');

depIdx = depIdx + 1;
subplot(3,3,depIdx);
for noiseIdx=noiseLevelsToPlot
    y = distCellM{noiseIdx+1,depIdx};
    x = 1:length(y);
    plot(x,y,'-*'); hold on;
end
for noiseIdx=noiseLevelsToPlot
    y = distCellW{noiseIdx+1,depIdx};
    x = 1:length(y);
    plot(x,y,'-.v'); hold on;
end
grid on;
xlabel('ii'); ylabel('Distance');
legendCell = cell(1,length(noiseLevelsToPlot)*2);
legendCellIdx = 1;
for jj=1:length(noiseLevelsToPlot)
    for ii=1:length(noiseLevelsToPlot)
        if(jj==1)
            legendCell{legendCellIdx} = sprintf('M(%d)',noiseLevelsToPlot(ii));
        else
            legendCell{legendCellIdx} = sprintf('W(%d)',noiseLevelsToPlot(ii));
        end
        legendCellIdx = legendCellIdx + 1;
    end
end
legend(legendCell);
title('Fourth-Root');

depIdx = depIdx + 1;
subplot(3,3,depIdx);
for noiseIdx=noiseLevelsToPlot
    y = distCellM{noiseIdx+1,depIdx};
    x = 1:length(y);
    plot(x,y,'-*'); hold on;
end
for noiseIdx=noiseLevelsToPlot
    y = distCellW{noiseIdx+1,depIdx};
    x = 1:length(y);
    plot(x,y,'-.v'); hold on;
end
grid on;
xlabel('ii'); ylabel('Distance');
legendCell = cell(1,length(noiseLevelsToPlot)*2);
legendCellIdx = 1;
for jj=1:length(noiseLevelsToPlot)
    for ii=1:length(noiseLevelsToPlot)
        if(jj==1)
            legendCell{legendCellIdx} = sprintf('M(%d)',noiseLevelsToPlot(ii));
        else
            legendCell{legendCellIdx} = sprintf('W(%d)',noiseLevelsToPlot(ii));
        end
        legendCellIdx = legendCellIdx + 1;
    end
end
legend(legendCell);
title('Circular');

depIdx = depIdx + 1;
subplot(3,3,depIdx);
for noiseIdx=noiseLevelsToPlot
    y = distCellM{noiseIdx+1,depIdx};
    x = 1:length(y);
    plot(x,y,'-*'); hold on;
end
for noiseIdx=noiseLevelsToPlot
    y = distCellW{noiseIdx+1,depIdx};
    x = 1:length(y);
    plot(x,y,'-.v'); hold on;
end
grid on;
xlabel('ii'); ylabel('Distance');
legendCell = cell(1,length(noiseLevelsToPlot)*2);
legendCellIdx = 1;
for jj=1:length(noiseLevelsToPlot)
    for ii=1:length(noiseLevelsToPlot)
        if(jj==1)
            legendCell{legendCellIdx} = sprintf('M(%d)',noiseLevelsToPlot(ii));
        else
            legendCell{legendCellIdx} = sprintf('W(%d)',noiseLevelsToPlot(ii));
        end
        legendCellIdx = legendCellIdx + 1;
    end
end
legend(legendCell);
title('Step-Function');

depIdx = depIdx + 1;
subplot(3,3,depIdx);
for noiseIdx=noiseLevelsToPlot
    y = distCellM{noiseIdx+1,depIdx};
    x = 1:length(y);
    plot(x,y,'-*'); hold on;
end
for noiseIdx=noiseLevelsToPlot
    y = distCellW{noiseIdx+1,depIdx};
    x = 1:length(y);
    plot(x,y,'-.v'); hold on;
end
grid on;
xlabel('ii'); ylabel('Distance');
legendCell = cell(1,length(noiseLevelsToPlot)*2);
legendCellIdx = 1;
for jj=1:length(noiseLevelsToPlot)
    for ii=1:length(noiseLevelsToPlot)
        if(jj==1)
            legendCell{legendCellIdx} = sprintf('M(%d)',noiseLevelsToPlot(ii));
        else
            legendCell{legendCellIdx} = sprintf('W(%d)',noiseLevelsToPlot(ii));
        end
        legendCellIdx = legendCellIdx + 1;
    end
end
legend(legendCell);
title('Independence');
%% Test how the empirical copula varies over different regions of co-monotonicity and counter-monotonicity

clear;
clc;

xMin = 0;
xMax = 1;
num_noise = 30;                    % The number of different noise levels used
noise = 3;                         % A constant to determine the amount of noise
K = 25;

MVecResults = 100:100:1000;
num_noise_test_min = 0;
num_noise_test_max = 30;
noiseVec = num_noise_test_min:num_noise_test_max;
numMCSim = 100;

M = MVecResults(5);
l = 1;

x = rand(M,1)*(xMax-xMin)+xMin;

y1 = x + noise*(l/num_noise)*randn(M,1);
y2 = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
y3 = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
y4 = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
y5 = sin(16*pi*x) + noise*(l/num_noise)*randn(M,1);
y6 = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
y7=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
y8 = (x > 0.5) + noise*5*l/num_noise*randn(M,1);

[y1Regions, E1_copula, u1] = cosf_experiments(x,y1);
[y2Regions, E2_copula, u2] = cosf_experiments(x,y2);
[y3Regions, E3_copula, u3] = cosf_experiments(x,y3);
[y4Regions, E4_copula, u4] = cosf_experiments(x,y4);
% y5Regions = cosf_experiments(x,y5);
% y6Regions = cosf_experiments(x,y6);
% y7Regions = cosf_experiments(x,y7);
% y8Regions = cosf_experiments(x,y8);

celldisp(y1Regions)
celldisp(y2Regions)
% celldisp(y3Regions)
% celldisp(y4Regions)
% celldisp(y5Regions)
% celldisp(y6Regions)
% celldisp(y7Regions)
% celldisp(y8Regions)

sz = 10;

subplot(2,2,1); plot(u1, E1_copula, '-.b*', 'LineWidth', .5); hold on; scatter(pobs(x),pobs(y1),sz,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5); grid on; title('y=x')
subplot(2,2,2); plot(u2, E2_copula, '-.b*', 'LineWidth', .5); hold on; scatter(pobs(x),pobs(y2),sz,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5); grid on; title('y=x^2');
subplot(2,2,3); plot(u3, E3_copula, '-.b*', 'LineWidth', .5); hold on; scatter(pobs(x),pobs(y3),sz,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5); grid on; title('y=x^3')
subplot(2,2,4); plot(u4, E4_copula, '-.b*', 'LineWidth', .5); hold on; scatter(pobs(x),pobs(y4),sz,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5); grid on; title('y=sin(x)');
figtitle(sprintf('Noise=%0.02f',l))


