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
MVec = 100:100:1000;
num_noise_test_min = 1;
num_noise_test_max = 30;

M = MVec(1);
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

%% simulate the difference between computing monotonic and detecting the region exactly

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

numMCSim = 100;

for M=MVec
    % Test Empirical copula distance as a distance measure under various
    % amounts of noise
    num_noise_test_min = 1;
    num_noise_test_max = 30;
    noiseVec = num_noise_test_min:num_noise_test_max;
    
    linearDep     = zeros(3,length(noiseVec));
    quadraticDep  = zeros(3,length(noiseVec));
    cubicDep      = zeros(3,length(noiseVec));
    sinusoidalDep = zeros(3,length(noiseVec));
                    % 1 - using tau directly
                    % 2 - segmenting into regions
                    
    for l=noiseVec
        dispstat(sprintf('Computing for noise level=%d >> M=%d',l,M),'keepthis', 'timestamp');
        t1 = 0; c1 = 0; c11 = 0;
        t2 = 0; c2 = 0; c22 = 0;
        t3 = 0; c3 = 0; c33 = 0;
        t4 = 0; c4 = 0; c44 = 0;
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
            
            t1 = t1 + corr(x,y1,'type','kendall');
            c1 = c1 + corr(x,y1,'type','kendall');
            c11 = c11 + cim(x,y1);
            
            t2 = t2 + corr(x,y2,'type','kendall');
            r1 = inBoundedPts(U2(:,1),U2(:,2),0,0.5,0,1); 
            r2 = inBoundedPts(U2(:,1),U2(:,2),0.5,1,0,1);
            c2 = c2 + ( abs(corr(r1(:,1),r1(:,2),'type','kendall'))*0.5 + ...
                        abs(corr(r2(:,1),r2(:,2),'type','kendall'))*0.5 );
            c22 = c22 + cim(x,y2);
            
            t3 = t3 + corr(x,y3,'type','kendall');
            r1 = inBoundedPts(U3(:,1),U3(:,2),0,0.1138,0,1); 
            r2 = inBoundedPts(U3(:,1),U3(:,2),0.1138,0.5768,0,1);
            r3 = inBoundedPts(U3(:,1),U3(:,2),0.5768,1,0,1);
            c3 = c3 + ( abs(corr(r1(:,1),r1(:,2),'type','kendall'))*0.1138 + ...
                        abs(corr(r2(:,1),r2(:,2),'type','kendall'))*0.4368 + ...
                        abs(corr(r3(:,1),r3(:,2),'type','kendall'))*0.4232);
            c33 = c33 + cim(x,y3);
                    
            t4 = t4 + corr(x,y4,'type','kendall');
            r1 = inBoundedPts(U4(:,1),U4(:,2),0,0.1138,0,1); 
            r2 = inBoundedPts(U4(:,1),U4(:,2),0.1138,0.3413,0,1);
            r3 = inBoundedPts(U4(:,1),U4(:,2),0.3413,.6208,0,1);
            r4 = inBoundedPts(U4(:,1),U4(:,2),0.6208,.8483,0,1);
            r5 = inBoundedPts(U4(:,1),U4(:,2),0.8483,1,0,1);
            c4 = c4 + ( abs(corr(r1(:,1),r1(:,2),'type','kendall'))*0.1138 + ...
                        abs(corr(r2(:,1),r2(:,2),'type','kendall'))*0.2275 + ...
                        abs(corr(r3(:,1),r3(:,2),'type','kendall'))*0.2795 + ...
                        abs(corr(r4(:,1),r4(:,2),'type','kendall'))*0.2275 + ...
                        abs(corr(r5(:,1),r5(:,2),'type','kendall'))*0.1517);
            c44 = c44 + cim(x,y4);
        end
        linearDep(1,l) = t1/numMCSim; linearDep(2,l) = c1/numMCSim; linearDep(3,l) = c11/numMCSim;
        quadraticDep(1,l) = t2/numMCSim; quadraticDep(2,l) = c2/numMCSim; quadraticDep(3,l) = c22/numMCSim;
        cubicDep(1,l) = t3/numMCSim; cubicDep(2,l) = c3/numMCSim; cubicDep(3,l) = c33/numMCSim;
        sinusoidalDep(1,l) = t4/numMCSim; sinusoidalDep(2,l) = c4/numMCSim; sinusoidalDep(3,l) = c44/numMCSim;
    end

    % save the data
    if(ispc)
        save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\regionDetection_M_%d.mat', M));
    elseif(ismac)
        save(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/clustering/regionDetection_M_%d.mat', M));
    else
        save(sprintf('/home/kiran/ownCloud/PhD/sim_results/clustering/regionDetection_M_%d.mat', M));
    end
end

%% Plot the results
clear;
clc;

M = 100;

if(ispc)
    load(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\clustering\\regionDetection_M_%d.mat', M));
elseif(ismac)
    load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/clustering/regionDetection_M_%d.mat', M));
else
    load(sprintf('/home/kiran/ownCloud/PhD/sim_results/clustering/regionDetection_M_%d.mat', M));
end

subplot(2,2,1);
hh1 = plot(noiseVec,linearDep(1,:),'o-.', ...
     noiseVec,linearDep(2,:),'+-.', ...
     noiseVec,linearDep(3,:),'d-.');
xlabel('noise');
legend('SSE[M(u,v)]','SSE[W(u,v)]','EIG');
title(sprintf('Cubic Dependency M=%d',M));
hh1(1).LineWidth = 1.5; 
hh1(2).LineWidth = 1.5; 
hh1(3).LineWidth = 1.5; 

plot(linearDep(1,:)); hold on; plot(linearDep(2,:)); plot(linearDep(3,:)); grid on;
subplot(2,2,2);
plot(quadraticDep(1,:)); hold on; plot(quadraticDep(2,:)); plot(quadraticDep(3,:)); grid on;
subplot(2,2,3);
plot(cubicDep(1,:)); hold on; plot(cubicDep(2,:)); plot(cubicDep(3,:)); grid on;
subplot(2,2,4);
plot(sinusoidalDep(1,:)); hold on; plot(sinusoidalDep(2,:)); plot(sinusoidalDep(3,:)); grid on;