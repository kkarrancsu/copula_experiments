% Test the Multivariate CoS against some baselines

%% basic testing w/ Clayton and Gumbel copulas

clear;
clc;

dVec = [2,3]; % the dimensionalities for which we will generate data
numMCSims = 100;
M = 300;
copulasToTest = {'Clayton', 'Gumbel'};
alphasToTest = linspace(1.5,30,10);
numMetricsToCompute = 3;  % we compute Td, tau_d, and CoS

resultsVec = zeros(length(dVec),length(copulasToTest),length(alphasToTest),numMetricsToCompute);

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');

for dIdx=1:length(dVec)
    d = dVec(dIdx);
    dispstat(sprintf('Computing for d=%d',d),'keepthis', 'timestamp');
    for copulaIdx=1:length(copulasToTest)
        copulaToTest = copulasToTest{copulaIdx};
        if(strcmpi(copulaToTest,'clayton'))
            coprng = @claytoncopularnd;
        elseif(strcmpi(copulaToTest,'gumbel'))
            coprng = @gumbelcopularnd;
        end
        for alphaIdx=1:length(alphasToTest)
            alphaVal = alphasToTest(alphaIdx);
            T_d_sum = 0; tau_d_sum = 0; cos_sum = 0;
            for mcSimNum=1:numMCSims
                U = coprng(M, d, alphaVal);
                % compute T_d
                if(d==2)
                    T_d_sum = T_d_sum + corr(U(:,1),U(:,2),'type','kendall');
                elseif(d==3)
                    sumVal = corr(U(:,1),U(:,2),'type','kendall') + ...
                             corr(U(:,1),U(:,3),'type','kendall') + ...
                             corr(U(:,2),U(:,3),'type','kendall');
                    T_d_sum = T_d_sum + 1/6*sumVal;
                end
                
                % compute tau_d
                if(d==2)
                    tau_d_sum = tau_d_sum + corr(U(:,1),U(:,2),'type','kendall');
                elseif(d==3)
                    % a vectorized but naive way to compute! -- research
                    % better ways perhaps ...
                    sumVal = 0;
                    for ii=1:M
                        sumVal = sumVal + sum(sum(U(ii,:)<U(ii+1:end,:),2)==d);
                    end
                    tau_d_sum = tau_d_sum + (2^d/(M*(M-1))*sumVal-1)*(1/(2^(d-1)-1));
                end
                
                % compute CoS
                if(d==2)
                    cos_sum = cos_sum + cosdv(U(:,1),U(:,2));
                elseif(d==3)
                    cos_sum = cos_sum + cos3d(U(:,1),U(:,2),U(:,3));
                end
            end
            % store average results
            resultsVec(dIdx,copulaIdx,alphaIdx,1) = T_d_sum/numMCSims;
            resultsVec(dIdx,copulaIdx,alphaIdx,2) = tau_d_sum/numMCSims;
            resultsVec(dIdx,copulaIdx,alphaIdx,3) = cos_sum/numMCSims;
        end
    end
end

if(ispc)
    % TODO!
elseif(ismac)
    % TODO!
elseif(isunix)
    save('/home/kiran/ownCloud/PhD/sim_results/multivariate/mv_tests.mat');
end

%%
if(ispc)
    % TODO!
elseif(ismac)
    % TODO!
elseif(isunix)
    load('/home/kiran/ownCloud/PhD/sim_results/multivariate/mv_tests.mat');
end

figure;
subplot(2,2,1); 
copulaIdx = 1;
dIdx = 1;
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,1)),'o-.'); hold on;
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,2)),'+-.');
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,3)),'d-.');
title({copulasToTest{copulaIdx}, 'D=2'});
xlabel('\alpha');
legend('T_d','\tau_d','CoS');
grid on;

subplot(2,2,2); 
copulaIdx = 2;
dIdx = 1;
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,1)),'o-.'); hold on;
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,2)),'+-.');
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,3)),'d-.');
title({copulasToTest{copulaIdx}, 'D=2'});
xlabel('\alpha');
legend('T_d','\tau_d','CoS');
grid on;

subplot(2,2,3);
copulaIdx = 1;
dIdx = 2;
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,1)),'o-.'); hold on;
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,2)),'+-.');
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,3)),'d-.');
title({copulasToTest{copulaIdx}, 'D=3'});
xlabel('\alpha');
legend('T_d','\tau_d','CoS');
grid on;

subplot(2,2,4);
copulaIdx = 2;
dIdx = 2;
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,1)),'o-.'); hold on;
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,2)),'+-.');
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,3)),'d-.');
title({copulasToTest{copulaIdx}, 'D=3'});
xlabel('\alpha');
legend('T_d','\tau_d','CoS');
grid on;
