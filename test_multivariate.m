% Test the Multivariate CoS against some baselines

%% basic testing w/ Clayton and Gumbel copulas

clear;
clc;

dVec = [2,3]; % the dimensionalities for which we will generate data
numMCSims = 100;
M = 500;
copulasToTest = {'Frank', 'Gumbel', 'Clayton'};
numAlphaSteps = 10;
alphasCell = {linspace(0.01,30,numAlphaSteps), ...
              linspace(1,30,numAlphaSteps), ...
              linspace(0.01,30,numAlphaSteps)};
numMetricsToCompute = 5;  % we compute sRho1, sRho2, sRho3, sRho4, and CoS

resultsVec = zeros(length(dVec),length(copulasToTest),numAlphaSteps,numMetricsToCompute);

% Setup the Spearman's Rho's estimators
mult = 1;
ds_2d = ones(2,1); ds_3d = ones(3,1);
dsCell = {ds_2d, ds_3d};

s1_obj = ASpearman1_initialization(mult);
s2_obj = ASpearman2_initialization(mult);
s3_obj = ASpearman3_initialization(mult);
s4_obj = ASpearman4_initialization(mult);

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');

for dIdx=1:length(dVec)
    d = dVec(dIdx);
    ds = dsCell{dIdx};
    dispstat(sprintf('Computing for d=%d',d),'keepthis', 'timestamp');
    for copulaIdx=1:length(copulasToTest)
        copulaToTest = copulasToTest{copulaIdx};
        alphasToTest = alphasCell{copulaIdx};
        if(strcmpi(copulaToTest,'frank'))
            coprng = @frankcopularnd;
        elseif(strcmpi(copulaToTest,'clayton'))
            coprng = @claytoncopularnd;
        elseif(strcmpi(copulaToTest,'gumbel'))
            coprng = @gumbelcopularnd;
        end
        for alphaIdx=1:length(alphasToTest)
            alphaVal = alphasToTest(alphaIdx);
            cos_sum = 0; srho1_sum = 0; srho2_sum = 0; srho3_sum = 0; srho4_sum = 0;
            for mcSimNum=1:numMCSims
                U = coprng(M, d, alphaVal);
                UU = U';  % the format needed by the ITE toolbox :/
                
                % compute Spearman1
                sRho1 = ASpearman1_estimation(UU,ds,s1_obj);
                srho1_sum = srho1_sum + sRho1;
                
                % compute Spearman2
                sRho2 = ASpearman2_estimation(UU,ds,s2_obj);
                srho2_sum = srho2_sum + sRho2;
                
                % compute Spearman3
                sRho3 = ASpearman3_estimation(UU,ds,s3_obj);
                srho3_sum = srho3_sum + sRho3;
                
                % compute Spearman4
                sRho4 = ASpearman4_estimation(UU,ds,s4_obj);
                srho4_sum = srho4_sum + sRho4;
                
                % compute CoS
                if(d==2)
                    cos_sum = cos_sum + cosdv(U(:,1),U(:,2));
                elseif(d==3)
                    cos_sum = cos_sum + cos3d(U(:,1),U(:,2),U(:,3));
                end
            end
            % store average results
            resultsVec(dIdx,copulaIdx,alphaIdx,1) = srho1_sum/numMCSims;
            resultsVec(dIdx,copulaIdx,alphaIdx,2) = srho2_sum/numMCSims;
            resultsVec(dIdx,copulaIdx,alphaIdx,3) = srho3_sum/numMCSims;
            resultsVec(dIdx,copulaIdx,alphaIdx,4) = srho4_sum/numMCSims;
            resultsVec(dIdx,copulaIdx,alphaIdx,5) = cos_sum/numMCSims;
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
subplot(2,3,1); 
copulaIdx = 1;
dIdx = 1;
alphasToTest = alphasCell{copulaIdx};
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,1)),'o-.'); hold on;
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,2)),'+-.');
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,3)),'d-.');
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,4)),'v-.');
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,5)),'s-.');
title({copulasToTest{copulaIdx}, sprintf('D=%d',dVec(dIdx))});
xlabel('\alpha');
legend('\rho_s_1','\rho_s_2','\rho_s_3','\rho_s_4','CoS', 'location', 'southeast');
grid on;

subplot(2,3,2); 
copulaIdx = 2;
dIdx = 1;
alphasToTest = alphasCell{copulaIdx};
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,1)),'o-.'); hold on;
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,2)),'+-.');
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,3)),'d-.');
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,4)),'v-.');
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,5)),'s-.');
title({copulasToTest{copulaIdx}, sprintf('D=%d',dVec(dIdx))});
xlabel('\alpha');
grid on;

subplot(2,3,3); 
copulaIdx = 3;
dIdx = 1;
alphasToTest = alphasCell{copulaIdx};
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,1)),'o-.'); hold on;
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,2)),'+-.');
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,3)),'d-.');
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,4)),'v-.');
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,5)),'s-.');
title({copulasToTest{copulaIdx}, sprintf('D=%d',dVec(dIdx))});
xlabel('\alpha');
grid on;

subplot(2,3,4); 
copulaIdx = 1;
dIdx = 2;
alphasToTest = alphasCell{copulaIdx};
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,1)),'o-.'); hold on;
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,2)),'+-.');
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,3)),'d-.');
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,4)),'v-.');
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,5)),'s-.');
title({copulasToTest{copulaIdx}, sprintf('D=%d',dVec(dIdx))});
xlabel('\alpha');
grid on;

subplot(2,3,5); 
copulaIdx = 2;
dIdx = 2;
alphasToTest = alphasCell{copulaIdx};
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,1)),'o-.'); hold on;
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,2)),'+-.');
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,3)),'d-.');
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,4)),'v-.');
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,5)),'s-.');
title({copulasToTest{copulaIdx}, sprintf('D=%d',dVec(dIdx))});
xlabel('\alpha');
grid on;

subplot(2,3,6); 
copulaIdx = 3;
dIdx = 2;
alphasToTest = alphasCell{copulaIdx};
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,1)),'o-.'); hold on;
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,2)),'+-.');
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,3)),'d-.');
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,4)),'v-.');
plot(alphasToTest, squeeze(resultsVec(dIdx,copulaIdx,:,5)),'s-.');
title({copulasToTest{copulaIdx}, sprintf('D=%d',dVec(dIdx))});
xlabel('\alpha');
grid on;
