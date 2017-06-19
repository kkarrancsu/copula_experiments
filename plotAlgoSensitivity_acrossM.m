function [] = plotAlgoSensitivity_acrossM(cimVersion,MVecToPlot)

switch cimVersion
    case 1
        fnameStr = 'cim';
    case 3
        fnameStr = 'cimv3';
    case 4
        fnameStr = 'cimv4';
    case 5
        fnameStr = 'cimv5';
    case 6
        fnameStr = 'cimv6';
    case 7
        fnameStr = 'cimv7';
    otherwise
end
CIMVECIDX = 3;
num_noise_test_min = 0;
num_noise_test_max = 30;
noiseVec = num_noise_test_min:num_noise_test_max;

figure;

linearDep = zeros(3,length(noiseVec));          % 1 - min
                                                % 2 - mean
                                                % 3 - max
quadraticDep = zeros(3,length(noiseVec));
cubicDep = zeros(3,length(noiseVec));
sinusoidalDep = zeros(3,length(noiseVec));
hiFreqSinDep = zeros(3,length(noiseVec));
fourthRootDep = zeros(3,length(noiseVec));
circleDep = zeros(3,length(noiseVec));
stepDep = zeros(3,length(noiseVec));
indep = zeros(3,length(noiseVec));

numLoopsRun = 1;
for MIdx=1:length(MVecToPlot)
    M = MVecToPlot(MIdx);
    % load the data
    if(ispc)
        load(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\%s_algoSensitivity_M_%d.mat', fnameStr, M));
    elseif(ismac)
        load(sprintf('/Users/Kiran/ownCloud/PhD/sim_results/independence/%s_algoSensitivity_M_%d.mat', fnameStr, M));
    else
        load(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/%s_algoSensitivity_M_%d.mat', fnameStr, M));
    end
    
    % collect all the data for each of the scanincrs we tested and store in
    % a matrix, for each dependency type, we plot the difference between
    % the minimum and the maximum
    for ii=1:length(scanincrsToTest)
        scanincrVal = scanincrsToTest(ii);
        rawData = algoSensitivityData{ii};  % algoSensitivityData is loaded from the file above
        linearData = rawData.linearDep(CIMVECIDX,:);
        quadraticData = rawData.quadraticDep(CIMVECIDX,:);
        cubicData = rawData.cubicDep(CIMVECIDX,:);
        sinusoidalData = rawData.sinusoidalDep(CIMVECIDX,:);
        hiFreqSinData = rawData.hiFreqSinDep(CIMVECIDX,:);
        fourthRootData = rawData.fourthRootDep(CIMVECIDX,:);
        circleData = rawData.circleDep(CIMVECIDX,:);
        stepData = rawData.stepDep(CIMVECIDX,:);
        indepData = rawData.indep(CIMVECIDX,:);
        
        if(MIdx==1)
            % seed the data
            if(ii==1)
                for jj=1:3
                    linearDep(jj,:) = linearData;
                    quadraticDep(jj,:) = quadraticData;
                    cubicDep(jj,:) = cubicData;
                    fourthRootDep(jj,:) = fourthRootData;
                    circleDep(jj,:) = circleData;
                    stepDep(jj,:) = stepData;
                    indep(jj,:) = indepData;
                end
            end
            if(scanincrVal<0.1138)
                for jj=1:3
                    sinusoidalDep(jj,:) = sinusoidalData;
                end
            end
            if(scanincrVal<0.029)
                for jj=1:3
                    hiFreqSinDep(jj,:) = hiFreqSinData;
                end
            end
        else
            for jj=1:length(noiseVec)
                % process linear
                if(linearData(jj)<linearDep(1,jj))
                    linearDep(1,jj) = linearData(jj);
                end
                linearDep(2,jj) = linearDep(2,jj) + linearData(jj);
                if(linearData(jj)>linearDep(3,jj))
                    linearDep(3,jj) = linearData(jj);
                end
                
                % process quadratic
                if(quadraticData(jj)<quadraticDep(1,jj))
                    quadraticDep(1,jj) = quadraticData(jj);
                end
                quadraticDep(2,jj) = quadraticDep(2,jj) + quadraticData(jj);
                if(quadraticData(jj)>quadraticDep(3,jj))
                    quadraticDep(3,jj) = quadraticData(jj);
                end
                
                % process cubic
                if(cubicData(jj)<cubicDep(1,jj))
                    cubicDep(1,jj) = cubicData(jj);
                end
                cubicDep(2,jj) = cubicDep(2,jj) + cubicData(jj);
                if(cubicData(jj)>cubicDep(3,jj))
                    cubicDep(3,jj) = cubicData(jj);
                end
                
                % process sinusoidal
                if(scanincrVal<0.1138)
                    if(sinusoidalData(jj)<sinusoidalDep(1,jj))
                        sinusoidalDep(1,jj) = sinusoidalData(jj);
                    end
                    sinusoidalDep(2,jj) = sinusoidalDep(2,jj) + sinusoidalData(jj);
                    if(sinusoidalData(jj)>sinusoidalDep(3,jj))
                        sinusoidalDep(3,jj) = sinusoidalData(jj);
                    end
                end
                
                % process hi-freq sine
                if(scanincrVal<0.029)
                    if(hiFreqSinData(jj)<hiFreqSinDep(1,jj))
                        hiFreqSinDep(1,jj) = hiFreqSinData(jj);
                    end
                    hiFreqSinDep(2,jj) = hiFreqSinDep(2,jj) + hiFreqSinData(jj);
                    if(hiFreqSinData(jj)>hiFreqSinDep(3,jj))
                        hiFreqSinDep(3,jj) = hiFreqSinData(jj);
                    end
                end
                
                % process fourth-root
                if(fourthRootData(jj)<fourthRootDep(1,jj))
                    fourthRootDep(1,jj) = fourthRootData(jj);
                end
                fourthRootDep(2,jj) = fourthRootDep(2,jj) + fourthRootData(jj);
                if(fourthRootData(jj)>fourthRootDep(3,jj))
                    fourthRootDep(3,jj) = fourthRootData(jj);
                end
                
                % process circular
                if(circleData(jj)<circleDep(1,jj))
                    circleDep(1,jj) = circleData(jj);
                end
                circleDep(2,jj) = circleDep(2,jj) + circleData(jj);
                if(circleData(jj)>circleDep(3,jj))
                    circleDep(3,jj) = circleData(jj);
                end
                
                % process step
                if(stepData(jj)<stepDep(1,jj))
                    stepDep(1,jj) = stepData(jj);
                end
                stepDep(2,jj) = stepDep(2,jj) + stepData(jj);
                if(stepData(jj)>stepDep(3,jj))
                    stepDep(3,jj) = stepData(jj);
                end
                
                % process indep
                if(indepData(jj)<indep(1,jj))
                    indep(1,jj) = indepData(jj);
                end
                indep(2,jj) = indep(2,jj) + indepData(jj);
                if(indepData(jj)>indep(3,jj))
                    indep(3,jj) = indepData(jj);
                end
            end
        end
        numLoopsRun = numLoopsRun + 1;
    end 
end

% scale the mean
linearDep(2,:) = linearDep(2,:)/numLoopsRun;
quadraticDep(2,:) = quadraticDep(2,:)/numLoopsRun;
cubicDep(2,:) = cubicDep(2,:)/numLoopsRun;
sinusoidalDep(2,:) = sinusoidalDep(2,:)/numLoopsRun;
hiFreqSinDep(2,:) = hiFreqSinDep(2,:)/numLoopsRun;
fourthRootDep(2,:) = fourthRootDep(2,:)/numLoopsRun;
circleDep(2,:) = circleDep(2,:)/numLoopsRun;
stepDep(2,:) = stepDep(2,:)/numLoopsRun;
indep(2,:) = indep(2,:)/numLoopsRun;

subplot(3,3,1);
jbfill(noiseVec,linearDep(3,:),linearDep(1,:));
grid on; xlabel('Noise'); title('Linear');

subplot(3,3,2);
jbfill(noiseVec,quadraticDep(3,:),quadraticDep(1,:));
grid on; xlabel('Noise');
title('Quadratic');

subplot(3,3,3);
jbfill(noiseVec,cubicDep(3,:),cubicDep(1,:));
grid on; xlabel('Noise');
title('Cubic');

subplot(3,3,4);
jbfill(noiseVec,sinusoidalDep(3,:),sinusoidalDep(1,:));
grid on; xlabel('Noise');
title('Sinusoidal');

subplot(3,3,5);
jbfill(noiseVec,hiFreqSinDep(3,:),hiFreqSinDep(1,:));
grid on; xlabel('Noise');
title('Hi-Freq Sin');

subplot(3,3,6);
jbfill(noiseVec,fourthRootDep(3,:),fourthRootDep(1,:));
grid on; xlabel('Noise');
title('Fourth-Root');

subplot(3,3,7);
jbfill(noiseVec,circleDep(3,:),circleDep(1,:));
grid on; xlabel('Noise');
title('Circular');

subplot(3,3,8);
jbfill(noiseVec,stepDep(3,:),stepDep(1,:));
grid on; xlabel('Noise');
title('Step');

subplot(3,3,9);
jbfill(noiseVec,indep(3,:),indep(1,:));
grid on; xlabel('Noise');
title('Independence');

end