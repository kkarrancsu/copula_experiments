%% Test how tau-kl is affected by skewed distributions for continuous vs. binary
% Note, we only care about binary here, so we generate bernoulli
% distribution, but in the future, we also want to simulate non-binary
% skewed output variables

clear;
clc;
dbstop if error;

scenarios = {'left-skew','no-skew','right-skew'};
tauVec = linspace(0.01,0.99,15);                    
copulas = {'Gaussian','Frank','Gumbel','Clayton'};
M = 500;
numMCSims = 100;

% manually generate a left-skewed and right-skewed data, from which we
% construct an empirical cdf
leftSkewData = pearsrnd(0,1,-1,3,5000,1);
rightSkewData = pearsrnd(0,1,1,3,5000,1);

[fLeftSkew,xiLeftSkew] = emppdf(leftSkewData,0);
FLeftSkew = empcdf(xiLeftSkew,0);
leftSkewContinuousDistInfo = rvEmpiricalInfo(xiLeftSkew,fLeftSkew,FLeftSkew,0);
[fRightSkew,xiRightSkew] = emppdf(rightSkewData,0);
FRightSkew = empcdf(xiRightSkew,0);
rightSkewContinuousDistInfo = rvEmpiricalInfo(xiRightSkew,fRightSkew,FRightSkew,0);

resVecTauB = zeros(numMCSims,length(copulas),length(tauVec),length(scenarios),length(scenarios));
resVecTauN = zeros(numMCSims,length(copulas),length(tauVec),length(scenarios),length(scenarios));
resVecTauKL = zeros(numMCSims,length(copulas),length(tauVec),length(scenarios),length(scenarios));

sampleDataMat = zeros(M,2,length(copulas),length(tauVec),length(scenarios),length(scenarios));

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');

aa = 1; 
for continuousDistScenario=scenarios
    bb = 1; 
    for discreteDistScenario=scenarios
        cc = 1; 
        for tau=tauVec
            dd = 1; 
            for cop=copulas
                dispstat(sprintf('Simulating for {%s,%s} tau=%0.02f >> %s',cell2str(continuousDistScenario),cell2str(discreteDistScenario),tau,cell2str(cop)), ...
                         'keepthis', 'timestamp');
        
                iTau = copulaparam(cop,tau);
                parfor mcSimNum=1:numMCSims
                    % generate U
                    U = copularnd(cop,iTau,M);
                    
                    % generate F_X
                    if(strcmpi('left-skew',continuousDistScenario))
                        X = zeros(length(U),1);
                        % TODO: vectorize
                        for ii=1:length(U)
                            X(ii) = leftSkewContinuousDistInfo.icdf(U(ii,1));
                        end
                    elseif(strcmpi('no-skew',continuousDistScenario))
                        distObj = makedist('Normal');
                        X = icdf(distObj,U(:,1));
                    elseif(strcmpi('right-skew',continuousDistScenario))
                        X = zeros(length(U),1);
                        % TODO: vectorize
                        for ii=1:length(U)
                            X(ii) = rightSkewContinuousDistInfo.icdf(U(ii,1));
                        end
                    end
                    
                    % generate F_Y
                    if(strcmpi('left-skew',discreteDistScenario))
                        distObj = makedist('Multinomial','probabilities',[0.1,0.9]);
                    elseif(strcmpi('no-skew',discreteDistScenario))
                        distObj = makedist('Multinomial','probabilities',[0.5,0.5]);
                    elseif(strcmpi('right-skew',discreteDistScenario))
                        distObj = makedist('Multinomial','probabilities',[0.9,0.1]);
                    end
                    Y = icdf(distObj,U(:,2));

                    % compute tau, tau_kl, tau_N and record
                    resVecTauN(mcSimNum,dd,cc,bb,aa) = corr(X,Y,'type','kendall');
                    resVecTauB(mcSimNum,dd,cc,bb,aa) = ktaub([X Y], 0.05, 0);
                    resVecTauKL(mcSimNum,dd,cc,bb,aa) = taukl_cc_mex_interface(X,Y,0,1,0);
%                     [u,v] = pobs_sorted_cc(X,Y); 
%                     resVecTauKL(mcSimNum,dd,cc,bb,aa) = taukl_cc_mex(u,v,int32(0),int32(1),int32(0));
                end                
                
                %%%%%%%%%%%%%%%%%%%% MESSY CODE !!!!!! %%%%%%%%%%%%%%%%%%%%
                % COPIED FROM ABOVE, b/c parfor suks sometimes :/
                % generate U
                U = copularnd(cop,iTau,M);

                % generate F_X
                if(strcmpi('left-skew',continuousDistScenario))
                    X = zeros(length(U),1);
                    % TODO: vectorize
                    for ii=1:length(U)
                        X(ii) = leftSkewContinuousDistInfo.icdf(U(ii,1));
                    end
                elseif(strcmpi('no-skew',continuousDistScenario))
                    distObj = makedist('Normal');
                    X = icdf(distObj,U(:,1));
                elseif(strcmpi('right-skew',continuousDistScenario))
                    X = zeros(length(U),1);
                    % TODO: vectorize
                    for ii=1:length(U)
                        X(ii) = rightSkewContinuousDistInfo.icdf(U(ii,1));
                    end
                end

                % generate F_Y
                numIndepTrials = 1;
                if(strcmpi('left-skew',discreteDistScenario))
                    distObj = makedist('Multinomial','probabilities',[0.1,0.9]);
                elseif(strcmpi('no-skew',discreteDistScenario))
                    distObj = makedist('Multinomial','probabilities',[0.5,0.5]);
                elseif(strcmpi('right-skew',discreteDistScenario))
                    distObj = makedist('Multinomial','probabilities',[0.9,0.1]);
                end
                Y = icdf(distObj,U(:,2));
                sampleDataMat(:,1,dd,cc,bb,aa) = X;
                sampleDataMat(:,2,dd,cc,bb,aa) = Y;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                dd = dd + 1;
            end
            cc = cc + 1;
        end
        bb = bb + 1;
    end
    aa = aa + 1;
end

% save the results
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\skewed_data\\binary_output_class.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/skewed_data/binary_output_class.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/skewed_data/binary_output_class.mat');
end

%% plot the bias

clear;
clc;

if(ispc)
    load('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\skewed_data\\binary_output_class.mat');
elseif(ismac)
    load('/Users/Kiran/ownCloud/PhD/sim_results/skewed_data/binary_output_class.mat');
else
    load('/home/kiran/ownCloud/PhD/sim_results/skewed_data/binary_output_class.mat');
end


sampleData1_subplotIdxVec = [1,3,5,   13,15,17, 25,27,29];
sampleData2_subplotIdxVec = [2,4,6,   14,16,18, 26,28,30];
sampleData3_subplotIdxVec = [7,9,11,  19,21,23, 31,33,35];
results_subplotIdxVec     = [8,10,12, 20,22,24, 32,34,36];

for ii=1:length(copulas)
% for ii=1:1
    copToVis = copulas{ii};
    % find the appropriate indices
    dd = find(contains(copulas,copToVis));

    f = figure;
    subplotIdx = 1;
    for aa=1:length(scenarios)
        for bb=1:length(scenarios)
            tauNVec = mean(squeeze(resVecTauN(:,dd,:,bb,aa)));
            tauBVec = mean(squeeze(resVecTauB(:,dd,:,bb,aa)));
            tauKLVec = mean(squeeze(resVecTauKL(:,dd,:,bb,aa)));
            
            % plot the first sample data
            hh = subplot(6,6,sampleData1_subplotIdxVec(subplotIdx));
            tauIdx = 1;
            X = sampleDataMat(:,1,dd,tauIdx,bb,aa);
            Y = sampleDataMat(:,2,dd,tauIdx,bb,aa);
            outCell = separateOverlaps_binary(X,Y);
            nonOvlpPts = outCell{1}; ovlpPts = outCell{2};
            scatter(nonOvlpPts(:,1),nonOvlpPts(:,2),'+'); 
            hold on;
            scatter(ovlpPts(:,1),ovlpPts(:,2),'filled','d');
            set(hh,'XTick',[]); set(hh,'YTick',[]); 
            axis([min(X),max(X),0.8,2.2]);
%             xText = min(X)+0.1; yText = (min(Y)+max(Y))/2; 
%             txt = sprintf('%0.02f >> %0.02f | %0.02f', tauVec(tauIdx), length(nonOvlpPts)/M, length(ovlpPts)/M); 
%             text(xText,yText,txt);
%             xText = min(X)+0.1; yText = (min(Y)+max(Y))/2.4; 
%             txt = sprintf('\t %0.02f | %0.02f', skewness(X), skewness(Y)); 
%             text(xText,yText,txt);
%             taukl(X,Y);
            
            % plot the second sample data
            hh = subplot(6,6,sampleData2_subplotIdxVec(subplotIdx));
            tauIdx = floor(length(tauVec)/2);
            X = sampleDataMat(:,1,dd,tauIdx,bb,aa);
            Y = sampleDataMat(:,2,dd,tauIdx,bb,aa);
            outCell = separateOverlaps_binary(X,Y);
            nonOvlpPts = outCell{1}; ovlpPts = outCell{2};
            scatter(nonOvlpPts(:,1),nonOvlpPts(:,2),'+'); 
            hold on;
            scatter(ovlpPts(:,1),ovlpPts(:,2),'filled','d');
            set(hh,'XTick',[]); set(hh,'YTick',[]); 
            axis([min(X),max(X),0.8,2.2]);
%             xText = min(X)+0.1; yText = (min(Y)+max(Y))/2; 
%             txt = sprintf('%0.02f >> %0.02f | %0.02f', tauVec(tauIdx), length(nonOvlpPts)/M, length(ovlpPts)/M); 
%             text(xText,yText,txt);
%             xText = min(X)+0.1; yText = (min(Y)+max(Y))/2.4; 
%             txt = sprintf('\t %0.02f | %0.02f', skewness(X), skewness(Y)); 
%             text(xText,yText,txt);
%             taukl(X,Y);
            
            % plot the third sample data
            hh = subplot(6,6,sampleData3_subplotIdxVec(subplotIdx));
            tauIdx = length(tauVec);
            X = sampleDataMat(:,1,dd,tauIdx,bb,aa);
            Y = sampleDataMat(:,2,dd,tauIdx,bb,aa);
            outCell = separateOverlaps_binary(X,Y);
            nonOvlpPts = outCell{1}; ovlpPts = outCell{2};
            scatter(nonOvlpPts(:,1),nonOvlpPts(:,2),'+'); 
            hold on;
            scatter(ovlpPts(:,1),ovlpPts(:,2),'filled','d');
            set(hh,'XTick',[]); set(hh,'YTick',[]); 
            axis([min(X),max(X),0.8,2.2]);
%             xText = min(X)+0.1; yText = (min(Y)+max(Y))/2; 
%             txt = sprintf('%0.02f >> %0.02f | %0.02f', tauVec(tauIdx), length(nonOvlpPts)/M, length(ovlpPts)/M); 
%             text(xText,yText,txt);
%             xText = min(X)+0.1; yText = (min(Y)+max(Y))/2.4; 
%             txt = sprintf('\t %0.02f | %0.02f', skewness(X), skewness(Y)); 
%             text(xText,yText,txt);
%             taukl(X,Y);
            
            % plot the results
            subplot(6,6,results_subplotIdxVec(subplotIdx));
            plot(tauVec,tauBVec,tauVec,tauNVec,tauVec,tauKLVec);
            xlabel('\tau');
            title(sprintf('%s | %s', scenarios{bb}, scenarios{aa}));
            if(subplotIdx==9)
                legend('\tau_b','\tau_N','\tau_{KL}', 'location', 'northwest');
            end
            grid on;
            
            subplotIdx = subplotIdx + 1;
        end
    end
    tightfig;
    set(f, 'Name', copToVis);
end
