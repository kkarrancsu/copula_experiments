%% Test the mRMR algorithm on various estimators of MI for different datasets

clear;
clc;
dbstop if error;

% setup the estimators of MI
minScanIncr = 0.015625;

functionHandlesCell = {@cim_v2_hybrid_mi;};

functionArgsCell    = {{}};
fNames = {'cim'};

datasets = {'dexter','dorothea','arcene','gisette','madelon'};

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');

for dIdx=1:length(datasets)
    dataset = datasets{dIdx};
    if(ispc)
        folder = 'C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\feature_select_challenge';
    elseif(ismac)
        folder = '/Users/Kiran/ownCloud/PhD/sim_results/feature_select_challenge';
    else
        folder = '/home/kiran/ownCloud/PhD/sim_results/feature_select_challenge';
    end
    dispstat(sprintf('Processing %s',dataset),'keepthis', 'timestamp');

    load(fullfile(folder,dataset,'data.mat'));
    X = double(X_train);
    y = double(y_train);

    numFeaturesToSelect = 50;
    for ii=1:length(fNames)
        fs_outputFname = strcat(dataset,'_fs_',fNames{ii},'.mat');
        fOut = fullfile(folder,dataset,fs_outputFname);
        dispstat(sprintf('\t> Processing %s',fNames{ii}),'keepthis', 'timestamp');
        % if file exists, don't re-do it!
        if(~exist(fOut,'file'))
            tic;
            featureVec = mrmr_mid(X, y, numFeaturesToSelect, functionHandlesCell{ii}, functionArgsCell{ii});
            elapsedTime = toc;
            save(fOut,'featureVec','elapsedTime');
        end
    end
end