%% Search for non-monotonic relationships in the wine dataset

%% process the WINE quality
clear;
clc;

if(ispc)
    dataFolder = 'C:\\Users\\Kiran\\Documents\\data\\wine_quality';
elseif(ismac)
    dataFolder = '/Users/kiran/Documents/data/wine_quality';
else
    dataFolder = '/home/kiran/data/wine_quality';
end

redwine_dataFile = fullfile(dataFolder, 'winequality-red.csv');
whitewine_dataFile = fullfile(dataFolder, 'winequality-white.csv');

% ensure dataset is all numeric, and convert categorical data to numeric
redwine_data = importdata(redwine_dataFile,';');
whitewine_data = importdata(redwine_dataFile,';');

[R_redwine,RectanglesCell_redwine] = paircim_v4(redwine_data);
[R_whitewine,RectanglesCell_whitewine] = paircim_v4(whitewine_data);

% store the results
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\rwd\\wine.mat');
elseif(ismac)
    save('/Users/Kiran/ownCloud/PhD/sim_results/rwd/wine.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/rwd/wine.mat');
end

%% Analyze the WINE data
clear;
clc;
if(ispc)
    load('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\rwd\\wine.mat');
elseif(ismac)
    load('/Users/Kiran/ownCloud/PhD/sim_results/rwd/wine.mat');
else
    load('/home/kiran/ownCloud/PhD/sim_results/rwd/wine.mat');
end

%% process the CRIME data
clear;
clc;

if(ispc)
    dataFolder = 'C:\\Users\\Kiran\\Documents\\data\\crime';
elseif(ismac)
    dataFolder = '/Users/kiran/Documents/data/crime';
else
    dataFolder = '/home/kiran/data/crime';
end

crime_dataFile = fullfile(dataFolder, 'crime.csv');
crime_data = csvread(crime_dataFile,1,0);   % ignore the header