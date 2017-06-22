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

%% NMDR Data -  This data was retrieved from the Supplemental Material for the paper:
% Modeling non monotonic dose response relationships: Model evaluation
% and hormetic quantities exploration
% The supplementary material contained the data in Table's A.1, A.2, and
% was copy pasted into the Matlab vectors below.  The supplementary
% material which contains this data can be downloaded directly from:
% http://www.sciencedirect.com/science/article/pii/S0147651312004423#s0055
clear;
clc;

% OMIM data
c = [3.4240E-02 1.9520E-02 1.0960E-02 6.5100E-03 3.7700E-03 2.0500E-03 1.1600E-03 6.8482E-04 3.7665E-04 2.0545E-04 1.3697E-04 6.8482E-05];
e1 = [9.9444E-01 5.0751E-01 -1.9143E-01 -3.3849E-01 -3.1219E-01 -2.6279E-01 -2.0190E-01 -1.2020E-01 -6.2626E-02 -5.4037E-02 -4.8928E-02 3.2496E-02];
e2 = [9.9493E-01 5.5256E-01 -1.6175E-01 -2.9962E-01 -2.8679E-01 -2.3212E-01 -1.8031E-01 -1.0766E-01 -4.2618E-02 -2.3945E-02 2.0211E-02 7.5162E-02 ];
e3 = [9.9725E-01 5.3179E-01 -1.5113E-01 -3.4173E-01 -3.4713E-01 -2.9402E-01 -1.6422E-01 -1.6645E-01 -1.0187E-01 -8.3437E-02 -8.1991E-02 1.8784E-02 ];

c_data = [c c c];
e_data = [e1 e2 e3];
c_data = c_data'; e_data = e_data';

[z1, z2] = cim_v4(c_data,e_data)
subplot(2,2,1);
scatter(pobs(c_data),pobs(e_data)); title('OMIM'); grid on;

% HMIM data
c = [1.7711E-01 1.0627E-01 6.3760E-02 3.8970E-02 2.3020E-02 1.3810E-02 8.3200E-03 4.9600E-03 3.0100E-03 1.7700E-03 1.0600E-03 7.0846E-04 ];
e1 = [9.9826E-01 6.6800E-01 2.6039E-02 -3.2590E-01 -4.5497E-01 -4.1174E-01 -3.5150E-01 -2.9002E-01 -2.0483E-01 -1.2880E-01 -7.6196E-02 -5.9985E-02 ];
e2 = [9.9820E-01 7.2556E-01 4.5334E-02 -3.5571E-01 -4.8631E-01 -4.8037E-01 -4.0739E-01 -3.0845E-01 -2.4493E-01 -1.5849E-01 -1.0505E-01 -8.7289E-02 ];
e3 = [9.9813E-01 7.0787E-01 3.7437E-02 -3.6165E-01 -4.8802E-01 -4.7210E-01 -3.7956E-01 -3.0123E-01 -1.6900E-01 -1.2401E-01 -7.2647E-02 9.4150E-03 ];

c_data = [c c c];
e_data = [e1 e2 e3];
c_data = c_data'; e_data = e_data';

[z1, z2] = cim_v4(c_data,e_data)
subplot(2,2,2);
scatter(pobs(c_data),pobs(e_data)); title('HMIM'); grid on;

% acetonitrile data
c = [4.3563E+00 2.4806E+00 1.4087E+00 8.0388E-01 4.5936E-01 2.6030E-01 1.4546E-01 8.4220E-02 4.8230E-02 2.7560E-02 1.5310E-02 9.1900E-03 ];
e1 = [9.9997E-01 9.9997E-01 9.7583E-01 3.9872E-01 -6.1830E-02 -4.3318E-01 -4.0722E-01 -3.8752E-01 -1.0000E-01 -2.5230E-02 3.2800E-03 -1.8400E-03 ];
e2 = [9.9998E-01 9.9997E-01 9.7270E-01 3.2705E-01 -7.3760E-02 -2.6153E-01 -4.2807E-01 -2.9503E-01 -1.1850E-01 -2.6530E-02 4.4720E-02 3.5910E-02 ];
e3 = [9.9996E-01 9.9996E-01 9.1610E-01 3.1633E-01 -1.0094E-01 -2.6721E-01 -3.1688E-01 -2.6036E-01 -1.2523E-01 -2.5160E-02 1.3000E-02 2.5890E-02 ];
c_data = [c c c];
e_data = [e1 e2 e3];
c_data = c_data'; e_data = e_data';

[z1, z2] = cim_v4(c_data,e_data)
subplot(2,2,3);
scatter(pobs(c_data),pobs(e_data)); title('Acetonitrile'); grid on;

% Isopropyl Alcohol data
c = [3.0919E+00 1.7606E+00 9.9983E-01 5.7056E-01 3.2603E-01 1.8475E-01 1.0324E-01 5.9770E-02 3.4230E-02 1.9560E-02 1.0870E-02 6.5200E-03 ];
e1 = [9.9996E-01 9.9317E-01 6.8870E-01 1.9105E-01 -1.9568E-01 -4.8446E-01 -4.0812E-01 -3.0762E-01 -2.0580E-01 -2.0850E-02 6.8820E-02 4.5547E-04 ];
e2 = [9.9952E-01 9.9106E-01 5.8804E-01 2.1478E-01 -1.4135E-01 -2.9080E-01 -3.8510E-01 -2.3188E-01 -2.1087E-01 -1.4716E-01 -8.3310E-02 -8.3990E-02 ];
e3 = [9.9997E-01 9.8903E-01 6.5128E-01 3.4319E-01 -4.0040E-02 -2.9836E-01 -3.4702E-01 -1.4845E-01 -1.3920E-01 -8.5730E-02 2.4000E-03 -6.3610E-02 ];
c_data = [c c c];
e_data = [e1 e2 e3];
c_data = c_data'; e_data = e_data';

[z1, z2] = cim_v4(c_data,e_data)
subplot(2,2,4);
scatter(pobs(c_data),pobs(e_data)); title('Isopropyl Alcohol'); grid on;