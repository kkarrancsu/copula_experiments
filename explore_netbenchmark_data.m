%% explore the netbenchmark data
clear;
clc;

dataRepo = '/Users/kiran/data/netbenchmark/inputs/';
f = fullfile(dataRepo,'syntren300_1.mat');
load(f);

n = 10;
x = data(:,1);

for ii=2:n
    y = data(:,ii);
    cimVal = cim_v4_cc_mex(x,y);
    fprintf('cim=%0.02f\n',cimVal);
end

matlabOutputRepo = '/Users/kiran/data/netbenchmark/matlab_outputs/';
f = fullfile(matlabOutputRepo,'syntren300_1_output.mat');
load(f);

%% try to compute pairwise correlation w/ tau and see what R says about the inferred network
clear;
clc;

syntrenIdx = 1;

dataRepo = '/Users/kiran/data/netbenchmark/inputs/';
f = fullfile(dataRepo,sprintf('syntren300_%d.mat',syntrenIdx));
load(f);

n = size(data,2);
R = corr(data, 'type', 'kendall');
% R(1:n+1:n*n) = 0;

matlabOutputRepo = '/Users/kiran/data/netbenchmark/matlab_outputs/';
f = fullfile(matlabOutputRepo,sprintf('syntren300_%d_tauoutput.mat',syntrenIdx));
save(f,'R');

%% 
clear;
clc;

clear;
clc;

syntrenIdx = 3;

dataRepo = '/Users/kiran/data/netbenchmark/inputs/';
f = fullfile(dataRepo,sprintf('syntren300_%d.mat',syntrenIdx));
load(f);  % creates a variable called "data" we can use for access and debugging

matlabOutputRepo = '/Users/kiran/data/netbenchmark/matlab_outputs/';
load(fullfile(matlabOutputRepo,sprintf('syntren300_%d_tauoutput.mat',syntrenIdx)));
tauMat = R;
n = size(tauMat,2);
tauMat(1:n+1:n*n) = 0;
tauMat = tauMat.^2;

load(fullfile(matlabOutputRepo,sprintf('syntren300_%d_output.mat',syntrenIdx)));
cimMat = R.^2;

% identify where significant differneces are
diffMat = tauMat-cimMat;
subplot(1,2,1);
imagesc(diffMat(1:50,1:50))
colorbar
subplot(1,2,2);
plot(min(diffMat));