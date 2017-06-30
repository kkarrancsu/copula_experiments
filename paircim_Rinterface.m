function [] = paircim_Rinterface(fname)
% Function which is the R-Interface for paircim, to run the MRNET
% experiments

load(fname);

% this is the code required if the incoming data is a dataframe
% % data is stored in a variable called 'dat'
% keys = fieldnames(dat);
% X = zeros(length(getfield(dat,keys{1})),length(keys));
% for keyIdx=1:length(keys)
%     key = keys{keyIdx};
%     % create a matrix of data that we can then deploy the paircim against
%     X(:,keyIdx) = getfield(dat,key);
% end
% this is the code required if the incoming data from R is a matrix
X = dat;

% compute the pairwise-cim estimates
R = paircim_v4( X );

fnameOut = strcat(fname,'.matlab');
dlmwrite(fnameOut,R)

end