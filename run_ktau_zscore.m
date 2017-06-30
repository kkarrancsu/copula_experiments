%% Analysis of the Z-Score for Kendalls' Tau, and how it relates to the CIM algorithm

clear;
clc;

tauVec = 0.1:0.1:1;
nVec = 100:100:1000;

resVec = zeros(length(tauVec),length(nVec));
legendCell = cell(1,length(tauVec));
for tauIdx=1:length(tauVec)
    tauVal = tauVec(tauIdx);
    resVec(tauIdx,:) = ktau_zscore(tauVal,nVec);
    
    plot(nVec, 1./resVec(tauIdx,:)); hold on;
    legendCell{tauIdx} = sprintf('%0.2f',tauVal);
end

grid on;
ylabel('1/ZSC');
xlabel('M');
legend(legendCell);

%% experiments for counting the # of concordant and discordant samples
clear;
clc;

M = 500;
x = rand(M,1);
y = 4*(x-0.5).^2;

[x,I] = sort(x); y = y(I);

for ii=1:M-1
    for jj=1:ii
        valVec = sign(x(ii)-x(ii+1:M)).*sign(y(ii)-y(ii+1:M));
        nc = sum(valVec>0); nd = sum(valVec<0);
    end
    
    if(ii>(M/2))
        tauActual = corr(x(1:ii+1),y(1:ii+1),'type','kendall');
        tauCompute1 = (sum(ncnd_vec(1:ii,1))-sum(ncnd_vec(1:ii,2)))/nchoosek(ii,2);
        fprintf('[%d] -- actual=%0.02f compute1=%0.02f\n', ii, tauActual, tauCompute1);
    end
    
end