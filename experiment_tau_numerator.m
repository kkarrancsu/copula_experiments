%% experiment with different ways of computing the numerator of kendall's tau
% trade speed for memory efficiency

clear;
clc;

u = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]';
v = u;

uu = repmat(u,1,length(u));
vv = repmat(v,1,length(v));
% uuMask = tril(uu,0);

n = size(uu,1);
opMat = -1*ones(size(uu));
% opMask = tril(opMat,0);
opMat(1:n+1:n*n) = 1;

U = uu.*opMat;
V = vv.*opMat;
% X = tril(X)
cmpMatrixU = repmat(diag(U)',length(u),1);
cmpMatrixV = repmat(diag(V)',length(v),1);
% cmpMatrix = tril(cmpMatrix)
sumMatU = U+cmpMatrixU;
sgnMatU = tril(sign(sumMatU),-1);
sumMatV = V+cmpMatrixV;
sgnMatV = tril(sign(sumMatV),-1);

K = sgnMatU.*sgnMatV
sum(sum(K))

%% prototype gpu version
clear;
clc;

x = rand(5000,1); y = rand(5000,1);
xGpu = gpuArray(x); yGpu = gpuArray(y);
K = kendallsTauNumer2_GPU(xGpu,yGpu)

%% Profile our different ways of computing K (without counting memory allocation/deallocation costs)

clear;
clc;

numMCSims = 500;
M = 6100;

k1Times = zeros(1,numMCSims);
k2Times = zeros(1,numMCSims);
k2GpuTimes = zeros(1,numMCSims);

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');

for mcSimNum=1:numMCSims
    dispstat(sprintf('%d/%d',mcSimNum, numMCSims),'timestamp');
        
    x = rand(M,1); y = rand(M,1);
    
    tic;
    K1 = kendallsTauNumer1(x,y);
    k1Times(mcSimNum) = toc;
    
    tic;
    K2 = kendallsTauNumer2(x,y);
    k2Times(mcSimNum) = toc;
    
    % try the GPU version
    tic;
    xGpu = gpuArray(x); yGpu = gpuArray(y);
    K2GPU = kendallsTauNumer2_GPU(xGpu,yGpu);
    k2GpuTimes(mcSimNum) = toc;
    
    assert(K1==K2); assert(K1==K2GPU);
end

subplot(1,3,1); histogram(k1Times); grid on; title('K1');
subplot(1,3,2); histogram(k2Times); grid on; title('K2');
subplot(1,3,3); histogram(k2GpuTimes); grid on; title('K2[GPU]');
suplabel(sprintf('M=%d',M),'t');
suplabel('processing time (sec)', 'x');