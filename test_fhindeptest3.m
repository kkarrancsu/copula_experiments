clear;
clc;
close all;
dbstop if error;

% rng(1);

plot_pobs = 1;

xi = linspace(-sqrt(2)/2,sqrt(2)/2,100);

num_noise = 30;
noise = 3;
l = 10;

num_shuffle = 50;

% Test the Frechet-Hoeffing M/W based independence test
M = 500;

filterLen = M/20;
kernel_bw = 0.1;
kernel_type = 'epanechnikov';
a = 1;
b = 1/filterLen*ones(1,filterLen);

dep_types = {'indep', 'linear', 'quadratic', 'cubic', 'sinu', ...
             'gaussian_copula'};

x = rand(M,1); y = rand(M,1);
dist_indep = fhindeptest3(x,y,num_shuffle);

y = x+noise*(l/num_noise)*randn(M,1);
dist_linear = fhindeptest3(x,y,num_shuffle);

y = 4*(x-0.5).^2+noise*(l/num_noise)*randn(M,1);
dist_quadratic = fhindeptest3(x,y,num_shuffle);

y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
dist_cubic = fhindeptest3(x,y,num_shuffle);

y = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
dist_sinu = fhindeptest3(x,y,num_shuffle);

uv = copularnd('Gaussian', 0.2, M);
dist_gauss_copula = fhindeptest3(uv(:,1),uv(:,2),num_shuffle);

% numPlots = min(6,num_shuffle+1);
% for ii=1:numPlots
%     
%     dist_M = dist_indep(:,ii*2-1);
%     dist_W = dist_indep(:,ii*2);
%     f_M = ksdensity(dist_M,xi,'Bandwidth',kernel_bw,'Kernel',kernel_type);
%     f_W = ksdensity(dist_W,xi,'Bandwidth',kernel_bw,'Kernel',kernel_type);
%     
%     subplot(4,numPlots,ii);
% %     plot(1:M,filter(b,a,dist_M),1:M,filter(b,a,dist_W)); grid on; ylabel('distance');
%     plot(xi,f_M,xi,f_W); grid on;
%     title(sprintf('Indep. / %0.02f / %0.02f / %0.02f', ...
%                emd(f_M,f_W,xi), ...
%                abs(skewness(dist_M)), abs(skewness(dist_W)) ));
%     
%     dist_M = dist_linear(:,ii*2-1);
%     dist_W = dist_linear(:,ii*2);
%     f_M = ksdensity(dist_M,xi,'Bandwidth',kernel_bw,'Kernel',kernel_type);
%     f_W = ksdensity(dist_W,xi,'Bandwidth',kernel_bw,'Kernel',kernel_type);
%     subplot(4,numPlots,ii+numPlots);
%     plot(xi,f_M,xi,f_W); grid on;
%     title(sprintf('Linear / %0.02f / %0.02f / %0.02f', ...
%                emd(f_M,f_W,xi), ...
%                abs(skewness(dist_M)), abs(skewness(dist_W)) ));
%     
%     dist_M = dist_quadratic(:,ii*2-1);
%     dist_W = dist_quadratic(:,ii*2);
%     f_M = ksdensity(dist_M,xi,'Bandwidth',kernel_bw,'Kernel',kernel_type);
%     f_W = ksdensity(dist_W,xi,'Bandwidth',kernel_bw,'Kernel',kernel_type);
%     subplot(4,numPlots,ii+2*numPlots);
%     plot(xi,f_M,xi,f_W); grid on;
%     title(sprintf('Quadratic / %0.02f / %0.02f / %0.02f', ...
%                emd(f_M,f_W,xi), ...
%                abs(skewness(dist_M)), abs(skewness(dist_W)) ));
%     
%     dist_M = dist_cubic(:,ii*2-1);
%     dist_W = dist_cubic(:,ii*2);
%     f_M = ksdensity(dist_M,xi,'Bandwidth',kernel_bw,'Kernel',kernel_type);
%     f_W = ksdensity(dist_W,xi,'Bandwidth',kernel_bw,'Kernel',kernel_type);
%     subplot(4,numPlots,ii+3*numPlots);
%     plot(xi,f_M,xi,f_W); grid on;
%     title(sprintf('Cubic / %0.02f / %0.02f / %0.02f', ...
%                emd(f_M,f_W,xi), ...
%                abs(skewness(dist_M)), abs(skewness(dist_W)) )); 
% end

% Run KS-Tests of distribution of distances between the original and
% shuffled versions for each dependency type
for depType=dep_types
    if(strcmpi(depType,'indep'))
        distMat = dist_indep;
    elseif(strcmpi(depType, 'linear'))
        distMat = dist_linear;
    elseif(strcmpi(depType, 'quadratic'))
        distMat = dist_quadratic;
    elseif(strcmpi(depType, 'cubic'))
        distMat = dist_cubic;
    elseif(strcmpi(depType, 'sinu'))
        distMat = dist_sinu;
    elseif(strcmpi(depType, 'gaussian_copula'))
        distMat = dist_gauss_copula;
    end
    
    ref_dist_M = distMat(:,1);
    ref_dist_W = distMat(:,2);
    ref_dist_V = distMat(:,3);
    ref_dist_H = distMat(:,4);
    
    ksResults_M = zeros(1,num_shuffle);
    ksResults_W = zeros(1,num_shuffle);
    ksResults_V = zeros(1,num_shuffle);
    ksResults_H = zeros(1,num_shuffle);
    for ii=2:num_shuffle+1
        dist_M = distMat(:,ii*4-3);
        dist_W = distMat(:,ii*4-2);
        dist_V = distMat(:,ii*4-1);
        dist_H = distMat(:,ii*4-0);
        
        ksResults_M(ii-1) = kstest2(ref_dist_M,dist_M);
        ksResults_W(ii-1) = kstest2(ref_dist_W,dist_W);
        ksResults_V(ii-1) = kstest2(ref_dist_V,dist_V);
        ksResults_H(ii-1) = kstest2(ref_dist_H,dist_H);
    end
    
    fprintf('*** %s ***\n', char(depType));
    fprintf('M >> '); fprintf('%d ', any(ksResults_M==1)); fprintf('\n');
    fprintf('W >> '); fprintf('%d ', any(ksResults_W==1)); fprintf('\n');
    fprintf('V >> '); fprintf('%d ', any(ksResults_V==1)); fprintf('\n');
    fprintf('H >> '); fprintf('%d ', any(ksResults_H==1)); fprintf('\n');
    fprintf('***********************\n\n');
    
end