clear;
clc;

M = 1000;
x = rand(M,1);

num_noise = 30; noise = 3;
l = 5;

subplot(2,3,1);
y_noise = rand(M,1);
u = pobs(x); v = pobs(y_noise);
C1 = corner([u v]);
scatter(u,v); grid on;

subplot(2,3,2);
y_noise = x + noise*(l/num_noise)*randn(M,1); 
u = pobs(x); v = pobs(y_noise);
C2 = corner([u v]);
scatter(u,v); grid on;

subplot (2,3,3);
y_noise = 4*(x-0.5).^2 + noise*(l/num_noise)*randn(M,1);
u = pobs(x); v = pobs(y_noise);
C3 = corner([u v]);
scatter(u,v); grid on;

subplot(2,3,4);
y_noise = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
u = pobs(x); v = pobs(y_noise);
C4 = corner([u v]);
scatter(u,v); grid on;

subplot(2,3,5);
y_noise=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
u = pobs(x); v = pobs(y_noise);
C5 = corner([u v]);
scatter(u,v); grid on;

subplot(2,3,6);
y_noise = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
u = pobs(x); v = pobs(y_noise);
C6 = corner([u v]);
scatter(u,v); grid on;

%% view the signal as a time-series
clear;
clc;

M = 1000; divFactor = 10;
x = rand(M,1);

num_noise = 30; noise = 3;
l = 0;

subplot(2,3,1);
y_noise = rand(M,1); u = pobs(x); v = pobs(y_noise);
signal = v-u;
plot((1:M)/M,signal); grid on; title('indep'); hold on;
plot((1:M)/M,movmean(signal,M/divFactor),'r.-', 'LineWidth', 4);
[f,xi] = ksdensity(pdist(signal));
plot(xi,f,'g.-','LineWidth',3);

subplot(2,3,2);
y_noise = x + noise*(l/num_noise)*randn(M,1); u = pobs(x); v = pobs(y_noise);
signal = v-u;
plot((1:M)/M,signal); grid on; title('linear'); hold on;
plot((1:M)/M,movmean(signal,M/divFactor),'r.-', 'LineWidth', 4);
[f,xi] = ksdensity(pdist(signal));
plot(xi,f,'g.-','LineWidth',3);

subplot (2,3,3);
y_noise = 4*(x-0.5).^2 + noise*(l/num_noise)*randn(M,1); u = pobs(x); v = pobs(y_noise);
signal = v-u;
plot((1:M)/M,signal); grid on; title('quadratic'); hold on;
plot((1:M)/M,movmean(signal,M/divFactor),'r.-', 'LineWidth', 4);
[f,xi] = ksdensity(pdist(signal));
plot(xi,f,'g.-','LineWidth',3);

subplot(2,3,4);
y_noise = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1); u = pobs(x); v = pobs(y_noise);
signal = v-u;
plot((1:M)/M,signal); grid on; title('sinusoidal'); hold on;
plot((1:M)/M,movmean(signal,M/divFactor),'r.-', 'LineWidth', 4);
[f,xi] = ksdensity(pdist(signal));
plot(xi,f,'g.-','LineWidth',3);

subplot(2,3,5);
y_noise=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
u = pobs(x); v = pobs(y_noise);
signal = v-u;
plot((1:M)/M,signal); grid on; title('circular'); hold on;
plot((1:M)/M,movmean(signal,M/divFactor),'r.-', 'LineWidth', 4);
[f,xi] = ksdensity(pdist(signal));
plot(xi,f,'g.-','LineWidth',3);

subplot(2,3,6);
y_noise = x.^(1/4) + noise*(l/num_noise)*randn(M,1); 
u = pobs(x); v = pobs(y_noise);
signal = v-u;
plot((1:M)/M,signal); grid on; title('fourth-root'); hold on;
plot((1:M)/M,movmean(signal,M/divFactor),'r.-', 'LineWidth', 4);
[f,xi] = ksdensity(pdist(signal));
plot(xi,f,'g.-','LineWidth',3);

%% Plot the tau profile for different dependencies

clear;
clc;

M = 1000; divFactor = 10;
x = rand(M,1);

num_noise = 30; noise = 3;
l = 10;

subplot(2,3,1);
y_noise = rand(M,1); tau_profile = zeros(1,M-1);
kso_tau1 = taukl_s(x, y_noise);
for ii=2:M
    tau_profile(ii-1) = kso_tau1.consume(1);
end
plot(tau_profile); grid on; hold on; title('indep');

subplot(2,3,2);tau_profile = zeros(2,M-1);
y_no_noise = x;
y_noise = y_no_noise + noise*(l/num_noise)*randn(M,1); 
kso_tau1 = taukl_s(x, y_no_noise);
kso_tau2 = taukl_s(x, y_noise);
for ii=2:M
    tau_profile(1,ii-1) = kso_tau1.consume(1);
    tau_profile(2,ii-1) = kso_tau2.consume(1);
end
plot(1:M-1,tau_profile(1,:), 1:M-1, tau_profile(2,:)); 
grid on; hold on; title('linear');

subplot (2,3,3);tau_profile = zeros(2,M-1);
y_no_noise = 4*(x-0.5).^2;
y_noise = y_no_noise + noise*(l/num_noise)*randn(M,1); 
kso_tau1 = taukl_s(x, y_no_noise);
kso_tau2 = taukl_s(x, y_noise);
for ii=2:M
    tau_profile(1,ii-1) = kso_tau1.consume(1);
    tau_profile(2,ii-1) = kso_tau2.consume(1);
end
plot(1:M-1,tau_profile(1,:), 1:M-1, tau_profile(2,:)); 
grid on; hold on; title('quadratic');

subplot(2,3,4);tau_profile = zeros(2,M-1);
y_no_noise = sin(4*pi*x);
y_noise = y_no_noise + 2*noise*(l/num_noise)*randn(M,1);
kso_tau1 = taukl_s(x, y_no_noise);
kso_tau2 = taukl_s(x, y_noise);
for ii=2:M
    tau_profile(1,ii-1) = kso_tau1.consume(1);
    tau_profile(2,ii-1) = kso_tau2.consume(1);
end
plot(1:M-1,tau_profile(1,:), 1:M-1, tau_profile(2,:)); 
grid on; hold on; title('sinu');

subplot(2,3,5);tau_profile = zeros(2,M-1);
y_no_noise=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2));
y_noise=y_no_noise + noise/4*l/num_noise*randn(M,1);
kso_tau1 = taukl_s(x, y_no_noise);
kso_tau2 = taukl_s(x, y_noise);
for ii=2:M
    tau_profile(1,ii-1) = kso_tau1.consume(1);
    tau_profile(2,ii-1) = kso_tau2.consume(1);
end
plot(1:M-1,tau_profile(1,:), 1:M-1, tau_profile(2,:)); 
grid on; hold on; title('circle');

subplot(2,3,6);tau_profile = zeros(2,M-1);
y_no_noise = x.^(1/4);
y_noise = y_no_noise + noise*(l/num_noise)*randn(M,1);
kso_tau1 = taukl_s(x, y_no_noise);
kso_tau2 = taukl_s(x, y_noise);
for ii=2:M
    tau_profile(1,ii-1) = kso_tau1.consume(1);
    tau_profile(2,ii-1) = kso_tau2.consume(1);
end
plot(1:M-1,tau_profile(1,:), 1:M-1, tau_profile(2,:)); 
grid on; hold on; title('fourth-root');
