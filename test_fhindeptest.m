clear;
clc;
close all;
dbstop if error;

figure;

% Test the Frechet-Hoeffing M/W based independence test
M = 500;

x = rand(M,1);
y = x;
[dist_M1, dist_W1] = fhindeptest(x,y);
xi = linspace(0,1,100);
f_P = ksdensity(dist_M1,xi);
f_Q = ksdensity(dist_W1,xi);
kldiv_1 = kldivergence(f_P,f_Q,xi);
subplot(2,4,1); scatter(x,y); xlabel('x'); ylabel('y'); grid on;
subplot(2,4,5); plot(xi,f_P, xi, f_Q); grid on; title(sprintf('div=%0.02f',abs(kldiv_1)));

x = rand(M,1);
y = rand(M,1);
[dist_M2, dist_W2] = fhindeptest(x,y);
f_P = ksdensity(dist_M2,xi);
f_Q = ksdensity(dist_W2,xi);
kldiv_2 = kldivergence(f_P,f_Q,xi);
subplot(2,4,2); scatter(x,y); xlabel('x'); ylabel('y'); grid on;
subplot(2,4,6); plot(xi,f_P, xi, f_Q); grid on; title(sprintf('div=%0.02f',abs(kldiv_2)));

x = rand(M,1);
y = 4*(x-0.5).^2;
[dist_M3, dist_W3] = fhindeptest(x,y);
f_P = ksdensity(dist_M3,xi);
f_Q = ksdensity(dist_W3,xi);
kldiv_3 = kldivergence(f_P,f_Q,xi);
subplot(2,4,3); scatter(x,y); xlabel('x'); ylabel('y'); grid on;
subplot(2,4,7); plot(xi,f_P, xi, f_Q); grid on; title(sprintf('div=%0.02f',abs(kldiv_3)));

x = rand(M,1);
y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3);
[dist_M4, dist_W4] = fhindeptest(x,y);
f_P = ksdensity(dist_M4,xi);
f_Q = ksdensity(dist_W4,xi);
kldiv_4 = kldivergence(f_P,f_Q,xi);
subplot(2,4,4); scatter(x,y); xlabel('x'); ylabel('y'); grid on;
subplot(2,4,8); plot(xi,f_P, xi, f_Q); grid on; title(sprintf('div=%0.02f',abs(kldiv_4)));
