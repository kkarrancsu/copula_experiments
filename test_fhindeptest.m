clear;
clc;
close all;
dbstop if error;

plot_pobs = 1;

xi = linspace(-sqrt(2)/2,sqrt(2)/2,100);

num_noise = 30;
noise = 3;
l = 10;

% Test the Frechet-Hoeffing M/W based independence test
M = 500;

filterLen = M/20;
a = 1;
b = 1/filterLen*ones(1,filterLen);

x = rand(M,1); y = rand(M,1);
u = pobs(x); v = pobs(y);
[dist_M1, dist_W1, dist_V1, dist_H1] = fhindeptest2(x,y);
f_M = ksdensity(dist_M1,xi);
f_W = ksdensity(dist_W1,xi);
f_V = ksdensity(dist_V1,xi);
f_H = ksdensity(dist_H1,xi);
subplot(3,6,1); 
if(plot_pobs)
    scatter(u,v); xlabel('u'); ylabel('v');
else
    scatter(x,y); xlabel('x'); ylabel('y');
end
grid on;
subplot(3,6,7); plot(1:M,filter(b,a,dist_M1),1:M,filter(b,a,dist_W1),1:M,filter(b,a,dist_V1),1:M,filter(b,a,dist_H1));
grid on; title('x indep. y'); ylabel('distance');
legend('M','W');
subplot(3,6,13); plot(xi,f_M,xi,f_W); grid on;
title(sprintf('%0.02f / %0.02f', ...
               kldivergence(f_M,f_W,xi), ...
               min(abs(skewness(dist_M1)), abs(skewness(dist_W1))) ));

x = rand(M,1); y = x+noise*(l/num_noise)*randn(M,1);
u = pobs(x); v = pobs(y);
[dist_M2, dist_W2, dist_V2, dist_H2] = fhindeptest2(x,y);
f_M = ksdensity(dist_M2,xi);
f_W = ksdensity(dist_W2,xi);
f_V = ksdensity(dist_V2,xi);
f_H = ksdensity(dist_H2,xi);
subplot(3,6,2); 
if(plot_pobs)
    scatter(u,v); xlabel('u'); ylabel('v');
else
    scatter(x,y); xlabel('x'); ylabel('y');
end
grid on;
subplot(3,6,8); plot(1:M,filter(b,a,dist_M2),1:M,filter(b,a,dist_W2),1:M,filter(b,a,dist_V2),1:M,filter(b,a,dist_H2));
grid on; title('y=x'); ylabel('distance');
subplot(3,6,14); plot(xi,f_M,xi,f_W); grid on;
title(sprintf('%0.02f / %0.02f', ...
               kldivergence(f_M,f_W,xi), ...
               min(abs(skewness(dist_M2)), abs(skewness(dist_W2))) ));
           
x = rand(M,1); y = 4*(x-0.5).^2+noise*(l/num_noise)*randn(M,1);
u = pobs(x); v = pobs(y);
[dist_M3, dist_W3, dist_V3, dist_H3] = fhindeptest2(x,y);
f_M = ksdensity(dist_M3,xi);
f_W = ksdensity(dist_W3,xi);
f_V = ksdensity(dist_V3,xi);
f_H = ksdensity(dist_H3,xi);
subplot(3,6,3); 
if(plot_pobs)
    scatter(u,v); xlabel('u'); ylabel('v');
else
    scatter(x,y); xlabel('x'); ylabel('y');
end
grid on;
subplot(3,6,9); plot(1:M,filter(b,a,dist_M3),1:M,filter(b,a,dist_W3),1:M,filter(b,a,dist_V3),1:M,filter(b,a,dist_H3));
grid on; title('squared'); ylabel('distance');
subplot(3,6,15); plot(xi,f_M,xi,f_W); grid on;
title(sprintf('%0.02f / %0.02f', ...
               kldivergence(f_M,f_W,xi), ...
               min(abs(skewness(dist_M3)), abs(skewness(dist_W3))) ));
           
x = rand(M,1); y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
u = pobs(x); v = pobs(y);
[dist_M4, dist_W4, dist_V4, dist_H4] = fhindeptest2(x,y);
f_M = ksdensity(dist_M4,xi);
f_W = ksdensity(dist_W4,xi);
f_V = ksdensity(dist_V4,xi);
f_H = ksdensity(dist_H4,xi);
subplot(3,6,4); 
if(plot_pobs)
    scatter(u,v); xlabel('u'); ylabel('v');
else
    scatter(x,y); xlabel('x'); ylabel('y');
end
grid on;
subplot(3,6,10); plot(1:M,filter(b,a,dist_M4),1:M,filter(b,a,dist_W4),1:M,filter(b,a,dist_V4),1:M,filter(b,a,dist_H4));
grid on; title('cubic'); ylabel('distance');
subplot(3,6,16); plot(xi,f_M,xi,f_W); grid on;
title(sprintf('%0.02f / %0.02f', ...
               kldivergence(f_M,f_W,xi), ...
               min(abs(skewness(dist_M4)), abs(skewness(dist_W4))) ));

x = rand(M,1); y = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
u = pobs(x); v = pobs(y);
[dist_M5, dist_W5, dist_V5, dist_H5] = fhindeptest2(x,y);
f_M = ksdensity(dist_M5,xi);
f_W = ksdensity(dist_W5,xi);
f_V = ksdensity(dist_V5,xi);
f_H = ksdensity(dist_H5,xi);
subplot(3,6,5); 
if(plot_pobs)
    scatter(u,v); xlabel('u'); ylabel('v');
else
    scatter(x,y); xlabel('x'); ylabel('y');
end
grid on;
subplot(3,6,11); plot(1:M,filter(b,a,dist_M5),1:M,filter(b,a,dist_W5));
grid on; title('low-freq sin'); ylabel('distance');
subplot(3,6,17); plot(xi,f_M,xi,f_W); grid on;
title(sprintf('%0.02f / %0.02f', ...
               kldivergence(f_M,f_W,xi), ...
               min(abs(skewness(dist_M5)), abs(skewness(dist_W5))) ));

x = rand(M,1); y = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
u = pobs(x); v = pobs(y);
[dist_M6, dist_W6, dist_V6, dist_H6] = fhindeptest2(x,y);
f_M = ksdensity(dist_M6,xi);
f_W = ksdensity(dist_W6,xi);
f_V = ksdensity(dist_V6,xi);
f_H = ksdensity(dist_H6,xi);
subplot(3,6,6); 
if(plot_pobs)
    scatter(u,v); xlabel('u'); ylabel('v');
else
    scatter(x,y); xlabel('x'); ylabel('y');
end
grid on;
subplot(3,6,12); plot(1:M,filter(b,a,dist_M6),1:M,filter(b,a,dist_W6));
grid on; title('fourth-root'); ylabel('distance');
subplot(3,6,18); plot(xi,f_M,xi,f_W); grid on;
title(sprintf('%0.02f / %0.02f', ...
               kldivergence(f_M,f_W,xi), ...
               min(abs(skewness(dist_M6)), abs(skewness(dist_W6))) ));

