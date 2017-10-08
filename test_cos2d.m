% Test CoS-2D

clear;
clc;

x = rand(500,1);
y1 = (x-0.5).^2;
y2 = sin(4*pi*x);

y = y2;

subplot(1,3,1);
scatter(x,y);
grid on;

subplot(1,3,2);
scatter(pobs(x),pobs(y))
grid on;

% [cosValue1,RR] = cos2d_v2(x,y1);
[cosValue2,RR] = cos2d_v2(x,y);

% fprintf('CoS(x,y1)=%0.2f CoS(x,y2)=%0.02f\n',cosValue1,cosValue2);

% plot the empirical copula -- color coded by the domain
colors = {'r', 'b'};
numDomains = unique(RR(:,end));
numDomains = numDomains';

subplot(1,3,3);
hold on;
for domain=numDomains
    if(mod(domain,2)==1)
        c = 'r';
    else
        c = 'b';
    end
    II = find(RR(:,end)==domain);
    u = RR(II,2);
    vals = RR(II,end-1);
    plot(u,vals,c);
end
grid on;