%% Test POBS Generation
clear;
clc;

M = 500;
x = rand(M,1); y = 4*(x-0.5).^2; [u1,v1] = pobs_sorted(x,y,1); UV1 = [u1 v1];

% second way of generating pobs not using tiedrank method b/c it is not
% supported by matlab coder
data_sorted = sort(x);
[~, u2] = ismember(x,data_sorted);
data_sorted = sort(y);
[~, v2] = ismember(y,data_sorted);
u2 = u2/M; v2 = v2/M; UV2 = [u2 v2];

subplot(1,2,1); scatter(u1,v1); subplot(1,2,2); scatter(u2,v2);