%% CIM algorithm speed test for differnet mex configurations

clear;
clc;

numMCSims = 200;
timeVec = zeros(3,numMCSims);  % 1 - normal for-loop
                               % 2 - inline
                               % 3 - inline-all

minScanIncr = 0.015625;
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');
M = 500;
for ii=1:numMCSims
    dispstat(sprintf('%d/%d',ii, numMCSims),'timestamp');
    x = rand(M,1); y = rand(M,1);
    tic;
    cc = cim_v4_cc_mex_plain(x,y,minScanIncr);
    timeVec(1,ii) = toc;
    
    tic;
    cc = cim_v4_cc_mex_inline(x,y,minScanIncr);
    timeVec(2,ii) = toc;
    
    tic;
    cc = cim_v4_cc_mex_numer2(x,y,minScanIncr);
    timeVec(3,ii) = toc;
    
end

p = numSubplots(size(timeVec,1));
subplot(p(1),p(2),1);
histogram(timeVec(1,:)); grid on; title('PLAIN');
subplot(p(1),p(2),2);
histogram(timeVec(2,:)); grid on; title('INLINE');
subplot(p(1),p(2),3);
histogram(timeVec(3,:)); grid on; title('K2');
% subplot(p(1),p(2),4);
% histogram(timeVec(4,:)); grid on; title('INLINE');
% subplot(p(1),p(2),5);
% histogram(timeVec(5,:)); grid on; title('INLINE FOR');

suplabel(sprintf('%d',M),'t');
suplabel('Processing Time (sec)', 'x');

median(timeVec,2)

%% CIM v4/8a/8b mex speed tests

clear;
clc;

numMCSims = 200;
timeVec = zeros(6,numMCSims);  % 1 - v4
                               % 2 - v8a
                               % 3 - v8b
minScanIncr = 0.015625;
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');
M = 500;
for ii=1:numMCSims
    dispstat(sprintf('%d/%d',ii, numMCSims),'timestamp');
    x = rand(M,1); y = rand(M,1);
    tic;
    cc = cim_v4_cc_mex(x,y,minScanIncr);
    timeVec(1,ii) = toc;
    
    tic;
    cc = cim_v8a_cc_mex(x,y,minScanIncr);
    timeVec(2,ii) = toc;
    
    tic;
    cc = cim_v8b_cc_mex(x,y,minScanIncr);
    timeVec(3,ii) = toc;
end

p = numSubplots(size(timeVec,1));
subplot(p(1),p(2),1);
histogram(timeVec(1,:)); grid on; title('v4');
subplot(p(1),p(2),2);
histogram(timeVec(2,:)); grid on; title('v8a');
subplot(p(1),p(2),3);
histogram(timeVec(3,:)); grid on; title('v8b');

suplabel(sprintf('%d',M),'t');
suplabel('Processing Time (sec)', 'x');

median(timeVec,2)