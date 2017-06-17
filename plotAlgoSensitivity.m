function [] = plotAlgoSensitivity(sensitivityData, scanincrs, noiseVec)

plotStyle = {'o-.', '+-.', 'd-.', 'v-.', 's-.', 'p-.'};

cimVecIdx = 3;

subplot(3,3,1);
legendCell = cell(1,length(scanincrs));
hhCell = cell(1,length(scanincrs));
for ii=1:length(scanincrs)
    rawData = sensitivityData{ii};
    hh1 = plot(noiseVec,rawData.linearDep(cimVecIdx,:),plotStyle{ii});
    hhCell{ii} = hh1;
    legendCell{ii} = num2str(scanincrs(ii));
    hold on;
end
grid on;
xlabel('noise');
legend(legendCell);
title('Linear Dependency');
for ii=1:length(hhCell)
    hhCell{ii}.LineWidth = 1.5;
end

subplot(3,3,2);
hhCell = cell(1,length(scanincrs));
for ii=1:length(scanincrs)
    rawData = sensitivityData{ii};
    hh1 = plot(noiseVec,rawData.quadraticDep(cimVecIdx,:),plotStyle{ii});
    hhCell{ii} = hh1;
    hold on;
end
grid on;
xlabel('noise');
title('Quadratic Dependency');
for ii=1:length(hhCell)
    hhCell{ii}.LineWidth = 1.5;
end

subplot(3,3,3);
hhCell = cell(1,length(scanincrs));
for ii=1:length(scanincrs)
    rawData = sensitivityData{ii};
    hh1 = plot(noiseVec,rawData.cubicDep(cimVecIdx,:),plotStyle{ii});
    hhCell{ii} = hh1;
    hold on;
end
grid on;
xlabel('noise');
title('Cubic Dependency');
for ii=1:length(hhCell)
    hhCell{ii}.LineWidth = 1.5;
end

subplot(3,3,4);
hhCell = cell(1,length(scanincrs));
for ii=1:length(scanincrs)
    rawData = sensitivityData{ii};
    hh1 = plot(noiseVec,rawData.sinusoidalDep(cimVecIdx,:),plotStyle{ii});
    hhCell{ii} = hh1;
    hold on;
end
grid on;
xlabel('noise');
title('Sinusoidal Dependency');
for ii=1:length(hhCell)
    hhCell{ii}.LineWidth = 1.5;
end

subplot(3,3,5);
hhCell = cell(1,length(scanincrs));
for ii=1:length(scanincrs)
    rawData = sensitivityData{ii};
    hh1 = plot(noiseVec,rawData.hiFreqSinDep(cimVecIdx,:),plotStyle{ii});
    hhCell{ii} = hh1;
    hold on;
end
grid on;
xlabel('noise');
title('HF-Sin Dependency');
for ii=1:length(hhCell)
    hhCell{ii}.LineWidth = 1.5;
end

subplot(3,3,6);
hhCell = cell(1,length(scanincrs));
for ii=1:length(scanincrs)
    rawData = sensitivityData{ii};
    hh1 = plot(noiseVec,rawData.fourthRootDep(cimVecIdx,:),plotStyle{ii});
    hhCell{ii} = hh1;
    hold on;
end
grid on;
xlabel('noise');
title('Fourth-Root Dependency');
for ii=1:length(hhCell)
    hhCell{ii}.LineWidth = 1.5;
end

subplot(3,3,7);
hhCell = cell(1,length(scanincrs));
for ii=1:length(scanincrs)
    rawData = sensitivityData{ii};
    hh1 = plot(noiseVec,rawData.circleDep(cimVecIdx,:),plotStyle{ii});
    hhCell{ii} = hh1;
    hold on;
end
grid on;
xlabel('noise');
title('Circular Dependency');
for ii=1:length(hhCell)
    hhCell{ii}.LineWidth = 1.5;
end

subplot(3,3,8);
hhCell = cell(1,length(scanincrs));
for ii=1:length(scanincrs)
    rawData = sensitivityData{ii};
    hh1 = plot(noiseVec,rawData.stepDep(cimVecIdx,:),plotStyle{ii});
    hhCell{ii} = hh1;
    hold on;
end
grid on;
xlabel('noise');
title('Step-Function Dependency');
for ii=1:length(hhCell)
    hhCell{ii}.LineWidth = 1.5;
end

subplot(3,3,9);
hhCell = cell(1,length(scanincrs));
for ii=1:length(scanincrs)
    rawData = sensitivityData{ii};
    hh1 = plot(noiseVec,rawData.indep(cimVecIdx,:),plotStyle{ii});
    hhCell{ii} = hh1;
    hold on;
end
grid on;
xlabel('noise');
title('Independence');
for ii=1:length(hhCell)
    hhCell{ii}.LineWidth = 1.5;
end

end