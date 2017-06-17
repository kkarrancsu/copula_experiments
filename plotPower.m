function [] = plotPower(powerMat, M, labels, noiseVec, noiseMin, noiseMax, style)

if(style==1)
    numDepTypes = 8;
    % the original Simon/Tibshirani style of power-plots
    plotStyle = {'o-.', '+-.', 'd-.', 'v-.', 's-.', 'p-.'};
    % inlet plot configuration
    M_inlet = 200;
    if(M<=500)
        inset_bufX = 0.0005; inset_bufY = 0.002;
    else
        inset_bufX = 0.15; inset_bufY = 0.26;
    end
    inset_width = 0.1; inset_height = 0.08;
    
    % generate the inlet data
    inletX = linspace(0,1,M_inlet);
    inletData = zeros(numDepTypes,M_inlet);
    inletData(1,:) = inletX;
    inletData(2,:) = 4*(inletX-.5).^2;
    inletData(3,:) = 128*(inletX-1/3).^3-48*(inletX-1/3).^3-12*(inletX-1/3);
    inletData(4,:) = sin(4*pi*inletX);
    inletData(5,:) = sin(16*pi*inletX);
    inletData(6,:) = inletX.^(1/4);
    inletData(7,:) = (sqrt(1 - (2*inletX - 1).^2));
    inletData(8,:) = (inletX > 0.5);

    % we break it up into 2 figures, 2x2 subplots on each figure for
    % readability
    for figNum=[1,2]
        figure(figNum);
        for subplotNum=1:4
            depTypeIdx = (figNum-1)*4+subplotNum;
            hhCell = cell(1,length(labels));
            for ii=1:length(labels)
                h = subplot(2,2,subplotNum);
                hh = plot(noiseVec, squeeze(powerMat(ii,depTypeIdx,noiseMin:noiseMax)), plotStyle{ii});
                hhCell{ii} = hh;
                hold on;
            end
            axis([min(noiseVec) max(noiseVec) 0 1]);
            xlabel({'Noise Level','(a)'}, 'FontSize', 20); ylabel('Power', 'FontSize', 20); grid on;
            h.FontSize = 20; 
            loc_inset = [h.Position(1)+inset_bufX h.Position(2)+inset_bufY inset_width inset_height];
            ax = axes('Position',loc_inset);
            tmp1 = linspace(0,1,M_inlet);
            tmp2 = tmp1;
            plot(inletX,inletData(depTypeIdx,:), 'k', 'LineWidth', 2);
            ax.Box = 'on'; ax.XTick = []; ax.YTick = [];
            ax.XLim = [min(inletX) max(inletX)];
            ax.YLim = [min(inletData(depTypeIdx,:)) max(inletData(depTypeIdx,:))];
            for ii=1:length(labels)
                hhCell{ii}.LineWidth = 1.5; 
            end
        end
        if(figNum==2)
            legend(h,labels);
        end
    end
    
else
    error('Current style not supported!!');
end

end