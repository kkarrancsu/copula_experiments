%% Test Tau prediction with region shift

clear;
clc;

M = 500;
x = rand(M,1);
noise = 3; num_noise = 30; l = 1;
y1 = x + noise*(l/num_noise)*randn(M,1);
y2 = 4*(x-0.5).^2 + noise*(l/num_noise)*randn(M,1);
y3 = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
y4 = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
y5 = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
y6 = (x > 0.5) + noise*5*l/num_noise*randn(M,1);

% select
y = y2;

[x,y] = pobs_sorted(x,y,1);

numPtsIncr = 25; startIdx = 1; zActualPrev = 0; tauActualPrev = 0;
f = figure;
regionIdxs = [];
for ii=numPtsIncr:numPtsIncr:M-numPtsIncr
    x_cur = x(startIdx:ii); y_cur = y(startIdx:ii);
    x_all = x(1:ii+numPtsIncr); y_all = y(1:ii+numPtsIncr);
    x_actual = x(startIdx:ii+numPtsIncr); y_actual = y(startIdx:ii+numPtsIncr);
    x_buf = x(ii-numPtsIncr+1:ii); x_buf = linspace(min(x_buf),max(x_buf),numPtsIncr)';
    y_buf = y(ii-numPtsIncr+1:ii); 
    
    tauActual = corr(x_actual,y_actual,'type','kendall'); 
    if(tauActual>-1 && tauActual<1)
        rho = copulaparam('Gaussian',tauActual,'type','kendall');
        sigma = 1-abs(rho);
        dither =  randn(numPtsIncr,1)*sigma;
    else
        dither = 0;
    end
    % check the direction
    if((y_buf(1)-y_buf(end))<0)
        y_buf = linspace(min(y_buf),max(y_buf),numPtsIncr)' + dither;
    else
        y_buf = linspace(max(y_buf),min(y_buf),numPtsIncr)' + dither;
    end
    
    y_bufRef1 = reflect_horiz(flipud(y_buf),y_buf(end),1);  % reflect below
    x_bufRef1 = reflect_vert(flipud(x_buf),x_buf(end),1);
    xy_ref1 = [x_bufRef1 y_bufRef1];
    xy_ref2 = [x_bufRef1 flipud(y_buf)];
    
    % find which one we are closer to
    numPtsActual = length(x_actual);
    zActual = tauActual/( (2*(2*numPtsActual))/(9*numPtsActual*(numPtsActual-1)) );
    uv_sameRegion = [[x_cur; xy_ref1(:,1)] [y_cur; xy_ref1(:,2)]]; numPtsRegionTest = size(uv_sameRegion,1);
    uv_diffRegion = [[x_cur; xy_ref2(:,1)]  [y_cur; xy_ref2(:,2)]];
    tauPredSameRegion = corr(uv_sameRegion(:,1),uv_sameRegion(:,2),'type','kendall'); zPredSameRegion = tauPredSameRegion/( (2*(2*numPtsRegionTest))/(9*numPtsRegionTest*(numPtsRegionTest-1)) );
    tauPredDiffRegion = corr(uv_diffRegion(:,1),uv_diffRegion(:,2),'type','kendall'); zPredDiffRegion = tauPredDiffRegion/( (2*(2*numPtsRegionTest))/(9*numPtsRegionTest*(numPtsRegionTest-1)) );
    
    sameRegionZDiff = abs(abs(zActual)-abs(zPredSameRegion));
    diffRegionZDiff = abs(abs(zActual)-abs(zPredDiffRegion));
        
    if((zActualPrev==0 || abs(tauActual)>=abs(tauActualPrev)))
        % means we are continuing in the same direction, don't change
        % anything due to noise variations
    else
        if(diffRegionZDiff<sameRegionZDiff)
            % start a new region
            startIdx = ii;
            x_actual = x(startIdx:ii+numPtsIncr); y_actual = y(startIdx:ii+numPtsIncr);
            regionIdxs = [regionIdxs startIdx];
        end
    end
    
    % plot all the region dileaneations
    for regionIdx=regionIdxs
        plot([x(regionIdx) x(regionIdx)], [0 1], '--','LineWidth',5,'Color','m'); hold on;
    end
    
    dotSize=100;
    scatter(x_all,y_all,dotSize,'b','filled'); hold on;
    dotSize=30;
    scatter(x_actual,y_actual,dotSize,'y','filled'); grid on;
    scatter(uv_sameRegion(:,1),uv_sameRegion(:,2),dotSize,'r+');
    scatter(uv_diffRegion(:,1),uv_diffRegion(:,2),dotSize,'d','g');
    legend({'All','Actual','SameRegion','NewRegion'},'location','SouthWest');
    title({sprintf('\\tau(actual)=%0.02f \\tau(actualPrev)=%0.02f',tauActual,tauActualPrev), ...
           sprintf('sameRegionDiff=%0.02f diffRegionDiff=%0.02f',sameRegionZDiff,diffRegionZDiff)});
    axis([0 1 0 1]);
    
    pause;
    clf(f);
    
    % TODO: move before plots for final algorithm!
    zActualPrev = zActual;
    tauActualPrev = tauActual;
    
end

close all;