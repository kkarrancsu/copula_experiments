function [ dist ] = empCopulaDistance(C1,C2,metric)
if(strcmpi(metric,'sse'))
    errC = C1-C2;
    errC = errC(:);
    dist = sum(errC.^2)/length(errC);
end
end