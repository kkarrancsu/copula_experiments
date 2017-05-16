function [ dist ] = empCopulaDistance(C1,C2,metric)
if(strcmpi(metric,'sse'))
    errC = C1-C2;
    dist = sum(sum(errC.^2));
end
end