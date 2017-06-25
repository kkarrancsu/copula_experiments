function [y2] = reflect_horiz(y1,yval,flag)

if(flag)
    distVec = y1-yval;
    % reflect below
    y2 = yval-distVec;
else
    distVec = yval-y1;
    % reflect above
    y2 = yval+distVec;
end

end