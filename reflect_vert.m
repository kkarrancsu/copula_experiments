function [x2] = reflect_vert(x1,xval,flag)

if(flag)
    distVec = xval-x1;
    % reflect to the right
    x2 = distVec+xval;
else
    distVec = x1-xval;
    % reflect to the left
    x2 = xval-distVec;
end

end