function [m,b4] = numUnique(a)

b4 = find(accumarray(a(:,1)+1,1)-1);
b4(b4==0) = [];
m = length(b4);
end