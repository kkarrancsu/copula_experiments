function d = point_to_line(P, Q1, Q2)
% a = v1 - v2;
% b = pt - v2;
% d = norm(cross(a,b)) / norm(a);

d = ( ( (Q2(:,1)-Q1(:,1)).*(Q1(:,2)-P(:,2)) )- ( (Q1(:,1)-P(:,1)).*(Q2(:,2)-Q1(:,2)) ) ) ./ sqrt( (Q2(:,1)-Q1(:,1)).^2 + (Q2(:,2)-Q1(:,2)).^2 );

end