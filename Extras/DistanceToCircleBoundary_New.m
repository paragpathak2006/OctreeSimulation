function [d] = DistanceToCircleBoundary_New( StartPoint, EndPoint, r, T )


m = length(T);
d = zeros(1,m);

for i = 1 : m

t = T(i);
EdgeDirection = AW_2D.GetEdgeDirection(StartPoint,EndPoint);
NormalDirection = AW_2D.GetEdgeNormal(StartPoint,EndPoint);

p = StartPoint' + t.*EdgeDirection';

nx = NormalDirection(1);
ny = NormalDirection(2);
px = p(1);
py = p(2);

a = nx^2 + ny^2;
b = 2*nx*px + 2*ny*py;
c = px^2 + py^2 - r^2;

D = abs(b^2 - 4*a*c);

d1 = (-b + sqrt(D))/(2*a);
d2 = (-b - sqrt(D))/(2*a);

if( d1 >= 0 )
    d(i) = d1;
else
    d(i) = d2;
end
end
end