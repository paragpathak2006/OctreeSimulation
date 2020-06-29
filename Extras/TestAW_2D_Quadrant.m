clear;
close all;

%f = @(x,y,z,r) abs(x).^(0.5) + abs(y).^(0.5) - r;
%f = @(x,y,z) (x.*(x-1).^2.*(x-2) + y.^2).^2 + z.^2-0.01;
%alpha =	0.75;
%k = 7;
%exponent = 1;
%f = @(x,y) x.^4+exp(x) - y.^exponent + 0.5 + alpha.*sin(k*pi*x);
%yFun = @(x) (x.^4 + exp(x) +  alpha.*sin(k*pi*x) + 0.5).^(1/exponent);
%FBeta = @(beta,px,py,nx,ny) (((px + beta.*nx).^4)+exp(px + beta.*nx) + 0.5 + alpha*sin(k*pi*(px + beta.*nx))).^(1/exponent) - py - beta.*ny;

r = 1;
exponent = 2;
FBeta = @(beta,px,py,nx,ny) (px + beta*nx).^exponent + (py+beta*ny).^exponent - r^exponent;
f = @(x,y) r.^exponent - x.^exponent + y.^exponent;
yFun = @(x) real(( r.^exponent - x.^exponent).^(1/exponent));

MyDistanceFunction = @(QP, Direction) DistanceFunction([QP(1) QP(2)], FBeta, Direction, yFun);

StartPt = 0;
EndPt = 1;

j = 1;


Segment = [1:2:20 100 1000]
%Segment = 1;
nSize = length(Segment);
%nSize = 1;
for z = 1 : nSize

nSegments = Segment(z);
dx = (EndPt - StartPt)/(nSegments);

XPoly = StartPt:dx:EndPt;
YPoly = yFun(XPoly);
XCurve = StartPt:0.001:EndPt;
YCurve = yFun(XCurve);

%XPoly = [XPoly EndPt StartPt XPoly(1)];
%YPoly = [YPoly 0 0 YPoly(1)];

XPoly = [XPoly StartPt XPoly(1)];
YPoly = [YPoly 0 YPoly(1)];

XPoly = fliplr(XPoly)';
YPoly = fliplr(YPoly)';

figure;
hold on;
plot(XPoly,YPoly,'Color','red','LineWidth',3);
plot(XCurve,YCurve,'Color','blue');
hold off;


nVertices = size(XPoly,1);

hold on;

ComputeFlag = ones( nVertices-1, 1);

ComputeFlag(1:2) = zeros(2,1);

Vmax = 0;



for i = 1 : nVertices-1
    if( ~ComputeFlag(i) )
        continue;
    end
%px = (XPoly(i) + XPoly(i+1))/2;
%py = (YPoly(i) + YPoly(i+1))/2;



x1 = XPoly(i);
x2 = XPoly(i+1);
y1 = YPoly(i);
y2 = YPoly(i+1);


ex = (x2 - x1);
ey = (y2 - y1);
e = [ex ey]';
e = e / norm(e);
ex = e(1);
ey = e(2);
nx = ey;
ny = -ex;

[P V] = GetMaxVelocityPt( XPoly(i), YPoly(i), XPoly(i+1), YPoly(i+1), MyDistanceFunction, [nx ny] );
px = P(1);
py = P(2);

if( Vmax < V )
    Vmax = V;
    Pmax = P;
    Nx = nx;
    Ny = ny;
end

[d,p] = MyDistanceFunction([px py], [nx ny]);

%scatter(px,py);
lx = [px p(1)];
ly = [py p(2)];

plot(lx,ly,'Color','magenta');

end
Vmax
px = Pmax(1);
py = Pmax(2);
[d,p] = MyDistanceFunction([px py], [Nx Ny]);
lx = [px p(1)];
ly = [py p(2)];

plot(lx,ly,'Color','black');

polynomialDegree = 0;
nSamplingPts = 15;
nPt = AW_2D.GetNumberOfPoints_2D( polynomialDegree );
hx = (EndPt - StartPt)/(nSamplingPts+1);
XTestPts = StartPt+hx:hx:EndPt-hx;
minY = 0;
maxY = max(YCurve);
hy = (maxY - minY)/(nSamplingPts+1);
YTestPts = minY + hy : hy : maxY - hy;

[XTestPts, YTestPts] = meshgrid(XTestPts,YTestPts); 



fVal = f(XTestPts,YTestPts);
[m, n] = size( XTestPts);
X= reshape(XTestPts,[m*n,1]);
Y = reshape(YTestPts,[m*n,1]);
fVal = reshape(fVal,[m*n,1]);

PassPtsX = X(fVal >0);
PassPtsY = Y(fVal >0);

[Index] = inpolygon( PassPtsX, PassPtsY, XPoly, YPoly);

PassPtsX = PassPtsX(Index);
PassPtsY = PassPtsY(Index);

r = randperm(size(PassPtsX,1));

PassPtsX = PassPtsX(r(1:nPt));
PassPtsY = PassPtsY(r(1:nPt));

scatter(PassPtsX, PassPtsY);
hold off;

IntegrationPoints_XY = [PassPtsX PassPtsY];

AdaptiveQuadrature = true;
ComputeFlag = zeros( nVertices-1, 1);
domainDegree = 13;
A = AW_2D([XPoly YPoly], MyDistanceFunction,domainDegree, AdaptiveQuadrature, yFun );
W = A.GetAdaptiveWeights( IntegrationPoints_XY, polynomialDegree, ComputeFlag );
Area_AW = sum(W)
Area_Analytical = pi/4;
Error(j) = abs(Area_AW - Area_Analytical);
Velocity_max(j) = Vmax;
nBoundaryPts(j) = A.GetBoundaryIntegrationOrder(polynomialDegree);
j = j + 1;
end

loglog(Velocity_max,Error);
Table = [Segment' nBoundaryPts' Velocity_max' Error']
%xlswrite('Quadrant_New_First_Order_New.xls', Table);
%xlswrite('Quadrant_New_First_Order_New_Ordinary.xls', Table);

fileName = 'Quadrant_Area_SSA_0.xls';
xlswrite(fileName, Table,7);
%xlswrite('Quadrant_New_Zero_Order.xls', Table);




