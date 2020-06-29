clear;
close all;

%f = @(x,y,z,r) abs(x).^(0.5) + abs(y).^(0.5) - r;
%f = @(x,y,z) (x.*(x-1).^2.*(x-2) + y.^2).^2 + z.^2-0.01;

sheet = 4;
for k = 2 : 2 : 10
alpha =	0.75;
exponent = 1;
xExponent = 4;
zeta = 0.5;
syms x;
f = @(x,y) x.^xExponent - y.^exponent + exp(x) + 0.5 + alpha.*sin(zeta*k*pi*x);
yFun = @(x) (x.^xExponent + exp(x) + alpha.*sin(zeta*k*pi*x) + 0.5).^(1/exponent);
FBeta = @(beta,px,py,nx,ny) (((px + beta.*nx).^4)+exp(px + beta.*nx) + 0.5 + alpha*sin(zeta*k*pi*(px + beta.*nx))).^(1/exponent) - py - beta.*ny;
ySym = (x.^xExponent + exp(x) + alpha.*sin(zeta*k*pi*x) + 0.5).^(1/exponent);
%FBeta = @(beta,px,py,nx,ny) yFun(px) - py;

% r = 2;
% FBeta = @(beta,px,py,nx,ny) (px + beta*nx).^2 + (py+beta*ny).^2 - r^2;
% f = @(x,y) r.^2 - x.^2 + y.^2;
% yFun = @(x) sqrt( r.^2 - x.^2);

MyDistanceFunction = @(QP, Direction) DistanceFunction([QP(1) QP(2)], FBeta, Direction, yFun);

StartPt = 0;
EndPt = 2;

j = 1;


%Segment = [3 10 100]
%Segment = [100];
Segment = [1000];
%Segment = [40 45 50 55 60]
%Segment = [2 4 6 8 10 11 12 13 14 15]
%Segment = [1]

nSize = length(Segment);

for seg = 1 : nSize
s = int2str(seg);
kstr = int2str(k);
str = strcat('k = ', kstr, ', seg = ',s)
nSegments = Segment(seg);
dx = (EndPt - StartPt)/(nSegments);

XPoly = StartPt:dx:EndPt;
YPoly = yFun(XPoly);
%YPoly = [YPoly(2) YPoly(2)];
XCurve = StartPt:0.001:EndPt;
YCurve = yFun(XCurve);

XPoly = [XPoly EndPt StartPt XPoly(1)];
YPoly = [YPoly 0 0 YPoly(1)];

%XPoly = [XPoly StartPt XPoly(1)];
%YPoly = [YPoly 0 YPoly(1)];

XPoly = fliplr(XPoly)';
YPoly = fliplr(YPoly)';

figure;
hold on;
%plot(XPoly,YPoly,'Color','red','LineWidth',3);
plot(XCurve,YCurve,'Color','blue');
axis equal tight;
hold off;


nVertices = size(XPoly,1);

hold on;

ComputeFlag = ones( nVertices-1, 1);

nStraighEdges = 3;
ComputeFlag(1:nStraighEdges) = zeros(nStraighEdges,1);

Vmax = 0;



for i = 1 : nVertices-1
x1 = XPoly(i);
x2 = XPoly(i+1);
y1 = YPoly(i);
y2 = YPoly(i+1);


%Query Pt
cx = (x2 + x1)/2;
cy = (y2 + y1)/2;
%Edge Direction
e = AW_2D.GetEdgeDirection([x1 y1],[x2,y2]);
ex = e(1);
ey = e(2);
%Normal Direction
n = AW_2D.GetEdgeNormal([x1 y1],[x2,y2]);
nx = n(1);
ny = n(2);

%scatter(x1,y1,'MarkerFaceColor','blue');
%scatter(x2,y2,'MarkerFaceColor','black');


cx1 = cx + 0.3*nx;
cy1 = cy + 0.3*ny;

lx = [cx cx1];
ly = [cy cy1];
[nx ny]';
theta = rad2deg( atan(ny/nx) );
plot(lx,ly,'Color','magenta');
    if( ~ComputeFlag(i) )
        continue;
    end
px = (XPoly(i) + XPoly(i+1))/2;
py = (YPoly(i) + YPoly(i+1))/2;

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

%plot(lx,ly,'Color','magenta');

end
Vmax
px = Pmax(1);
py = Pmax(2);
[d,p] = MyDistanceFunction([px py], [Nx Ny]);
lx = [px p(1)];
ly = [py p(2)];

plot(lx,ly,'Color','black');

kStr = int2str(k);
kStr = strcat({'k = '}, kStr);
title(kStr);


%ComputeFlag = zeros( nVertices-1, 1);

polynomialDegree = 3;
nSamplingPts = 100;
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

Adaptive = true;

%axis equal tight;
%ComputeFlag = zeros( nVertices-1, 1);
A = AW_2D([XPoly YPoly], MyDistanceFunction, 20, Adaptive, yFun );
W = A.GetAdaptiveWeights( IntegrationPoints_XY, polynomialDegree, ComputeFlag );
Area_AW = sum(W)
Area_Analytical = EndPt.^(xExponent+1) / (xExponent+1) + (exp(EndPt) - exp(StartPt)) + (-alpha.*cos(zeta*k*pi*EndPt)/(zeta*k*pi) + alpha.*cos(zeta*k*pi*StartPt)/(zeta*k*pi)) + 0.5*(EndPt - StartPt);
monomials = GetMonomialFunctions();

Area_AW = zeros(10,1);
Area_Analytical = zeros(10,1);

for z = 1 : 10
    f_z = monomials{z};
    for w = 1 : size(W)
        Area_AW(z) = Area_AW(z) + W(w)*f_z(PassPtsX(w),PassPtsY(w));
    end
    Area_Analytical(z) = Analytical_Integral_NonConvex_Symbolic( z, ySym, StartPt , EndPt );
end

Error(j,:) = abs(Area_AW - Area_Analytical);
Velocity_max(j) = Vmax;
j = j + 1;
end

figure;
loglog(Velocity_max,Error);
Table = [Segment' Velocity_max' Error]

kStr = int2str(k);
%FName = strcat('NonCovex_K_',kStr);
%FName = strcat(FName,'.xls');
%FName = 'NonConvex_Final_ArbitraryPolygon.xls';
FName = 'NonConvex_Final_Order_Of_Convergence_Basis_1000.xls';
%xlswrite('Quadrant_New_First_Order.xls', Table);
xlswrite(FName, Table, sheet);
sheet = sheet + 1;
end