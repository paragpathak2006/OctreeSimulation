clear;
clc;
close all;

a = 0;
b = 1;

nDivisions = 100;
dt = (b-a)/nDivisions;
StartPoint = [b,0];
EndPoint = [0,b];
r = 1;

nPoints = nDivisions + 1;
Y = zeros( nPoints, 1);
X = zeros( nPoints, 1);
i = 1;

length = AW_2D.GetEdgeLength(StartPoint,EndPoint);
for t = a : dt : b
    [d,p,pCurve] = DistanceToCircleBoundary(StartPoint,EndPoint,r,t);
    X(i) = t*length;
    Y(i) = d;
    i = i + 1;
end

figure;
plot(X,Y);

nPtRule = [1:12 14 16 21 30 40 50];
%nPtRule = [1:12];
nMaxPt = size(nPtRule,2);
%Integral_Analytical = (pi/4) - 0.5;


a = [1,0];
b = [0,1];



L = AW_2D.GetEdgeLength(a,b);
x = L;
n = 1/2;
%Integral_Analytical = x.^(n+1)/(n+1);
Integral_Analytical = (pi/4) - 0.5;


Integral_Computed = 0;

nElements =2;

Error = zeros( nMaxPt , 1 );

dh = norm((b - a)/ (nElements));

e = AW_2D.GetEdgeDirection(a,b);

for i = 1 : nMaxPt
   Integral_Computed = 0;
for h = 1 : nElements
   StartPoint(1) = a(1) + (h-1)*dh*e(1);
   StartPoint(2) = a(2) + (h-1)*dh*e(2);
   EndPoint = StartPoint + dh*e;
   L = AW_2D.GetEdgeLength(StartPoint,EndPoint);
   nPt = nPtRule(i);
   [ W1 , t1 ] = AW_2D.Gauss_W_X_Scaled( 0,1, nPt );
   d = zeros(nPt,1);
   for j = 1 : nPt
        [d1,~,~] = DistanceToCircleBoundary(StartPoint,EndPoint,r,t1(i));
        %X = t1(j);
        %d1 = X.^n;
        d(j) = d1;
    end
     
Integral_Computed = Integral_Computed + W1'*d*L;
end
Error(i) = 100*abs(Integral_Computed - Integral_Analytical)/Integral_Analytical;
end

figure;
plot(nPtRule,Error);
Error



