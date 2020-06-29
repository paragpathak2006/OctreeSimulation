clear;
clc;
close all;


StartPoint = [1 0];
EndPoint = [0 1];
r = 1;
d = @(t) DistanceToCircleBoundary_New(StartPoint,EndPoint,r,t);
L = AW_2D.GetEdgeLength(StartPoint,EndPoint);

Integral_C = integral(d,0,L,'RelTol',0,'AbsTol',1e-12)

Analytical_Integral = pi/4 - 0.5

Error = 100*abs( (Analytical_Integral - Integral_C)/Analytical_Integral)

