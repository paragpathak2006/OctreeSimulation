clear;
clc;
close all;

C = CircularHoleFeature(5,5,10);
C.Polygonize(10);
C.Display('r');
C.DisplayPolygon('b');
C.DisplayNormal('y');
P = C.m_Polygon;
distanceFunction = @(QP,Direction) C.DistanceToBoundary(QP,Direction);
domainDegree = 2;

C2D = Cell_2D(P,distanceFunction,domainDegree,C.m_ComputeSSACorrection);
[Lhs_1,error] = C2D.GetLHS( 3, C.m_ComputeSSACorrection);