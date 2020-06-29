function [ val ] = Analytical_Integral_NonConvex( integrand, yFun, a , b )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

c = 0;
d = @(x) yFun(x);

val = quad2d(integrand,a,b,c,d);

end

