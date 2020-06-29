function [ f ] = GetMonomialFunctions(  )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

f1 = @(x,y) ones(size(x));
f2 = @(x,y) x;
f3 = @(x,y) y;
f4 = @(x,y) x.^2;
f5 = @(x,y) x.*y;
f6 = @(x,y) y.^2;
f7 = @(x,y) x.^3;
f8 = @(x,y) x.^2.*y;
f9 = @(x,y) x.*y.^2;
f10 = @(x,y) y.^3;

f = {f1, f2, f3, f4, f5, f6, f7, f8, f9, f10};


end

