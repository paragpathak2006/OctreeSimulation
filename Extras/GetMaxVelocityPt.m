function [P,V] = GetMaxVelocityPt( x1, y1, x2, y2, MyDistanceFunction, n )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


fun = @(t) -(MyDistanceFunction( [x1 + t.*(x2-x1) y1+t.*(y2-y1)], n )).^2;

t0 = 0.5;

t = fmincon(fun,t0,[],[],[],[],0,1);

P = [x1 + t*(x2-x1) y1+t*(y2-y1)];
V = abs(MyDistanceFunction( [x1 + t*(x2-x1) y1+t*(y2-y1)], n ));

end

