function [d,pt] = DistanceFunction_Circle_General(QueryPt,  MinimizerFn, Direction, x_c,y_c,radius )
%UNTITLED11 Summary of this function goes here, y
%   Detailed explanation goes here

px = QueryPt(1);
py = QueryPt(2);
nx = Direction(1);
ny = Direction(2);
FBeta = @(Beta) MinimizerFn(Beta,px,py,nx,ny);

% options = optimset('Display','iter'); 

%y options = optimset(options,'TolX',1e-15);
%initialGuess = QueryPt+ 0.1*[nx ny]; 
%initialGuess = 0;
%initialGuess = (yFun(px)-py );
% initialGuess = 0;
% 
% options = optimoptions(@fsolve,'TolFun',1e-16);
% d = fsolve(FBeta,initialGuess,options);
%d = (yFun(px)-py);
a = nx^2 + ny^2;
b = 2*nx*(px-x_c) + 2*ny*(py-y_c);
c = (px-x_c)^2 + (py-y_c)^2 - radius.^2;

d1 = (-b + sqrt(b^2 - 4*a*c))/(2*a);
d2 = (-b - sqrt(b^2 - 4*a*c))/(2*a);
if( d1 > 0 )
    d = d1;
else
    d = d2;
end

pt = [px; py]+ d*[nx;ny];
%beta_i = 0.1;

end

