function [d,pt] = DistanceFunction_Circle(QueryPt,  MinimizerFn, Direction )
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
d = (yFun(px)-py);
a = nx^2 + ny^2;
b = 2*nx*px + 2*ny*py;
c = px^2 + py^2 - 1;

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

