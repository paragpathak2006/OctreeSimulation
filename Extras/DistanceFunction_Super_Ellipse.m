function [d,pt] = DistanceFunction_Super_Ellipse(QueryPt,  MinimizerFn, Direction )
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
initialGuess = 0;
%initialGuess = (yFun(px)-py );
%initialGuess = 0;

%options = optimoptions(@fsolve,'TolFun',1e-16);
%d = fzero(FBeta,initialGuess,options);
%d = fzero(FBeta,initialGuess);
%d = (yFun(px)-py);
% a = nx^2 + ny^2;
% b = 2*nx*px + 2*ny*py;
% c = px^2 + py^2 - 1;
% 
% d1 = (-b + sqrt(b^2 - 4*a*c))/(2*a);
% d2 = (-b - sqrt(b^2 - 4*a*c))/(2*a);
% if( d1 > 0 )
%     d = d1;
% else
%     d = d2;
% end

diff = 1e8;
tol = 1e-6;
n = 0.5;
alpha_old = 0;
while( diff > tol )
    alpha_new = ((1- abs( py + alpha_old*ny)^n )^(1/n) - px)/nx;
    diff = abs(alpha_new - alpha_old);
    alpha_old = alpha_new;
end

d = alpha_new;
pt = [px; py]+ d*[nx;ny];
%beta_i = 0.1;

end


