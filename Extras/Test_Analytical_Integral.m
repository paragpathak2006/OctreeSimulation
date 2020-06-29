
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
alpha =	0.75;
exponent = 1;
xExponent = 4;
zeta = 0.5;
syms x;

yFun = @(x) (x.^xExponent + exp(x) + alpha.*sin(zeta*k*pi*x) + 0.5).^(1/exponent);
ySym = (x.^xExponent + exp(x) + alpha.*sin(zeta*k*pi*x) + 0.5).^(1/exponent);
f = {f1, f2, f3, f4, f5, f6, f7, f8, f9, f10};

a = 0;
b = 2;
Integral_Matlab = zeros(10,1);
Integral_Symbolic = zeros(10,1);

for i = 1 : 10
    f_i = f{i};
    Integral_Matlab(i) = Analytical_Integral_NonConvex( f_i, yFun, a , b );
    Integral_Symbolic(i) =  Analytical_Integral_NonConvex_Symbolic(i,ySym,a,b);
end
Integral_Matlab - Integral_Symbolic