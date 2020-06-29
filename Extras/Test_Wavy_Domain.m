clear;
clc;

theta = 0:0.0001:2*pi;
rad = 1;
r = @(t) rad;
k = 2;
alpha = 0.1;
f2 = @(x,y) alpha*exp(x)*sin(k*pi*x);
g1 = @(x) (x-1);
g2 = @(x) (x+1);
g3 = @(y) y;
g4 = @(x,y) ( y - ( (sqrt(1-x.^2)+alpha.*exp(x).*sin(k*pi*x))));
a = 1;
Rfun = @(x1,x2) x1 + x2 - sqrt( x1.^2  + x2.^2 - 2*a*x1.*x2); 
implicitFun = @(x,y) y.*( y - ( (sqrt(1-x.^2)+alpha.*exp(x).*sin(k*pi*x))));
%implicitFun = @(x,y) ( y - ( (sqrt(1-x.^2)+alpha.*exp(x).*sin(k*pi*x))));
%implicitFun = @(x,y) x.^2 + y.^2 - 1
%implicitFun = @(x,y) Rfun(Rfun(Rfun(g1(x),g2(x)),g3(y)),g4(x,y));
axis equal tight;

x = -1:0.01:1;

y = 0:0.01:1;

[X Y] = meshgrid(x,y);
[m,n] = size(X);
XY = zeros(m*n,2);

fval = implicitFun(X,Y);
Z = reshape(fval,m,n);

surf(X,Y,Z);
figure;
contour(X,Y,Z,'ShowText','on');



ezplot(implicitFun,[-1,1,0,2]);

[m,n] = size(theta);

x = zeros(size(theta));
y = zeros(size(theta));
for i = 1 : n/2
x_1 = rad.*cos(theta(i));
y_1 = rad.*sin(theta(i));
x_2 = 0;
y_2 = f2(x_1,y_1);
x(i) = x_1 + x_2;
y(i) = y_1 + y_2;
end


yFun1 = @(x) sqrt(rad.^2 - x.^2) + alpha*exp(x).*sin(k*pi*x);
yFun2 = @(x) - sqrt(rad.^2 - x.^2) + alpha*exp(x).*sin(k*pi*x);
noise1 = @(x) alpha*exp(x).*sin(k*pi*x);
noise2 = @(x) -alpha*exp(x).*sin(k*pi*x);

i1 = integral(yFun1,-1,1)
i2 = integral(yFun2,1,-1)
noise1 = integral(noise1,-1,1)
noise2 = integral(noise2,-1,1)
i = i1 + i2
hold on;
plot(x,y);
hold off;