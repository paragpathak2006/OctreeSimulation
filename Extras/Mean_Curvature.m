clear;
clc;
syms x y;
a = 1;
b = 1;

D = sqrt( (a.^2 - x.^2).^2 + (b.^2 - y.^2)^2);
f = a^2 + b^2 - x^2 - y^2 - D;

fsym = '-(1 + 1 - x.^2 - y.^2 - sqrt( (1 - x.^2)^2 + (1 - y.^2).^2))';
ezplot(fsym,[-3 3 -3 3]);

df_x = diff(f,x,2);
df_y = diff(f,y,2);

divergence = df_x + df_y;

div = matlabFunction(divergence);

% nPoints = 10000;
% 
% ndx = a/nPoints;
% ndy = a/nPoints;
% 
% e1 = -1:ndx:a;*
% e2 = -1:ndy:b;
% e1_xy = [e1' -ones(size(e1'))];
% e2_xy = [ones(size(e1')) e2'];
% e3_xy = [e1' ones(size(e1'))];
% e4_xy = [-ones(size(e1')) e1' ];
% 
% XY = [e1_xy;e2_xy;e3_xy;e4_xy];
% 
% div_val = subs(divergence,'x',XY(:,1));
% div_val = subs(div_val,'y',XY(:,2));
% for i = 1 : size(div_val)
%     try
%     div_val(i) = subs(div_val(i),'y',XY(i,2));
%     catch
%         disp('infinite value at i');
%          XY(i,:)
%         div_val(i) = 0;
%     end
% end

px1 = -1;
py1 = -1;
nx1 = 0;
ny1 = 1;
I1 = @(t) div(px1 + t*nx1,py1 + t*ny1);
px2 = 1;
py2 = -1;
nx2 = 1;
ny2 = 0;
I2 = @(t) div(px2 + t*nx2,py2 + t*ny2);
px3 = 1;
py3 = 1;
nx3 = 0;
ny3 = -1;
I3 = @(t) div(px3 + t*nx3,py3 + t*ny3);
px4 = -1;
py4 = 1;
nx4 = -1;
ny4 = 0;
I4 = @(t) div(px4 + t*nx4,py4 + t*ny4);



dx = 0.001;
x = -0.9:dx:0.9;
y = x;
[X,Y] = meshgrid(x,y);

[m,n] = size(X);

XY = zeros(m*n,2);
k = 1;
for i = 1:m
    for j = 1 :n
        XY(k,1) = X(i,j);
        XY(k,2) = Y(i,j);
        k = k + 1;
    end
end

figure;
z = div(XY(:,1),XY(:,2));

z = reshape(z,m,n);
surf(X,Y,z);


startVal = 1e-9;
endVal = 1 - startVal;
I = @(startVal,endVal) integral(I1,startVal,endVal) + integral(I2,startVal,endVal) + integral(I3,startVal,endVal) + integral(I4,startVal,endVal)

index = 1:size(div_val);

plot(index,div_val);

%+ (y.^2 - ones(size(x))).^2).^(3/2) - (4*y.^2)/((x.^2 - ones(size(x))).^2 + (y.^2 - ones(size(x))).^2).^(1/2) - (2*(x.^2 - ones(size(x))))/((x.^2 - ones(size(x))).^2 + (y.^2 - ones(size(x))).^2).^(1/2) - (2*(y.^2 - ones(size(x))))/((x.^2 - ones(size(x))).^2 + (y.^2 - ones(size(x))).^2).^(1/2) - (4*x.^2)/((x.^2 - ones(size(x))).^2 + (y.^2 - ones(size(x))).^2).^(1/2) + (4*y.^2*(y.^2 - ones(size(x))).^2)/((x.^2 - ones(size(x))).^2 + (y.^2 - ones(size(x))).^2).^(3/2) - 4*ones(size(x));



