clear;
clc;
alpha =	0.75;
k = 2;
exponent = 1;
yFun = @(x) (x.^4 + exp(x) + alpha.*sin(k*pi*x) + 0.5).^(1/exponent);

StartPt = 0;
EndPt = 2;

nSegments =3;
dx = (EndPt - StartPt)/(nSegments);

XPoly = StartPt:dx:EndPt;
YPoly = yFun(XPoly);
XPoly = [XPoly EndPt StartPt XPoly(1)];
YPoly = [YPoly 0 0 YPoly(1)];

XPoly = fliplr(XPoly)';
YPoly = fliplr(YPoly)';

XY = [0 1.500000000000000;
      0 0;
      2 0;
      2 23.889056098930652;
      1.333333333333334   8.103680774682003;
      2/3   1.995745852413878;
      0   1.5;];
 
  X = XPoly;
  Y = YPoly;
  
  b1 = X(6) - X(1);
  b2 = X(5) - X(6);
  b3 = X(4) - X(5);

  h1 = Y(1) - Y(2);        
  h2 = Y(6) - Y(2);
  h3 = Y(5) - Y(2);
  h4 = Y(4) - Y(2);
  
  A1 = 0.5*b1*(h1 + h2 )
  A2 = 0.5*b2*(h2 + h3 )
  A3 = 0.5*b3*(h3 + h4 )
  
  T1 = integral(yFun,X(1),X(6));
  T2 = integral(yFun,X(6),X(5));
  T3 = integral(yFun,X(5),X(4));
  
  B1 = T1 - A1
  B2 = T2 - A2
  B3 = T3 - A3
  
  A1 + A2 + A3