clear;
clc;
close all;

p = 0.5;
q = 0.5;

f1 = @(x,y) (x/7)^2 + (y/3)^2 - 1*ones(size(x));
f2 = @(x,y) abs(x/2)-(3*sqrt(33)-7) * x^2/112 - 3*ones(size(x)) + sqrt(1-(abs(abs(x)-2)-1)^2) - y;
f3 = @(x,y) 9*ones(size(x))- 8*abs(x) - y;
f4 = @(x,y) 3*abs(x) + 0.75*ones(size(x)) - y;
f5 = @(x,y) 2.25*ones(size(x)) + 0*x - y;
f6 = @(x,y) 6*sqrt(10)/7*ones(size(x)) + (1.5-0.5*abs(x)) - 6*sqrt(10)/14*ones(size(x)) * sqrt(4-(abs(x)-1)^2) - y;

h=@(x) 0.5*(
h=12[f(|x+12|+|x?12|+6)?11(x+34)+|x?34|]
 

f = @(x,y) (h?l)H(x+1)+(r?h)H(x?1)+(l?w)H(x+3)+(w?r)H(x?3)+w

implicitFun = @(x,y) f1(x,y)*f2(x,y)*f3(x,y)*f4(x,y)*f5(x,y)*f6(x,y);
fy = @(X) sqrt((1-X.^2)/2) + p*abs(X)^q;
fy1 = @(X) -sqrt((1-X.^2)/2) + p*abs(X)^q;
FBeta = @(beta,px,py,nx,ny) (px + beta*nx).^2 + 2*((py + beta*ny)-p*abs(px + beta*nx)^q).^2  - (px + beta*nx);

distanceToBoundary = @(QP,direction) DistanceFunction(QP, FBeta, direction, fy);
syms x;
ySym =  sqrt((1-x.^2)/2) + p*abs(x)^q;
ySym2 =  -sqrt((1-x.^2)/2) + p*abs(x)^q;


StartX = -6;
EndX = 6;
StartY = -5;
EndY = fy(0);
bBox.left = [StartX;StartY];
bBox.right = [EndX;EndY];

MaxSpaceTreeDepth = 4;

polynomialDegree = 3;
monomials = GetMonomialFunctions();
nBasis = 10;
Analytical_Integral = zeros(nBasis,1);
Error_Actual = zeros(MaxSpaceTreeDepth,1);
Integral_AW = zeros(MaxSpaceTreeDepth,nBasis);
Integral_Error = zeros(MaxSpaceTreeDepth,nBasis);
Error_Predicted = zeros( MaxSpaceTreeDepth, 1);


for z = 1 : 10
    A1 = Analytical_Integral_NonConvex_Symbolic( z, ySym, StartX , EndX);
    A2 = Analytical_Integral_NonConvex_Symbolic( z, ySym2, StartX , EndX);
    Analytical_Integral(z) = A1 + A2;
end

f_integrand = @(x,y) 10*ones(size(x)) + 0.1*x + 0.4*y - 1*x.^2 + 5*x.*y  + 2*y.^2 + 9*x.^3 - 10*x.^2.*y  + 10*x.*y.^2 - 10*y.^3 ; 
Analytical_Integral_f_Integrand = [10 0.1 0.4 -1 5 2 9 -10 10 -10]*Analytical_Integral;

for SpaceTreeDepth = 1 : MaxSpaceTreeDepth
figure;
x = StartX:0.01:EndX;
y1 =  fy(x);
XY = [x' y1'];
y2 =  fy1(x);
XY2 = [x' y2'];
hold on;
plot(XY(:,1),XY(:,2),'b');
plot(XY2(:,1),XY2(:,2),'b');
hold off;



qTree = Quadtree( f,implicitFun, bBox );
qTree.Partition( SpaceTreeDepth );
qTree.DrawTree('red');
[interiorNodes,boundaryNodes] = qTree.GetInterior_BoundaryNodes();

nSamplingPts = 100;

nPt = Cell_2D.GetNumberOfPoints_2D(polynomialDegree);
for i = 1 : length(boundaryNodes)

    node = boundaryNodes(i);
    %node.Draw('green');
    [ApproxPolygon, ComputeSSATerm] = node.GetApproximatePolygon();
    if( ~isempty(ApproxPolygon ) )
        Node_2D.DrawPolygon( ApproxPolygon, 'black');
       m = length( ComputeSSATerm );
       for j = 1 : m
           if( ComputeSSATerm(j) )
               
            if( j+1 < m )
               lineXY = [ApproxPolygon(j,:); ApproxPolygon(j+1,:)];
            else
                 lineXY = [ApproxPolygon(j,:); ApproxPolygon(1,:)];
            end
           line( lineXY(:,1), lineXY(:,2), 'Color','m','LineWidth',2);
           end
       end
    end
end

leafNodes = [interiorNodes; boundaryNodes];
domainDegree = 2;

axis tight equal;

%Error in Basis functions

Integral_f = zeros( MaxSpaceTreeDepth, 1);
for i = 1 : length(leafNodes)
   
   [I,e,XY] = leafNodes(i).GetIntegral( polynomialDegree, distanceToBoundary, domainDegree, f_integrand  );
   Integral_f(SpaceTreeDepth) =  Integral_f(SpaceTreeDepth) + I;
   Error_Predicted(SpaceTreeDepth) = Error_Predicted(SpaceTreeDepth) + e;
   hold on;
      scatter(XY(:,1),XY(:,2));
   hold off;
end


% for i = 1 : length(leafNodes)
%    
%     for z = 1 : 10
%         f_z = monomials{z};
%         Integral_AW(SpaceTreeDepth,z) =  Integral_AW(SpaceTreeDepth,z) + f_z(XY(:,1),XY(:,2))'*W;
%         [I,e,XY] = leafNodes(i).GetIntegral( polynomialDegree, distanceToBoundary, domainDegree, f_z  );
%         Integral_AW(SpaceTreeDepth,z) =  Integral_AW(SpaceTreeDepth,z) + I;
%         Integral_Error(SpaceTreeDepth,z) = Integral_Error(SpaceTreeDepth,z) + e;
%         hold on;
%           scatter(XY(:,1),XY(:,2));
%         hold off;
%     end
% end

Error_Actual(SpaceTreeDepth) = abs(Analytical_Integral_f_Integrand - Integral_f(SpaceTreeDepth))
end
Error_Actual
Error_Predicted
xlswrite('Error_Analysis_NonConvex_Error_Estimation.xls', [(1:MaxSpaceTreeDepth)' Error_Actual abs(Error_Predicted)], 1);
