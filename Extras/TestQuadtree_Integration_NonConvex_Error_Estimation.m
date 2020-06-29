clear;
clc;
close all;


k = 3;
alpha =	0.75;
exponent = 1;
xExponent = 4;
zeta = 0.5;
f = @(X,Y) 1*(X.^xExponent - Y.^exponent + exp(X) + 0.5*ones(size(X)) + alpha.*sin(zeta*k*pi*X) >= zeros(size(X)));
implicitFun = @(X,Y) X.^xExponent - Y.^exponent + exp(X) + 0.5*ones(size(X)) + alpha.*sin(zeta*k*pi*X);
fy = @(x) (x.^xExponent + exp(x) + alpha.*sin(zeta*k*pi*x) + 0.5).^(1/exponent);
FBeta = @(beta,px,py,nx,ny) (((px + beta.*nx).^4)+exp(px + beta.*nx) + 0.5 + alpha*sin(zeta*k*pi*(px + beta.*nx))).^(1/exponent) - py - beta.*ny;


distanceToBoundary = @(QP,direction) DistanceFunction(QP, FBeta, direction, fy);
syms x;
ySym =  (x.^xExponent + exp(x) + alpha.*sin(zeta*k*pi*x) + 0.5).^(1/exponent);

StartX = 0;
EndX = 1;
StartY = 0;
EndY = fy(EndX);
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
    Analytical_Integral(z) = A1;
end

f_integrand = @(x,y) 10*ones(size(x)) + 0.1*x + 0.4*y - 1*x.^2 + 5*x.*y  + 2*y.^2 + 9*x.^3 - 10*x.^2.*y  + 10*x.*y.^2 - 10*y.^3 ; 
Analytical_Integral_f_Integrand = [10 0.1 0.4 -1 5 2 9 -10 10 -10]*Analytical_Integral;

for SpaceTreeDepth = 2 : MaxSpaceTreeDepth
figure;
x = StartX:0.01:EndX;
y1 =  fy(x);
XY = [x' y1'];
hold on;
plot(XY(:,1),XY(:,2),'b');
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
   if( abs(e) < eps )
      scatter(XY(:,1),XY(:,2),'b','filled');
   else
       scatter(XY(:,1),XY(:,2),'+','r');
   end
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
