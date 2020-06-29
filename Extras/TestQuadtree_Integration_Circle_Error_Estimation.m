clear;
clc;
close all;

f = @(X,Y) 1*(((-X.^2 - Y.^2 + ones(size(X)) )) >= zeros(size(X)));
implicitFun = @(X,Y) -X.^2 - Y.^2 + ones(size(X));
FBeta = @(beta,px,py,nx,ny) (px + beta.*nx).^2 + (py + beta.*ny).^2 - 1;
distanceToBoundary = @(QP,direction) DistanceFunction_Circle(QP, FBeta, direction);
syms x;
ySym = sqrt(1 - x.^2);

fy = @(x) sqrt(1 - x.^2);
bBox.left = [-1;-1];
bBox.right = [1;1];

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
    A1 = Analytical_Integral_NonConvex_Symbolic( z, ySym, -1 , 1);
    A2 = Analytical_Integral_NonConvex_Symbolic( z, -ySym, 1 , -1);
    Analytical_Integral(z) = A1 + A2;
end

f_integrand = @(x,y) 10*ones(size(x)) + 0.1*x + 0.4*y - 1*x.^2 + 5*x.*y  + 2*y.^2 + 9*x.^3 - 10*x.^2.*y  + 10*x.*y.^2 - 10*y.^3 ; 
Analytical_Integral_f_Integrand = [10 0.1 0.4 -1 5 2 9 -10 10 -10]*Analytical_Integral;

for SpaceTreeDepth = 2 : MaxSpaceTreeDepth
figure;
x = -1:0.01:1;
y1 =  fy(x);
y2 = -y1;
XY = [x' y1'];
hold on;
plot(XY(:,1),XY(:,2),'b');
XY = [x' y2'];
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

Error_Actual(SpaceTreeDepth) = abs(Analytical_Integral_f_Integrand - Integral_f(SpaceTreeDepth));
end
Error_Actual
Error_Predicted
xlswrite('Error_Analysis_Circle_Cubic_function.xls', [(1:MaxSpaceTreeDepth)' Error_Actual abs(Error_Predicted)], 1);
