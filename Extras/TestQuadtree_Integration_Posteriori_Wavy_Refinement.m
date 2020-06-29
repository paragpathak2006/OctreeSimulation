clear;
clc;
close all;


rad = 1;
r = @(t) rad;
k = 10;
alpha = 0.1;
implicitFun = @(x,y) (x-1).*(x+1).*y.*( y - ( (sqrt(1-x.^2)+alpha.*exp(x).*sin(k*pi*x))));
f = @(X,Y) 1*(implicitFun(X,Y)  >= zeros(size(X)));
FBeta = @(beta,px,py,nx,ny) implicitFun(px+beta*nx,py+beta*ny);
fy = [];
distanceToBoundary = @(QP,direction) DistanceFunction(QP, FBeta, direction, fy);

StartX = -1;
EndX = 1;
StartY = 0;
EndY = 1.5;
bBox.left = [StartX;StartY];
bBox.right = [EndX;EndY];


MaxSpaceTreeDepth = 4;

polynomialDegree = 0;
nBasis = 1;
Analytical_Integral = zeros(nBasis,1);
Error_Actual = zeros(MaxSpaceTreeDepth,1);
Integral_AW = zeros(MaxSpaceTreeDepth,nBasis);
Integral_Error = zeros(MaxSpaceTreeDepth,nBasis);
Error_Predicted = zeros( MaxSpaceTreeDepth, 1);

f_integrand = @(x,y) 1; 
Analytical_Integral_f_Integrand = 1.5*pi;

for SpaceTreeDepth = 1 : MaxSpaceTreeDepth
figure;
hold on;
ezplot('(X.^2 + Y.^2 -1).^3 - X.^2.*Y.^3');
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
Area_Actual = (1:MaxSpaceTreeDepth)*Analytical_Integral_f_Integrand;
Error_Actual
Error_Predicted
xlswrite('Error_Analysis_NonConvex_Error_Estimation.xls', [(1:MaxSpaceTreeDepth)' Area_Actual' Integral_f], 1);
