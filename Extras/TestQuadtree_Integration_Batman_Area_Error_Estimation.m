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

syms x y;
eq1 = '(x/7)^2 + (y/3)^2 - 1';
eq2 = 'abs(x/2)-(3*sqrt(33)-7) * x^2/112 - 3 + sqrt(1-(abs(abs(x)-2)-1)^2) - y';
eq3 = '9 -8*abs(x) - y';
eq4 = '3*abs(x) + 0.75 - y';
eq5 = '2.25+ 0*x - y';
eq6 = '6*sqrt(10)/7 + (1.5-0.5*abs(x)) - 6*sqrt(10)/14  * sqrt(4-(abs(x)-1)^2) - y';

implicitFun = @(x,y) f1(x,y)*f2(x,y)*f3(x,y)*f4(x,y)*f5(x,y)*f6(x,y);
f = @(x,y) 1*(implicitFun(x,y)>=0);
FBeta = @(beta,px,py,nx,ny) implicitFun(px+beta*nx,py+beta*ny);

distanceToBoundary = @(QP,direction) DistanceFunction(QP, FBeta, direction, fy);


StartX = -7.25;
EndX = 7.25;
StartY = -5;
EndY = 5;
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
Analytical_Integral_f_Integrand = 955/48 -(2/7)*(2*sqrt(33) + 7*pi + 3*sqrt(10)*(pi-1) ) + 21*(acos(3/7) + acos(4/7));

for SpaceTreeDepth = 4 : MaxSpaceTreeDepth
figure;
x = StartX:0.01:EndX;
axes('Xlim', [-7.25 7.25], 'Ylim', [-5 5]);
hold on
ezplot(eq1,[-8 8 -3*sqrt(33)/7 6-4*sqrt(33)/7]);
ezplot(eq2,[-4 4]);
ezplot(eq3,[-1 -0.75 -5 5]);
ezplot(eq3,[0.75 1 -5 5]);
ezplot(eq4,[-0.75 0.75 2.25 5]);
ezplot(eq5,[-0.5 0.5 -5 5]);
ezplot(eq6,[-3 -1 -5 5]);
ezplot(eq6,[1 3 -5 5]);
colormap([0 0 1])
title('Batman');
xlabel('');
ylabel('');
hold off



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
Actual_Area = ones(1:MaxSpaceTreeDepth)*Analytical_Integral_f_Integrand;
Error_Actual
Error_Predicted
xlswrite('Error_Analysis_Batman_Curve_Error_Estimation.xls', [(1:MaxSpaceTreeDepth)' Actual_Area' Integral_f], 1);
