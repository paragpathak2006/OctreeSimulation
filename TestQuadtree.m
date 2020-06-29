clear;
clc;
close all;
format compact;

load MarchingCubes.txt; global Cube;    Cube = MarchingCubes;
load Polyhedras.txt;    global FV;      FV = Polyhedras;

r=.2;
R=.8;
%   f = @(X,Y,Z) 1*(((-X.^2 - Y.^2 - Z.^2 + ones(size(X)))) >= zeros(size(X)));
%   implicitFun = @(X,Y,Z) -X.^2 - Y.^2  - Z.^2 + ones(size(X));


  f = @(X,Y,Z) r^2*ones(size(X))-(R*ones(size(X))-sqrt(X.^2+Y.^2)).^2-Z.^2 >= zeros(size(X));
  implicitFun = @(X,Y,Z) r^2*ones(size(X))-(R*ones(size(X))-sqrt(X.^2+Y.^2)).^2-Z.^2;

fz = @(x,y) sqrt(1 - x.^2 - y.^2);
bBox.left = [-1;-1;-1];
bBox.right = [1;1;1];

figure;
x = -1:0.01:1;
y = -1:0.01:1;
z1 =  fz(x,y);
z2 = -z1;
XYZ = [x' y' z1'];
hold on;
% plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'b');
XYZ = [x' y' z2'];
% plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'b');
hold off;
qTree = Quadtree( f,implicitFun, bBox );
qTree.Partition(3);
qTree.DrawTree('red');
[interiorNodes,boundaryNodes] = qTree.GetInterior_BoundaryNodes();

% nSamplingPts = 100;
% polynomialDegree = 3;
% nPt = Cell_2D.GetNumberOfPoints_2D(polynomialDegree);
% for i = 1 : length(boundaryNodes)
% 
%     node = boundaryNodes(i);
%     %node.Draw('green');
%     [ApproxPolygon, ComputeSSATerm] = node.GetApproximatePolygon();
%     if( ~isempty(ApproxPolygon ) )
%         Node_2D.DrawPolygon( ApproxPolygon, 'black');
%        m = length( ComputeSSATerm );
%        for j = 1 : m
%            if( ComputeSSATerm(j) )
%                
%             if( j+1 < m )
%                lineXY = [ApproxPolygon(j,:); ApproxPolygon(j+1,:)];
%             else
%                  lineXY = [ApproxPolygon(j,:); ApproxPolygon(1,:)];
%             end
%            line( lineXY(:,1), lineXY(:,2), 'Color','m','LineWidth',2);
%            end
%        end
%        Pts = node.GetRandomQuadraturePts( nPt, ApproxPolygon, nSamplingPts );
%        hold on;
%        scatter(Pts(:,1),Pts(:,2));
%        hold off;
%        
%     end
% end

% for i = 1 : length(interiorNodes)
%     node = interiorNodes(i);
%     CornerPts = node.GetCornerPts();
%     [XY,W] = Cell_2D.GetCartesianQuadrature_Box( CornerPts, 3 );
%     hold on;
%     scatter(XY(:,1),XY(:,2));
%     hold off;
%     node.Draw('black');
% end