clear;
clc;
close all;

f = @(X,Y) 1*(((-X.^2 - Y.^2 + ones(size(X)) )) >= zeros(size(X)));
implicitFun = @(X,Y) -X.^2 - Y.^2 + ones(size(X));
fy = @(x) sqrt(1 - x.^2);
bBox.left = [-1;-1];
bBox.right = [1;1];

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
qTree.Partition( 4 );
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