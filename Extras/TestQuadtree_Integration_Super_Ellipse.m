clear;
clc;
close all;

e = 0.5;
f = @(X,Y) 1*(((-abs(X).^e - abs(Y).^e + ones(size(X)) )) >= zeros(size(X)));
implicitFun = @(X,Y) -abs(X).^e - abs(Y).^e + ones(size(X));
FBeta = @(beta,px,py,nx,ny) -abs((px + beta.*nx)).^e - abs((py + beta.*ny)).^e + 1;
distanceToBoundary = @(QP,direction) DistanceFunction_Super_Ellipse(QP, FBeta, direction);
syms x;
ySym = (1 - abs(x).^e).^(1/e);

fy = @(x) (1 - abs(x).^e).^(1/e);
bBox.left = [-1;-1];
bBox.right = [1;1];

MaxSpaceTreeDepth = 1;


polynomialDegree = 3;
monomials = GetMonomialFunctions();
nBasis = 10;
Analytical_Integral = zeros(nBasis,1);
Error = zeros(MaxSpaceTreeDepth,nBasis+1);
Error(:,1) = (1:MaxSpaceTreeDepth)';
Integral_AW = zeros(MaxSpaceTreeDepth,nBasis);


for z = 1 : 10
    A1 = Analytical_Integral_NonConvex_Symbolic( z, ySym, -1 , 1);
    A2 = Analytical_Integral_NonConvex_Symbolic( z, -ySym, 1 , -1);
    Analytical_Integral(z) = A1 + A2;
end

for SpaceTreeDepth = 1 : MaxSpaceTreeDepth
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

for i = 1 : length(leafNodes)
    [XY, W] = leafNodes(i).GetQuadrature_AW( polynomialDegree, distanceToBoundary, domainDegree  );
    hold on;
    scatter(XY(:,1),XY(:,2));
    hold off;
    
    for z = 1 : 10
        f_z = monomials{z};
        Integral_AW(SpaceTreeDepth,z) =  Integral_AW(SpaceTreeDepth,z) + f_z(XY(:,1),XY(:,2))'*W;
    end
end

Error(SpaceTreeDepth,2:end) = abs(Analytical_Integral' - Integral_AW(SpaceTreeDepth,:));
end
Error
