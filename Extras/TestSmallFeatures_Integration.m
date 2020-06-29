clear;
clc;
close all;


rad = 1;
r = @(t) rad;
k = 2;
alpha = 0.1;
f2 = @(x,y) alpha*exp(x)*sin(k*pi*x);
implicitFun = @(x,y) ( y - ( (sqrt(1-x.^2)+alpha.*exp(x).*sin(k*pi*x))));

%implicitFun2 = @(x,y) (x+1).*(x-0.99).*y.*( y - ( (sqrt(1-x.^2)+alpha.*exp(x).*sin(k*pi*x))));

f = @(X,Y) 1*(-implicitFun(X,Y) >= zeros(size(X)));

fy = @(x) ( (sqrt(1-x.^2)+alpha.*exp(x).*sin(k*pi*x)));
FBeta = @(beta,px,py,nx,ny) implicitFun(px + beta*nx, py + beta*ny);
distanceToBoundary = @(QP,direction) DistanceFunction(QP, FBeta, direction, fy);
syms x;
ySym =  ( (sqrt(1-x.^2)+alpha.*exp(x).*sin(k*pi*x)));

StartX = -0.8;
EndX = 0.99;
StartY = 0;
EndY = 1.3;
bBox.left = [StartX;StartY];
bBox.right = [EndX;EndY];

SmallFeatureRadius = [0.025; 0.025/2; 0.025/3; 0.025/4; 0.025/5; 0.025/6];
nEdges = [4 6 8 10 15 20 50 100 1000];
nR = length( SmallFeatureRadius );
nE = length(nEdges);

polynomialDegree = 3;
monomials = GetMonomialFunctions();
nBasis = 10;
Analytical_Integral = zeros(nBasis,1);

for z = 1 : 10
    A1 = Analytical_Integral_NonConvex_Symbolic( z, ySym, StartX , EndX);
    Analytical_Integral(z) = A1;
end

f_integrand = @(x,y) 10*ones(size(x)) + 0.1*x + 0.4*y - 1*x.^2 + 5*x.*y  + 2*y.^2 + 9*x.^3 - 10*x.^2.*y  + 10*x.*y.^2 - 10*y.^3 ; 
Analytical_Integral_f_Integrand = [10 0.1 0.4 -1 5 2 9 -10 10 -10]*Analytical_Integral

SpaceTreeDepth = 2;
qTree = Quadtree( f,implicitFun, bBox );
qTree.Partition( SpaceTreeDepth );
[interiorNodes,boundaryNodes] = qTree.GetInterior_BoundaryNodes();
nSmallFeatures = 100;
maxRadius = 0.025;

for nNodes = 1 : length(boundaryNodes)
    
    bbox = boundaryNodes(nNodes).m_bbox;
    dx = abs(bbox(2,1) - bbox(1,1));
    dy = abs(bbox(2,2) - bbox(1,2));
    
    x_s = bbox(1,1);
    y_s = bbox(1,2);
    
    dh = min(dx,dy);
    X_S = x_s + maxRadius : 2.5*maxRadius:x_s + dh - maxRadius;
    Y_S = y_s + maxRadius : 2.5*maxRadius:y_s + dh - maxRadius;
    [X,Y] = meshgrid(X_S,Y_S);
    fVal = boundaryNodes(nNodes).m_PMC(X,Y);
    fVal1 =  boundaryNodes(nNodes).m_PMC(X+maxRadius,Y+maxRadius);
    fVal2 =  boundaryNodes(nNodes).m_PMC(X-maxRadius,Y-maxRadius);
    fVal3 =  boundaryNodes(nNodes).m_PMC(X-maxRadius,Y+maxRadius);
    fVal4 =  boundaryNodes(nNodes).m_PMC(X+maxRadius,Y-maxRadius);
    fVal = and( fVal, fVal1);
    fVal = and( fVal, fVal2);
    fVal = and( fVal, fVal3);
    fVal = and( fVal, fVal4);
    X = X(fVal > 0);
    Y = Y(fVal > 0);
    
    maxFeaturesPossible = min(nSmallFeatures,length(X));
    thisRadius = ones(maxFeaturesPossible,1).*maxRadius;
    r = randperm(size(X,1));
    X = X(r);
    Y = Y(r);
    
    SmallFeatureList = cell(maxFeaturesPossible,1);
    for nFeatures = 1 : maxFeaturesPossible
        SmallFeatureList{nFeatures} = CircularHoleFeature(X(nFeatures),Y(nFeatures),thisRadius(nFeatures));
        S = SmallFeatureList{nFeatures};
        A_S_I = real(S.GetAnalyticalIntegral(polynomialDegree));
        Analytical_Integral = Analytical_Integral - A_S_I; 
    end
    
    boundaryNodes(nNodes).AddSmallFeatureList(SmallFeatureList);
    
end
Analytical_Integral_f_Integrand = [10 0.1 0.4 -1 5 2 9 -10 10 -10]*Analytical_Integral;
nSamplingPts = 1000;
for nT = 1 : 4
Error_Actual = zeros(nE,1);
Error_Predicted = zeros( nE, 1);
Integral_f = zeros( nE, 1);
    
for ie = 1 : nE
close all;
s = sprintf('nT = %d, iE = %d',nT,ie);
disp(s)
nPt = Cell_2D.GetNumberOfPoints_2D(polynomialDegree);
figure;
figure;
x = StartX:0.01:EndX;
y1 =  fy(x);
XY = [x' y1'];
hold on;
plot(XY(:,1),XY(:,2),'b');
hold off;

qTree.DrawTree('red');
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
    
    %node.m_SmallFeatureList = [];
    if( nT == 4 )
        node.m_SmallFeatureList = [];
    else
        for s = 1 : length(node.m_SmallFeatureList)
            node.m_SmallFeatureList{s}.Polygonize(nEdges(ie));
            node.m_SmallFeatureList{s}.DisplayPolygon('g');
            node.m_SmallFeatureList{s}.Display('b');
            node.m_SmallFeatureList{s}.SetSensitivityType(nT-1);
        end
    end
end

leafNodes = [interiorNodes; boundaryNodes];
domainDegree = 2;

axis tight equal;

%Error in Basis functions


for i = 1 : length(leafNodes)
   
   [I,e,XY] = leafNodes(i).GetIntegral( polynomialDegree, distanceToBoundary, domainDegree, f_integrand  );
   Integral_f(ie) =  Integral_f(ie) + I;
   Error_Predicted(ie) = Error_Predicted(ie) + e;
   hold on;
   if(~isempty(XY))
       if( leafNodes(i).IsBoundary())
        scatter(XY(:,1),XY(:,2),'+','r');
       else
        scatter(XY(:,1),XY(:,2),'+','b');
       end
   end
   hold off;
end


% for i = 1 : length(leafNodes)
%    
%     for z = 1 : 10
%         f_z = monomials{z};
%         Integral_AW(SpaceTreeDepth,z) =  Integral_AW(SpaceTreeDepth,z) + f_z(XY(:,1),XY(:,2))'*W;
%     %    [I,e,XY] = leafNodes(i).GetIntegral( polynomialDegree, distanceToBoundary, domainDegree, f_z  );
%         Integral_AW(SpaceTreeDepth,z) =  Integral_AW(SpaceTreeDepth,z) + I;
%         Integral_Error(SpaceTreeDepth,z) = Integral_Error(SpaceTreeDepth,z) + e;
%         hold on;
%           scatter(XY(:,1),XY(:,2));
%         hold off;
%     end
% end

%qTree.PlotError();

Error_Actual(ie) = 100*abs((Analytical_Integral_f_Integrand - Integral_f(ie))/Analytical_Integral_f_Integrand);
end
Error_Actual
Error_Predicted
filename = 'SmallFeatures_CircularHole.xls';
xlswrite(filename, [nEdges' Error_Actual abs(Error_Predicted)], nT);
%xlswrite('SmallFeatures_CircularHole.xls', Analytical_Integral_f_Integrand, 5);
end

% % Open Excel Automation server
% Excel = actxserver('Excel.Application');
% % Make Excel visible
% Excel.Visible=1;
% % Open Excel file
% Workbook = Excel.Workbooks.Open(filename);
% % Get the list of sheets in the workbook
% Sheets = Excel.ActiveWorkbook.Sheets;
% % Rename the first sheet
% Sheets.Item(1).Name = 'Missing Features';
% Sheets.Item(2).Name = '1st Order SSA';
% Sheets.Item(3).Name = '2nd Order SSA';