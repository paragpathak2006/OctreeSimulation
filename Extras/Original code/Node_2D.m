classdef Node_2D < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        m_bbox;
        m_ChildNodes;
        m_PMC;
        m_PointA;
        m_PointB;
        m_PointC;
        m_PointD;
        m_cellType;
        m_ImplicitFn;
        m_Integral;
        m_Error;
        m_SmallFeatureList;
    end
    
    methods
        function obj = Node_2D( PMC, implicitFn, PointA, PointB, PointC, PointD )
            obj.m_PMC = PMC;
            obj.m_PointA = PointA;
            obj.m_PointB = PointB;
            obj.m_PointC = PointC;
            obj.m_PointD = PointD;
            obj.m_ImplicitFn = implicitFn;
            obj.m_cellType = obj.Classify();
            obj.m_bbox = [obj.m_PointA';obj.m_PointC'];
        end
        
        
        
        
        function Corners = GetCornerPts( obj )
            Corners = [obj.m_PointA'; obj.m_PointB'; obj.m_PointC'; obj.m_PointD';];
        end
        
        function [] = AddSmallFeatureList( obj, SmallFeatureList )
            obj.m_SmallFeatureList = SmallFeatureList;
        end
        
        function type = GetCellType( obj )
            type = obj.m_cellType;
        end
        function flag = IsBoundary( obj ) 
            if( obj.m_cellType == 0 || obj.m_cellType == 15 )
                flag = false;
            else
                flag = true;
            end
        end
        
        function flag = IsInterior( obj ) 
            if( obj.m_cellType == 15 )
                flag = true;
            else
                flag = false;
            end
        end
        
        function flag = IsOuter( obj ) 
            if( obj.m_cellType == 0 )
                flag = true;
            else
                flag = false;
            end
        end
        
%Begin 1---------------------------------------------------------
        function Partition( obj, depth )
            
            if( depth > 0 )
                obj.Partition_One_Node();
                for i = 1 : 4
                    if( obj.m_ChildNodes(i).IsBoundary() )
                        obj.m_ChildNodes(i).Partition( depth -1 );
                    end
                end
            end
        end
        
        function GetLeafNodes( obj, nodeList )
            
            if( isempty( obj.m_ChildNodes ) )
                nodeList.Add(obj);
                return;
            else
                for i = 1 : 4 
                    obj.m_ChildNodes(i).GetLeafNodes( nodeList );
                end
            end
        end
        
%Begin 2--------------------------------------------------------
        function Partition_One_Node( obj )
            PointF = (obj.m_PointA + obj.m_PointC)/2;
            PointE = [PointF(1);obj.m_PointA(2)];
            PointG = [obj.m_PointA(1); PointF(2)];
            PointH = [obj.m_PointB(1); PointF(2)];
            PointI = [PointF(1); obj.m_PointC(2)];
 
            Node1 = Node_2D( obj.m_PMC,obj.m_ImplicitFn, obj.m_PointA, PointE, PointF, PointG);
            Node2 = Node_2D( obj.m_PMC, obj.m_ImplicitFn, PointE, obj.m_PointB, PointH, PointF);
            Node3 = Node_2D( obj.m_PMC, obj.m_ImplicitFn, PointF, PointH, obj.m_PointC, PointI);
            Node4 = Node_2D( obj.m_PMC, obj.m_ImplicitFn, PointG, PointF, PointI, obj.m_PointD);
            obj.m_ChildNodes = [Node1; Node2; Node3; Node4];
        end
        
%Begin 3--------------------------------------------------------
        function cellType = Classify( obj )
            f1 = obj.m_PMC(obj.m_PointA(1),obj.m_PointA(2));
            f2 = bitsll(obj.m_PMC(obj.m_PointB(1),obj.m_PointB(2)),1);
            f3 = bitsll(obj.m_PMC(obj.m_PointC(1), obj.m_PointC(2)),2);
            f4 = bitsll(obj.m_PMC(obj.m_PointD(1),obj.m_PointD(2) ),3);
            
            f12 = bitor(f1,f2,'uint8');
            f34 = bitor(f3,f4,'uint8');
            cellType = bitor(f12,f34,'uint8');
        end
        
        function Bbox = GetBoundingBox( obj, XYPolygon )
            Xmin = min(XYPolygon(:,1));
            Ymin = min(XYPolygon(:,2));
            Xmax = max(XYPolygon(:,1));
            Ymax = max(XYPolygon(:,2));
            
            Bbox = [ Xmin Ymin; Xmax Ymax];
        end
        
        
        
%Begin 4--------------------------------------------------------
        
        function [ApproxPolygon, ComputeSSATerm] = GetApproximatePolygon( obj )
            ApproxPolygon = [];
            ComputeSSATerm = [];
            switch( obj.GetCellType() )
                case 0
                    %Outer Node
                case 1
                    %Interior Node
                    ApproxPolygon = zeros( 4, 2);
                    ComputeSSATerm = zeros( 4, 1);
                    ApproxPolygon(1,:) = obj.m_PointA';
                    ApproxPolygon(2,:) = obj.InterpolateZeroLevelSet( obj.m_PointA, obj.m_PointB );
                    ApproxPolygon(3,:) = obj.InterpolateZeroLevelSet( obj.m_PointA, obj.m_PointD ); 
                    ApproxPolygon(4,:) = ApproxPolygon(1,:);
                    ComputeSSATerm(2) = 1;
                case 2
                    ApproxPolygon = zeros( 4, 2);
                    ComputeSSATerm = zeros( 4, 1);
                    ApproxPolygon(1,:) = obj.InterpolateZeroLevelSet( obj.m_PointB, obj.m_PointA );
                    ApproxPolygon(2,:) = obj.m_PointB';
                    ApproxPolygon(3,:) = obj.InterpolateZeroLevelSet( obj.m_PointB, obj.m_PointC ); 
                    ApproxPolygon(4,:) = ApproxPolygon(1,:);
                    ComputeSSATerm(3) = 1;
                case 3
                    ApproxPolygon = zeros( 5, 2);
                    ComputeSSATerm = zeros( 5, 1);
                    ApproxPolygon(1,:) = obj.m_PointA';
                    ApproxPolygon(2,:) = obj.m_PointB';
                    ApproxPolygon(3,:) = obj.InterpolateZeroLevelSet( obj.m_PointB, obj.m_PointC ); 
                    ApproxPolygon(4,:) = obj.InterpolateZeroLevelSet( obj.m_PointA, obj.m_PointD );
                    ApproxPolygon(5,:) = ApproxPolygon(1,:);
                    ComputeSSATerm(3) = 1;
                case 4
                    ApproxPolygon = zeros( 4, 2);
                    ComputeSSATerm = zeros( 4, 1);
                    ApproxPolygon(1,:) = obj.InterpolateZeroLevelSet( obj.m_PointC, obj.m_PointB );
                    ApproxPolygon(2,:) = obj.m_PointC';
                    ApproxPolygon(3,:) = obj.InterpolateZeroLevelSet( obj.m_PointC, obj.m_PointD ); 
                    ApproxPolygon(4,:) = ApproxPolygon(1,:);
                    ComputeSSATerm(3) = 1;
                case 5
                    disp('Ambiguous Case');
                case 6
                    ApproxPolygon = zeros( 5, 2);
                    ComputeSSATerm = zeros( 5, 1);
                    ApproxPolygon(1,:) = obj.InterpolateZeroLevelSet( obj.m_PointA, obj.m_PointB );
                    ApproxPolygon(2,:) = obj.m_PointB';
                    ApproxPolygon(3,:) = obj.m_PointC'; 
                    ApproxPolygon(4,:) = obj.InterpolateZeroLevelSet( obj.m_PointC, obj.m_PointD );
                    ApproxPolygon(5,:) = ApproxPolygon(1,:);
                    ComputeSSATerm(4) = 1;
                case 7
                    ApproxPolygon = zeros( 6, 2);
                    ComputeSSATerm = zeros( 6, 1);
                    ApproxPolygon(1,:) = obj.m_PointA';
                    ApproxPolygon(2,:) = obj.m_PointB';
                    ApproxPolygon(3,:) = obj.m_PointC'; 
                    ApproxPolygon(4,:) = obj.InterpolateZeroLevelSet( obj.m_PointC, obj.m_PointD );
                    ApproxPolygon(5,:) = obj.InterpolateZeroLevelSet( obj.m_PointA, obj.m_PointD );
                    ApproxPolygon(6,:) = ApproxPolygon(1,:) ;
                    ComputeSSATerm(4) = 1;
                case 8 
                    ApproxPolygon = zeros( 4, 2);
                    ComputeSSATerm = zeros( 4, 1);
                    ApproxPolygon(1,:) = obj.m_PointD';
                    ApproxPolygon(2,:) = obj.InterpolateZeroLevelSet( obj.m_PointD, obj.m_PointA );
                    ApproxPolygon(3,:) = obj.InterpolateZeroLevelSet( obj.m_PointD, obj.m_PointC ); 
                    ApproxPolygon(4,:) = ApproxPolygon(1,:);
                    ComputeSSATerm(2) = 1;
                case 9
                    ApproxPolygon = zeros( 5, 2);
                    ComputeSSATerm = zeros( 5, 1);
                    ApproxPolygon(1,:) = obj.m_PointA';
                    ApproxPolygon(2,:) = obj.InterpolateZeroLevelSet( obj.m_PointA, obj.m_PointB );
                    ApproxPolygon(3,:) = obj.InterpolateZeroLevelSet( obj.m_PointC, obj.m_PointD ); 
                    ApproxPolygon(4,:) = obj.m_PointD';
                    ApproxPolygon(5,:) = ApproxPolygon(1,:);
                    ComputeSSATerm(2) = 1;

                case 10 
                    disp('Ambiguous Case');
                case 11
                    ApproxPolygon = zeros( 6, 2);
                    ComputeSSATerm = zeros( 6, 1);
                    ApproxPolygon(1,:) = obj.m_PointA';
                    ApproxPolygon(2,:) = obj.m_PointB';
                    ApproxPolygon(3,:) = obj.InterpolateZeroLevelSet( obj.m_PointB, obj.m_PointC ); 
                    ApproxPolygon(4,:) = obj.InterpolateZeroLevelSet( obj.m_PointD, obj.m_PointC );
                    ApproxPolygon(5,:) = obj.m_PointD';
                    ApproxPolygon(6,:) = ApproxPolygon(1,:);
                    ComputeSSATerm(3) = 1;
                case 12
                    ApproxPolygon = zeros( 5, 2);
                    ComputeSSATerm = zeros( 5, 1);
                    ApproxPolygon(1,:) = obj.InterpolateZeroLevelSet( obj.m_PointD, obj.m_PointA );
                    ApproxPolygon(2,:) = obj.InterpolateZeroLevelSet( obj.m_PointC, obj.m_PointB );
                    ApproxPolygon(3,:) = obj.m_PointC'; 
                    ApproxPolygon(4,:) = obj.m_PointD';
                    ApproxPolygon(5,:) = ApproxPolygon(1,:);
                    ComputeSSATerm(1) = 1;
                case 13
                    ApproxPolygon = zeros( 6, 2);
                    ComputeSSATerm = zeros( 6, 1);
                    ApproxPolygon(1,:) = obj.m_PointA';
                    ApproxPolygon(2,:) = obj.InterpolateZeroLevelSet( obj.m_PointA, obj.m_PointB );
                    ApproxPolygon(3,:) = obj.InterpolateZeroLevelSet( obj.m_PointC, obj.m_PointB ); 
                    ApproxPolygon(4,:) = obj.m_PointC';
                    ApproxPolygon(5,:) = obj.m_PointD';
                    ApproxPolygon(6,:) = ApproxPolygon(1,:);
                    ComputeSSATerm(2) = 1;
                case 14
                    ApproxPolygon = zeros( 6, 2);
                    ComputeSSATerm = zeros( 6, 1);
                    ApproxPolygon(1,:) = obj.InterpolateZeroLevelSet( obj.m_PointB, obj.m_PointA );
                    ApproxPolygon(2,:) = obj.m_PointB';
                    ApproxPolygon(3,:) = obj.m_PointC';
                    ApproxPolygon(4,:) = obj.m_PointD';
                    ApproxPolygon(5,:) = obj.InterpolateZeroLevelSet( obj.m_PointD, obj.m_PointA ); 
                    ApproxPolygon(6,:) = ApproxPolygon(1,:);
                    ComputeSSATerm(5) = 1;
                case 15
                    %Outer Node
                otherwise
                    disp('Error : Cell Type Out of Range');
            end
        end
        
%Begin 5--------------------------------------------------------
        function Pts = InterpolateZeroLevelSet( obj, ptA, ptB )
            ax = ptA(1);
            ay = ptA(2);
            d = ptB - ptA;
            dx = d(1);
            dy = d(2);
            fun = @(t) obj.m_ImplicitFn( ax + t*dx , ay + t*dy );
            f1 = abs(fun(0));
            f2 = abs(fun(1));
            
            %t0 = f1/(f1 + f2);
            t0 = 0.5;
            t = fsolve(fun,t0);
            %fMin = @(t) fun(t0).*fun(t0);
            %t = fminsearch(fMin,0);
            if( t >= 0 )
                Pts = [ax + t*dx; ay + t*dy];
            else
                Pts = ptA';
            end
        end
        
        function Draw( obj, color )
            obj.Draw_One_Node( color );
            
            if(~isempty(obj.m_ChildNodes) )
                for i = 1: 4
                    obj.m_ChildNodes(i).Draw(color);
                end
            end
            
        end
        function Draw_One_Node(obj, color)
            diff = (obj.m_PointC - obj.m_PointA);
            w = abs(diff(1));
            h = abs(diff(2));
            hold on;
            %rectangle('Position',[obj.m_PointA(1),obj.m_PointB(2),w,h],'EdgeColor',color);
            p1 = obj.m_PointA;
            p2 = obj.m_PointB;
            p3 = obj.m_PointC;
            p4 = obj.m_PointD;
            p5 = p1;
            p = [p1';p2';p3';p4';p5'];
            line(p(:,1),p(:,2),'Color',color);
            hold off;
        end
       
    end
    
    methods(Static)
        function DrawPolygon( PolygonCornerXY, color )
            line( PolygonCornerXY(:,1), PolygonCornerXY(:,2), 'Color', color, 'LineWidth',2 );
        end
        
%         function [Pts,W] = GetCartesianQuadrature_Box( CornerPtsXY, nPt )
%            
%            PtA = CornerPtsXY(1,:);
%            PtB = CornerPtsXY(2,:);
%            PtC = CornerPtsXY(3,:);
%            PtD = CornerPtsXY(4,:);
%            
%            a = PtA(1);
%            b = PtB(1);
%            c = PtA(2);
%            d = PtC(2);
%            
%            [Wx,X] = Cell_2D.Gauss_W_X_Scaled( a , b, nPt );
%            [Wy,Y] = Cell_2D.Gauss_W_X_Scaled( c , d, nPt );
%            
%            Pts = zeros(nPt*nPt,2);
%            W = zeros(nPt*nPt,1);
%            
%            count = 1;
%            for i = 1 : nPt
%                thisY = Y(i);
%                thisWy = Wy(i);
%                for j = 1 : nPt
%                    thisX = X(j);
%                    thisWx = Wx(j);
%                    Pts(count,:) = [thisX thisY];
%                    W(count) = thisWx * thisWy;
%                    count = count + 1;
%                end
%            end
%        end
       
    end
    
end

%          function Pts = GetRandomQuadraturePts_Boundary( obj, nPt, XYPolygon, nSamplingPts )
%             
%             XPoly = XYPolygon(:,1);
%             YPoly = XYPolygon(:,2);
%             
%             n = length(XPoly);
%             
%             
%             SamplingPts_PerEdge = ceil(nSamplingPts);
%             TestPts = [];
%             for i = 1 : n-1
%                 XY1 = [XPoly(i); YPoly(i)];
%                 XY2 = [XPoly(i+1); YPoly(i+1)];
%                 L = Cell_2D.GetEdgeLength(XY1,XY2);
%                 normal = Cell_2D.GetEdgeNormal(XY1,XY2);
%                 t = L / (SamplingPts_PerEdge-1);
%                 p = 0:t:L;
%                 Origin = repmat(XY1',length(p),1);
%                 n = repmat(normal,1,length(p));
%                 TestPts = [TestPts ;Origin + (n*diag(p))'];
%             end
%             XTestPts = TestPts(:,1); 
%             YTestPts = TestPts(:,2); 
% 
%             fVal = obj.m_PMC(XTestPts, YTestPts);
%             [m, n] = size( XTestPts);
%             X= reshape(XTestPts,[m*n,1]);
%             Y = reshape(YTestPts,[m*n,1]);
%             fVal = reshape(fVal,[m*n,1]);
% 
%             PassPtsX = X(fVal >0);
%             PassPtsY = Y(fVal >0);
% 
% %             [Index] = inpolygon( PassPtsX, PassPtsY, XPoly, YPoly);
% % 
% %             PassPtsX = PassPtsX(Index);
% %             PassPtsY = PassPtsY(Index);
% 
%             r = randperm(size(PassPtsX,1));
% 
%             PassPtsX = PassPtsX(r(1:nPt));
%             PassPtsY = PassPtsY(r(1:nPt));
%             Pts = [PassPtsX PassPtsY];
%          end
%         
%         function Pts = GetRandomQuadraturePts( obj, nPt, XYPolygon, nSamplingPts )
%             
%             XPoly = XYPolygon(:,1);
%             YPoly = XYPolygon(:,2);
%             BBox = obj.GetBoundingBox( XYPolygon );
%             StartPt = BBox(1,:);
%             EndPt = BBox(2,:);
%             hx = (EndPt - StartPt)/(nSamplingPts+1);
%             XTestPts = StartPt+hx:hx:EndPt-hx;
%             minY = obj.m_PointB(2);
%             maxY = obj.m_PointC(2);
%             hy = (maxY - minY)/(nSamplingPts+1);
%             YTestPts = minY + hy : hy : maxY - hy;
% 
%             [XTestPts, YTestPts] = meshgrid(XTestPts,YTestPts); 
% 
%             fVal = obj.m_PMC(XTestPts, YTestPts);
%             
%             nSmallFeatures = size(obj.m_SmallFeatureList,1);
% 
%             for i = 1 : nSmallFeatures
%                 if( i == 1 )
%                     fVal_Small_Features = ~obj.m_SmallFeatureList{i}.PMC(XTestPts,YTestPts);
%                 else
%                     fVal_Small_Features = and( fVal_Small_Features ,  ~obj.m_SmallFeatureList{i}.PMC(XTestPts,YTestPts));
%                 end
%                 
%             end
%             
%             [m, n] = size( XTestPts);
%             X= reshape(XTestPts,[m*n,1]);
%             Y = reshape(YTestPts,[m*n,1]);
%             fVal = reshape(fVal,[m*n,1]);
%             if( nSmallFeatures > 0 )
%                 fVal_Small_Features = reshape(fVal_Small_Features,[m*n,1]);
%                 fVal = and( fVal, fVal_Small_Features);
%             end
%             PassPtsX = X(fVal >0);
%             PassPtsY = Y(fVal >0);
% 
%             [Index] = inpolygon( PassPtsX, PassPtsY, XPoly, YPoly);
% 
%             PassPtsX = PassPtsX(Index);
%             PassPtsY = PassPtsY(Index);
% 
%             r = randperm(size(PassPtsX,1));
% 
%             PassPtsX = PassPtsX(r(1:nPt));
%             PassPtsY = PassPtsY(r(1:nPt));
%             Pts = [PassPtsX PassPtsY];
%         end
%         
%         function [XY,W] = GetQuadrature_AW( obj, polynomialDegree, distanceToBoundary, domainDegree  )
%             nPt = Cell_2D.GetNumberOfPoints_2D(polynomialDegree);
%             if( obj.IsInterior() )
%                 CornerPtsXY = [obj.m_PointA';obj.m_PointB'; obj.m_PointC'; obj.m_PointD'];
%                 numPts = ceil(sqrt(nPt));
%                 [XY,W] = Node_2D.GetCartesianQuadrature_Box( CornerPtsXY, numPts );
%             elseif( obj.IsBoundary() )
%                 nSamplingPts = 20*nPt;
%                 [ApproxPolygon, ComputeSSATerm] = obj.GetApproximatePolygon();
%                 XY = obj.GetRandomQuadraturePts( nPt, ApproxPolygon, nSamplingPts );
%                 thisCell = Cell_2D(ApproxPolygon, distanceToBoundary, domainDegree, true);
%                 W = thisCell.GetAdaptiveWeights( XY, polynomialDegree, ComputeSSATerm );
%                 
%                 nSmallFeatures = length(obj.m_SmallFeatureList);
%                 
%                 for i = 1 : nSmallFeatures
%                     [X,Y] = obj.m_SmallFeatureList{i}.Polygonize();
%                     ApproxPolygon = [X Y];
%                     ComputeSSATerm = obj.m_SmallFeatureList{i}.GetSSACorrection();
%                     dist = @(QP,Direction) obj.m_SmallFeatureList{i}.DistanceToBoundary(QP,Direction);
%                     AdaptiveQuadrature = true;
%                     domainDegree = obj.m_SmallFeatureList{i}.GetDomainDegree();
%                     thisCell = Cell_2D(ApproxPolygon, dist, domainDegree, AdaptiveQuadrature);
%                     W_f = thisCell.GetAdaptiveWeights( XY, polynomialDegree, ComputeSSATerm );
%                     W = W - W_f;
%                 end
%             else
%                 XY = [];
%                 W = [];
%             end
%         end
%         
%        
%          function [I, e, XY] = GetIntegral( obj, polynomialDegree, distanceToBoundary, domainDegree, integrand  )
%             nPt = Cell_2D.GetNumberOfPoints_2D(polynomialDegree);
%             if( obj.IsInterior() )
%                 CornerPtsXY = [obj.m_PointA';obj.m_PointB'; obj.m_PointC'; obj.m_PointD'];
%                 numPts = ceil(sqrt(nPt));
%                 [XY,W] = Node_2D.GetCartesianQuadrature_Box( CornerPtsXY, numPts );
%                 I = integrand(XY(:,1),XY(:,2))'*W;
%                 e = 0;
%             elseif( obj.IsBoundary() )
%                 nSamplingPts = 10*nPt;
%                 [ApproxPolygon, ComputeSSATerm] = obj.GetApproximatePolygon();
%                 polygonArea = abs(polyarea(ApproxPolygon(1:end-1,1), ApproxPolygon(1:end-1,2))) ;
%                 if( polygonArea < 1e-5)
%                     e = 0;
%                     I = polygonArea;
%                     XY = [];
%                 else
%                 XY = obj.GetRandomQuadraturePts( nPt, ApproxPolygon, nSamplingPts );
%                 %XY = obj.GetRandomQuadraturePts_Boundary(nPt, ApproxPolygon, nSamplingPts );
%                 thisCell = Cell_2D(ApproxPolygon, distanceToBoundary, domainDegree, true);
%                 [I, e] = thisCell.GetIntegral( XY, polynomialDegree, ComputeSSATerm, integrand, obj.m_SmallFeatureList );
%                 end
%             else
%                 I = 0;
%                 e = 0;
%             end
%             obj.m_Integral = I;
%             obj.m_Error = e;
%         end
%         
