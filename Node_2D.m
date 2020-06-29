

classdef Node_2D < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        m_bbox;
        m_ChildNodes;
        m_PMC;
        m_ImplicitFn;
        
        m_PointA;
        m_PointB;
        m_PointC;
        m_PointD;
        
        m_PointE;
        m_PointF;
        m_PointG;
        m_PointH;
        
        m_cellType;
        m_Integral;
        m_Error;
        m_SmallFeatureList;
        P;
        numTri;
    end
    
    methods
        function obj = Node_2D( PMC, implicitFn, PointA, PointB, PointC, PointD, PointE, PointF, PointG, PointH)
            
            obj.m_PMC = PMC;
            obj.m_ImplicitFn = implicitFn;
            
            obj.m_PointA = PointA;
            obj.m_PointB = PointB;
            obj.m_PointC = PointC;
            obj.m_PointD = PointD;
            
            obj.m_PointE = PointE;
            obj.m_PointF = PointF;
            obj.m_PointG = PointG;
            obj.m_PointH = PointH;
            
            %             obj.p = zeros(28);
            obj.P = zeros(28,3);
            
            obj.Classify();
            obj.m_bbox = [obj.m_PointA';obj.m_PointG'];
        end
        
        %         function Corners = GetCornerPts( obj )
        %             Corners = [obj.m_PointA'; obj.m_PointB'; obj.m_PointC'; obj.m_PointD';];
        %         end
        %
        %         function [] = AddSmallFeatureList( obj, SmallFeatureList )
        %             obj.m_SmallFeatureList = SmallFeatureList;
        %         end
        
        function type = GetCellType( obj )
            type = obj.m_cellType;
        end
        
        function flag = IsBoundary( obj )
            if( obj.m_cellType == 0 || obj.m_cellType == 255 )
                flag = false;
            else
                flag = true;
            end
        end
        
        function flag = IsInterior( obj )
            if( obj.m_cellType == 255 )
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
            global Cube;
            if( depth > 0 )
                obj.Partition_One_Node();
                for i = 1 : 8
                    if( obj.m_ChildNodes(i).IsBoundary() )
                        obj.m_ChildNodes(i).Partition( depth -1 );
                    end
                end
            else
                cellType = obj.m_cellType;
                isAmbigous = Cube(cellType*5+1,14);
                if(isAmbigous) obj.Partition_One_Node();end
            end
        end
        
        function GetLeafNodes( obj, nodeList )
            
            if( isempty( obj.m_ChildNodes ) )
                nodeList.Add(obj);
                return;
            else
                for i = 1 : 8
                    obj.m_ChildNodes(i).GetLeafNodes( nodeList );
                end
            end
        end
        
        %Begin 2--------------------------------------------------------
        function Partition_One_Node( obj )
            
            
            PointAB = 1.0*(obj.m_PointA + obj.m_PointB)/2.0;
            PointBC = 1.0*(obj.m_PointB + obj.m_PointC)/2.0;
            PointCD = 1.0*(obj.m_PointC + obj.m_PointD)/2.0;
            PointAD = 1.0*(obj.m_PointD + obj.m_PointA)/2.0;
            
            PointEF = 1.0*(obj.m_PointE +obj.m_PointF)/2.0;
            PointFG = 1.0*(obj.m_PointF +obj.m_PointG)/2.0;
            PointGH = 1.0*(obj.m_PointG +obj.m_PointH)/2.0;
            PointEH = 1.0*(obj.m_PointH +obj.m_PointE)/2.0;
            
            PointAE = 1.0*(obj.m_PointA +obj.m_PointE)/2.0;
            PointBF = 1.0*(obj.m_PointB +obj.m_PointF)/2.0;
            PointCG = 1.0*(obj.m_PointC +obj.m_PointG)/2.0;
            PointDH = 1.0*(obj.m_PointD +obj.m_PointH)/2.0;
            
            PointABCD = (PointAB +PointCD)/2;
            PointEFGH = (PointEF +PointGH)/2;
            
            PointABEF = (PointAB +PointEF)/2;
            PointCDGH = (PointCD +PointGH)/2;
            
            PointADEH = (PointAD +PointEH)/2;
            PointBCFG = (PointBC +PointFG)/2;
            
            PointCET = (PointABCD +PointEFGH)/2;
            
            %-----------------------------
            Node1 = Node_2D( obj.m_PMC,obj.m_ImplicitFn, obj.m_PointA,PointAB,PointABCD,PointAD,PointAE,PointABEF,PointCET,PointADEH);
            Node2 = Node_2D( obj.m_PMC,obj.m_ImplicitFn, PointAB,obj.m_PointB,PointBC,PointABCD,PointABEF,PointBF,PointBCFG,PointCET);
            Node3 = Node_2D( obj.m_PMC,obj.m_ImplicitFn, PointABCD,PointBC,obj.m_PointC,PointCD,PointCET,PointBCFG,PointCG,PointCDGH);
            Node4 = Node_2D( obj.m_PMC,obj.m_ImplicitFn, PointAD,PointABCD,PointCD,obj.m_PointD,PointADEH,PointCET,PointCDGH,PointDH);
            %-----------------------------
            Node5 = Node_2D( obj.m_PMC,obj.m_ImplicitFn, PointAE,PointABEF,PointCET,PointADEH,obj.m_PointE,PointEF,PointEFGH,PointEH);
            Node6 = Node_2D( obj.m_PMC,obj.m_ImplicitFn, PointABEF,PointBF,PointBCFG,PointCET,PointEF,obj.m_PointF,PointFG,PointEFGH);
            Node7 = Node_2D( obj.m_PMC,obj.m_ImplicitFn, PointCET,PointBCFG,PointCG,PointCDGH,PointEFGH,PointFG,obj.m_PointG,PointGH);
            Node8 = Node_2D( obj.m_PMC,obj.m_ImplicitFn, PointADEH,PointCET,PointCDGH,PointDH,PointEH,PointEFGH,PointGH,obj.m_PointH);
            
            obj.m_ChildNodes = [Node1; Node2; Node3; Node4;Node5; Node6; Node7; Node8];
        end
        
        %Begin 3--------------------------------------------------------
        function cellType = Classify( obj )
            %              f1 = obj.m_PMC(obj.m_PointA(1),obj.m_PointA(2),obj.m_PointA(3));
            %              f2 = bitsll(obj.m_PMC(obj.m_PointB(1),obj.m_PointB(2),obj.m_PointB(3)),1);
            %              f3 = bitsll(obj.m_PMC(obj.m_PointC(1), obj.m_PointC(2),obj.m_PointC(3)),2);
            %              f4 = bitsll(obj.m_PMC(obj.m_PointD(1),obj.m_PointD(2),obj.m_PointD(3) ),3);
            %
            %              f5 = bitsll(obj.m_PMC(obj.m_PointE(1),obj.m_PointE(2) ,obj.m_PointE(3)),4);
            %              f6 = bitsll(obj.m_PMC(obj.m_PointF(1),obj.m_PointF(2),obj.m_PointF(3) ),5);
            %              f7 = bitsll(obj.m_PMC(obj.m_PointG(1),obj.m_PointG(2) ,obj.m_PointG(3)),6);
            %              f8 = bitsll(obj.m_PMC(obj.m_PointH(1),obj.m_PointH(2) ,obj.m_PointH(3)),7);
            %
            %              f12 = bitor(f1,f2,'uint8');
            %              f34 = bitor(f3,f4,'uint8');
            %
            %             cellType = bitor(f12,f34,'uint8');
            
            f1 = obj.m_PMC(obj.m_PointA(1),obj.m_PointA(2),obj.m_PointA(3));
            f2 = obj.m_PMC(obj.m_PointB(1),obj.m_PointB(2),obj.m_PointB(3))*2;
            f3 = obj.m_PMC(obj.m_PointC(1), obj.m_PointC(2),obj.m_PointC(3))*4;
            f4 = obj.m_PMC(obj.m_PointD(1),obj.m_PointD(2),obj.m_PointD(3))*8;
            
            f5 = obj.m_PMC(obj.m_PointE(1),obj.m_PointE(2) ,obj.m_PointE(3))*16;
            f6 = obj.m_PMC(obj.m_PointF(1),obj.m_PointF(2),obj.m_PointF(3))*32;
            f7 = obj.m_PMC(obj.m_PointG(1),obj.m_PointG(2) ,obj.m_PointG(3))*64;
            f8 = obj.m_PMC(obj.m_PointH(1),obj.m_PointH(2) ,obj.m_PointH(3))*128;
            
            obj.m_cellType = f1 + f2 + f3 + f4 + f5 + f6 +f7 + f8;
            %             obj.numTri = Cube(cellType+1,1);
            %             obj.p(1) = obj.numTri ;
            %             for i = 2 : obj.numTri*3 + 1
            %                 obj.p(i) = Cube(cellType+1,i);
            %             end
        end
        
        function Pts = Intersection(obj, ptA, ptB )
            ax = ptA(1);            ay = ptA(2);            az = ptA(3);
            
            d = ptB - ptA;
            dx = d(1);            dy = d(2);            dz = d(3);
            
            fun = @(t) obj.m_ImplicitFn( ax + t*dx , ay + t*dy , az + t*dz );
            
            %             f1 = abs(fun(0));
            %             f2 = abs(fun(1));
            
            %t0 = f1/(f1 + f2);
            t0 = 0.5;
            t = fsolve(fun,t0);
            %fMin = @(t) fun(t0).*fun(t0);
            %t = fminsearch(fMin,0);
            if( t >= 0&&t<=1)
                Pts = [ax + t*dx; ay + t*dy;az + t*dz ];
            else
                if(t<0)Pts = ptA';end
                if(t>1)Pts = ptB';end
            end
        end
        % % % % %
        function [IntegrationPts] = AllocateIntegrationPoints_Random(obj, nPt,nsolid)
            startVertex = obj.m_PointA;
            endVertex = obj.m_PointG;
            
            [Vertices, Faces,numVertices, numFaces,isTrim ] = obj.getVerticesFaces(nsolid);
            
            totalNumberOfPoints = nPt*nPt*nPt;
            
            delta = (endVertex - startVertex)/12;
            delta = (endVertex - startVertex)/5;

            dx = delta(1);
            dy = delta(2);
            dz = delta(3);
            
            x = (startVertex(1)+dx):dx:(endVertex(1)-dx);
            y = (startVertex(2)+dy):dy:(endVertex(2)-dy);
            z = (startVertex(3)+dz):dz:(endVertex(3)-dz);

%             x = (startVertex(1)):dx:(endVertex(1));
%             y = (startVertex(2)):dy:(endVertex(2));
%             z = (startVertex(3)):dz:(endVertex(3));

            [X,Y,Z] = meshgrid(x,y,z);
            X = reshape(X,[],1);
            Y = reshape(Y,[],1);
            Z = reshape(Z,[],1);
            Pts = [X Y Z];
             inFlags = obj.m_PMC(Pts(:,1),Pts(:,2),Pts(:,3));
             ValidPts = Pts(inFlags,:);
             inFlags = inpolyhedron(Faces,Vertices,ValidPts);
             ValidPts = ValidPts(inFlags,:);
%             %IntegrationPts = ValidPts(1:totalNumberOfPoints,:);
%             ids = zeros(totalNumberOfPoints,1);
%             totalValidPoints = size(ValidPts,1);
%             for i = 1 : totalNumberOfPoints
%                 while(true)
%                     id = randi([1 totalValidPoints],1,1);
%                     if( ismember(ids,id) )
%                         continue;
%                     else
%                         ids(i) = id;
%                         break;
%                     end
%                 end
%             end
% %             IntegrationPts = ValidPts(ids,:);
            IntegrationPts = ValidPts';
        end
        
        function [numSolid, solid] = getSolids(obj)
            global Cube;
            
            cellType = obj.m_cellType;
            isConjugate = Cube(cellType*5+1,13);
            numSolid = Cube(cellType*5+1,15);
            if(numSolid==0) return;end
            if(isConjugate==0)
                for k = 1:numSolid-1
                    solid(k) = k;
                end
                numSolid = numSolid-1;
            else
                solid(1) = numSolid;
                numSolid = 1;
            end
        end
        %         function samplePointsOnBoundary(obj,n)
        %
        %             global Cube;
        %             global FV;
        %             k=0;
        %             for i = 0:10
        %             for j = 0:10-i
        %                 k++;
        %                 [x(k);y(k);z(k)] = [i/10;j/10;(10-i-j)/10];
        %             end
        %             end
        %
        %             cellType = obj.m_cellType;
        %             numSolids=Cube(cellType*5+1,15);
        %             if(numSolids==0) return;end
        %             isConjugate = Cube(cellType*5+1,13);
        %
        %             for k = 1:numSolids-1
        %                 if(isConjugate) k = numSolids;end
        %                 solid = 5*cellType+k;
        %
        %                 numpts = FV(solid,1);       V=zeros(numpts,3);
        %                 numface = FV(solid,2);      F=zeros(numface,3);
        %
        %                 for i=1:numface/3
        %
        %                     fIndex = FV(solid,18+3*i);     pIndex = FV(solid,2+fIndex);  p1 = P(pIndex); if(pIndex>15) continue;end
        %                     fIndex = FV(solid,18+3*i+1);   pIndex = FV(solid,2+fIndex);  p2 = P(pIndex); if(pIndex>15) continue;end
        %                     fIndex = FV(solid,18+3*i+2);   pIndex = FV(solid,2+fIndex);  p3 = P(pIndex); if(pIndex>15) continue;end
        %
        %                     pts(k) = p1*x(k) + p2*y(k) + p3*z(k);
        %                     n(k) = normal(p1,p2,p3);
        %                 end
        %             end
        %
        %             pin = insidePolyhedron(obj,pts);
        %             j=1;
        %             for i =1:n
        %                 if(pin(i))
        %                     pout(j,:) = pt(i,:);
        %                     nout(j,:) = n(i,:);j++;
        %                 end
        %             end
        %             return pout,nout;
        %         end
        %
        %         function  p = normal(p1,p2,p3)
        %             a=p2-p1;b=p3-p1;
        %             p(1) = a(2)*b(3) - a(3)*b(2);
        %             p(2) = a(3)*b(1) - a(1)*b(3);
        %             p(3) = a(1)*b(2) - a(2)*b(1);
        %             return p;
        %         end
        % % % % % % % % % %
        %         function [IntegrationPts] = AllocateIntegrationPoints_Random_Boundary( PolyhedraCoords, nPt, EmbeddedDomain, Vertices, Faces )
        %
        %             totalNumberOfPoints = nPt*nPt*nPt;
        %
        %             nPolygons = size(PolyhedraCoords,1);
        %
        %             nPointsPerPolygon = ceil((totalNumberOfPoints/nPolygons));
        %             IntegrationPts = [];
        %             for i = 1 : nPolygons
        %                 ThreeDPolygonPoints = PolyhedraCoords{i};
        %                 [TwoDPolyPoints,dir,startVertex,endVertex,Origin] = Helper.Shit2Dto3D( ThreeDPolygonPoints );
        %
        %                 [Pts, Pts3D] = Helper.SamplePoints2D( startVertex,endVertex,dir, 12);
        %                 [m,n] = size(Pts3D);
        %                 OriginArray = repmat(Origin,m,1);
        %                 Pts3D = Pts3D + OriginArray;
        %                 %Pts3D = Pts3D(inpolyhedron(Faces,Vertices,Pts3D),:);
        %                 inFlags = EmbeddedDomain.getDomainIndex_vectorized(Pts3D);
        %                 Pts3D = Pts3D(inFlags,:);
        %                 ValidPts = Pts(inFlags,:);
        %
        %                 inFlags = inpolygon(ValidPts(:,1),ValidPts(:,2),TwoDPolyPoints(:,1),TwoDPolyPoints(:,2));
        %                 ValidPts = Pts3D(inFlags,:);
        %                 %IntegrationPts = ValidPts(1:totalNumberOfPoints,:);
        %                 ids = zeros(nPointsPerPolygon,1);
        %                 totalValidPoints = size(ValidPts,1);
        %
        %                 for j = 1 : nPointsPerPolygon
        %                     while(true)
        %                         id = randi([1 totalValidPoints],1,1);
        %                         if( ismember(ids,id) )
        %                             continue;
        %                         else
        %                             ids(j) = id;
        %                             break;
        %                         end
        %                     end
        %                 end
        %                 IntegrationPts = [IntegrationPts;ValidPts(ids,:)];
        %             end
        %
        %         end
        
        %         function samplePointsInsidePolyhedra(obj,n)
        %
        %             m = 10;
        %
        %             xMin = m_PointA(1); xMax = m_PointG(1);
        %             yMin = m_PointA(2); yMax = m_PointG(2);
        %             zMin = m_PointA(3); zMax = m_PointG(3);
        %
        %             xp = (xMax + xMin)/2;     Lx = (xMax - xMin)/2;
        %             yp = (yMax + yMin)/2;     Ly = (yMax - yMin)/2;
        %             zp = (zMax + zMin)/2;     Lz = (zMax - zMin)/2;
        %
        %             x=xMin:Lx/m:xMax;
        %             y=yMin:Ly/m:yMax;
        %             z=zMin:Lz/m:zMax;
        %
        % %             for i = -m:m
        % %                x(m+1+i) = xp + i*Lx/m;
        % %                y(m+1+i) = yp + i*Ly/m;
        % %                z(m+1+i) = zp + i*Lz/m;
        % %             end
        %
        %             [X,Y,Z] = meshgrid(x,y,z);
        %             pt= [X(i);Y(j);Z(k)];    %generate array of points
        %             pin = insidePolyhedron(obj,pts);
        %             j=1;
        %             for i =1:n
        %                 if(pin(i)) pout(j,:) = pt(i,:);j++;end
        %             end
        %
        %             return pout;
        %         end
        
        function [dist , iPt] = distance(obj,pt,d)
            ax = pt(1);            ay = pt(2);            az = pt(3);
            dx = d(1);            dy = d(2);            dz = d(3);
            
            fun = @(t) obj.m_ImplicitFn( ax + t*dx , ay + t*dy , az + t*dz );
            t0 = 0.5;
            t = fsolve(fun,t0);
            iPt = [ax + t*dx; ay + t*dy;az + t*dz ];
            dist  = norm(pt-iPt);
        end
        
        function [V,F,numVertices,numfaces,isTrim] = getVerticesFaces(obj,nSolid)
            %             global Cube;
            global FV;
            
            cellType = obj.m_cellType;
            %             numSolids=Cube(cellType*5+1,15);
            %             if(numSolids==0) return;end
            %             for k = 1:numSolids-1
            %                 if(isConjugate) k = numSolids;end
            solid = 5*cellType + nSolid;
            numVertices = FV(solid,1);       V=zeros(numVertices,3);
            numfaces = FV(solid,2);      %F = zeros(numfaces*3)';

%             nSolid
%             solid
%             numVertices 
%             numfaces 
            
            for i=1:numVertices
                pIndex = FV(solid,2+i);
%                 pIndex
                
                V(i,:)= obj.P(pIndex,:)';
                if(pIndex>15) isTrim(i)=0;end
                if(pIndex<15)isTrim(i)=1;end
                    %                 pIndex
%                 obj.P(1,:)
%                 obj.P(2,:)
%                 obj.P(3,:)
%                 obj.P(4,:)
% 
%                 obj.P(5,:)
%                 obj.P(6,:)
%                 obj.P(7,:)
%                 obj.P(8,:)
% 
%                 obj.P(9,:)
%                 obj.P(10,:)
%                 obj.P(11,:)
%                 obj.P(12,:)
% 
%                 pIndex
%                 obj.P(:,1)
%                 obj.P(:,2)
%                 obj.P(:,3)
%                 obj.P(:,4)
% 
%                 obj.P(:,5)
%                 obj.P(:,6)
%                 obj.P(:,7)
%                 obj.P(:,8)
% 
%                 obj.P(:,9)
%                 obj.P(:,10)
%                 obj.P(:,11)
%                 obj.P(:,12)

                
            end
%             V
            for i=1:numfaces*3
                fIndex = FV(solid,18+i);
%                  fIndex
                F(i) = fIndex;
            end
            %             end
            
            
        end
        
        %         function [V,F] = insidePolyhedron(obj)
        %             global Cube;
        %             global FV;
        %             pout = zeros(size(pts(1,:)));
        %             pin = pout;
        %             cellType = obj.m_cellType;
        %             numSolids=Cube(cellType*5+1,15);
        %             if(numSolids==0) return;end
        %             isConjugate = Cube(cellType*5+1,13);
        %             for k = 1:numSolids-1
        %                 if(isConjugate) k = numSolids;end
        %                 solid = 5*cellType+k;
        %
        %                 numpts = FV(solid,1);       V=zeros(numpts,3);
        %                 numface = FV(solid,2);      F=zeros(numface,3);
        %
        %                 for i=1:numpts
        %                     pIndex = FV(solid,2+i)
        %                     V(i,:)= P(pIndex);
        %                 end
        %
        %                 for i=1:numface
        %                     fIndex = FV(solid,18+i);
        %                     F(i,:) = P(fIndex);
        %                 end
        %                 pin = Inpolyhedron(V,F,pts);
        %             end
        %
        %             for i=1:size(pts(1,:))
        %                 if(pin(i)==1)pout(i)=1;end
        %             end
        %
        %         end
        %
        function Polygonise( obj)
            global Cube;
            Type=obj.m_cellType;
            
            if(Cube(Type*5+1 ,1)) obj.P(1,:)= obj.Intersection(obj.m_PointA,obj.m_PointB);end
            if(Cube(Type*5+1 ,2)) obj.P(2,:) = obj.Intersection(obj.m_PointB,obj.m_PointC); end
            if(Cube(Type*5+1 ,3)) obj.P(3,:)  = obj.Intersection(obj.m_PointC,obj.m_PointD); end
            if(Cube(Type*5+1 ,4)) obj.P(4,:)  = obj.Intersection(obj.m_PointA,obj.m_PointD); end
            if(Cube(Type*5+1 ,5)) obj.P(5,:)  = obj.Intersection(obj.m_PointE,obj.m_PointF); end
            if(Cube(Type*5+1 ,6)) obj.P(6,:)  = obj.Intersection(obj.m_PointF,obj.m_PointG); end
            if(Cube(Type*5+1 ,7)) obj.P(7,:)  = obj.Intersection(obj.m_PointG,obj.m_PointH); end
            if(Cube(Type*5+1 ,8)) obj.P(8,:)  = obj.Intersection(obj.m_PointE,obj.m_PointH); end
            if(Cube(Type*5+1 ,9)) obj.P(9,:)  = obj.Intersection(obj.m_PointA,obj.m_PointE); end
            if(Cube(Type*5+1 ,10)) obj.P(10,:)  = obj.Intersection(obj.m_PointB,obj.m_PointF); end
            if(Cube(Type*5+1 ,11)) obj.P(11,:)  = obj.Intersection(obj.m_PointC,obj.m_PointG); end
            if(Cube(Type*5+1 ,12)) obj.P(12,:)  = obj.Intersection(obj.m_PointD,obj.m_PointH); end
            
            obj.P(21,:)  = obj.m_PointA ;
            obj.P(22,:)  = obj.m_PointB ;
            obj.P(23,:)  = obj.m_PointC ;
            obj.P(24,:)  = obj.m_PointD ;
            
            obj.P(25,:)  = obj.m_PointE ;
            obj.P(26,:)  = obj.m_PointF ;
            obj.P(27,:)  = obj.m_PointG ;
            obj.P(28,:)  = obj.m_PointH ;
            
        end
        
        %         function Bbox = GetBoundingBox( obj, XYPolygon )
        %             Xmin = min(XYPolygon(:,1));
        %             Ymin = min(XYPolygon(:,2));
        %             Xmax = max(XYPolygon(:,1));
        %             Ymax = max(XYPolygon(:,2));
        %
        %             Bbox = [ Xmin Ymin; Xmax Ymax];
        %         end
        
        
        
        %Begin 4--------------------------------------------------------
        
        %         function [ApproxPolygon, ComputeSSATerm] = GetApproximatePolygon( obj )
        %             ApproxPolygon = [];
        %             ComputeSSATerm = [];
        %             switch( obj.GetCellType() )
        %                 case 0
        %                     %Outer Node
        %                 case 1
        %                     %Interior Node
        %                     ApproxPolygon = zeros( 4, 2);
        %                     ComputeSSATerm = zeros( 4, 1);
        %                     ApproxPolygon(1,:) = obj.m_PointA';
        %                     ApproxPolygon(2,:) = obj.InterpolateZeroLevelSet( obj.m_PointA, obj.m_PointB );
        %                     ApproxPolygon(3,:) = obj.InterpolateZeroLevelSet( obj.m_PointA, obj.m_PointD );
        %                     ApproxPolygon(4,:) = ApproxPolygon(1,:);
        %                     ComputeSSATerm(2) = 1;
        %             end
        %         end
        %
        % %Begin 5--------------------------------------------------------
        %         function Pts = InterpolateZeroLevelSet( obj, ptA, ptB )
        %             ax = ptA(1);
        %             ay = ptA(2);
        %             d = ptB - ptA;
        %             dx = d(1);
        %             dy = d(2);
        %             fun = @(t) obj.m_ImplicitFn( ax + t*dx , ay + t*dy );
        %             f1 = abs(fun(0));
        %             f2 = abs(fun(1));
        %
        %             %t0 = f1/(f1 + f2);
        %             t0 = 0.5;
        %             t = fsolve(fun,t0);
        %             %fMin = @(t) fun(t0).*fun(t0);
        %             %t = fminsearch(fMin,0);
        %             if( t >= 0 )
        %                 Pts = [ax + t*dx; ay + t*dy];
        %             else
        %                 Pts = ptA';
        %             end
        %         end
        %
        function Draw( obj, color )
            
            %              if(isempty(obj.m_ChildNodes) && IsInterior(obj))
            %                if(isempty(obj.m_ChildNodes) && IsOuter(obj))
            if(isempty(obj.m_ChildNodes) && IsBoundary(obj))
                
%                 obj.Draw_One_Node( color );
                  obj.Draw_Points( color );
                  obj.Draw_IntegrationPts_Points( color );
            end
            
            if(~isempty(obj.m_ChildNodes) )
                for i = 1: 8
                    obj.m_ChildNodes(i).Draw(color);
                end
            end
        end
        %         function Draw_One_Node(obj, color)
        %             diff = (obj.m_PointC - obj.m_PointA);
        %             w = abs(diff(1));
        %             h = abs(diff(2));
        %             hold on;
        %             %rectangle('Position',[obj.m_PointA(1),obj.m_PointB(2),w,h],'EdgeColor',color);
        %             p1 = obj.m_PointA;
        %             p2 = obj.m_PointB;
        %             p3 = obj.m_PointC;
        %             p4 = obj.m_PointD;
        %             p5 = p1;
        %
        %             p6 = obj.m_PointE;
        %             p7 = obj.m_PointF;
        %             p8 = obj.m_PointG;
        %             p9 = obj.m_PointH;
        %             p10 = p6;
        %
        %             A1=p1(1)>=0;A2=p2(1)>=0;A3=p3(1)>=0;A4=p4(1)>=0;A5=p5(1)>=0;A6=p6(1)>=0;A7=p7(1)>=0;A8=p8(1)>=0;
        %             B1=p1(2)>=0;B2=p2(2)>=0;B3=p3(2)>=0;B4=p4(2)>=0;B5=p5(2)>=0;B6=p6(2)>=0;B7=p7(2)>=0;B8=p8(2)>=0;
        %             C1=p1(3)>=0;C2=p2(3)>=0;C3=p3(3)>=0;C4=p4(3)>=0;C5=p5(3)>=0;C6=p6(3)>=0;C7=p7(3)>=0;C8=p8(3)>=0;
        %
        %             D1=A1&&A2&&A3&&A4&&A5&&A6&&A7&&A8;
        %             E1=B1&&B2&&B3&&B4&&B5&&B6&&B7&&B8;
        %             F1=C1&&C2&&C3&&C4&&C5&&C6&&C7&&C8;
        %
        %             if(D1&&E1&&F1)
        %               p = [p1';p2';p3';p4';p5';p6';p7';p8';p9';p10'];
        % %             line(p(:,1),p(:,2),p(:,3),'Color',color);
        %               plot3(p(:,1),p(:,2),p(:,3),'Color',color);
        %             end
        %              hold off;
        %         end
        %     end
        
        %         function Draw_One_Node(obj, color)
        %
        % % if(obj.m_PointA(1)+obj.m_PointG(1)<0) return;end
        % % if(obj.m_PointA(2)+obj.m_PointG(2)<0) return;end
        % % if(obj.m_PointA(3)+obj.m_PointG(3)<0) return;end
        %             Polygonise( obj);
        %
        %             hold on
        %
        %             p1 = obj.m_PointA;
        %             p2 = obj.m_PointB;
        %             p3 = obj.m_PointC;
        %             p4 = obj.m_PointD;
        %             p5 = p1;
        %
        %             p6 = obj.m_PointE;
        %             p7 = obj.m_PointF;
        %             p8 = obj.m_PointG;
        %             p9 = obj.m_PointH;
        %             p10 = p6;
        %
        % %               p = [p1';p2';p3';p4';p5';p6';p7';p8';p9';p10'];
        % %             line(p(:,1),p(:,2),p(:,3),'Color',color);
        % %               plot3(p(:,1),p(:,2),p(:,3),'Color',color);
        %
        %             for i=0:obj.p(1)-1
        %                 P1 = [obj.P(3*i+2,:)',obj.P(3*i+3,:)',obj.P(3*i+4,:)',obj.P(3*i+2,:)'];
        % %               line(p(:,1),p(:,2),p(:,3),'Color',color);
        %                 plot3(P1(1,:),P1(2,:),P1(3,:),'Color',color);
        % %                 plot3(P1(:,1),P1(:,2),P1(:,3),'Color',color);
        %             end
        %              hold off;
        %         end
        %     end
        
        function Draw_One_Node(obj, color)
            global Cube;
            if(obj.m_PointA(1)+obj.m_PointG(1)<0) return;end
            if(obj.m_PointA(2)+obj.m_PointG(2)<0) return;end
            if(obj.m_PointA(3)+obj.m_PointG(3)<0) return;end
            Polygonise( obj);
            
            hold on
            Type = obj.m_cellType;
            [numSolid, solid] = obj.getSolids();
            numSolid = Cube( Type*5+1 ,14);
            for i = 1 : numSolid-1
                %                  if(Cube( Type*5+1 ,13 )==0) i = numSolid;end % outer
                if(Cube( Type*5+1 ,13 )==1) i = numSolid;end % inner
                %                 i=1;
                numFace =   Cube( Type*5+i ,15);
                for j = 0 : numFace-1
                    numP =   Cube( Type*5+i ,16+j*7);
                    if(numP>6||numP<3) numP; end
                    switch( numP)
                        case 3
                            p1 =   Cube( Type*5+i ,16+j*7+1);
                            p2 =   Cube( Type*5+i ,16+j*7+2);
                            p3 =   Cube( Type*5+i ,16+j*7+3);
                            P3 = [obj.P(p1,:)',obj.P(p2,:)',obj.P(p3,:)',obj.P(p1,:)'];
                            plot3(P3(1,:),P3(2,:),P3(3,:),'Color',color);
                            
                        case 4
                            p1 =   Cube( Type*5+i ,16+j*7+1);
                            p2 =   Cube( Type*5+i ,16+j*7+2);
                            p3 =   Cube( Type*5+i ,16+j*7+3);
                            p4 =   Cube( Type*5+i ,16+j*7+4);
                            P4 = [obj.P(p1,:)',obj.P(p2,:)',obj.P(p3,:)',obj.P(p4,:)',obj.P(p1,:)'];
                            plot3(P4(1,:),P4(2,:),P4(3,:),'Color',color);
                        case 5
                            p1 =   Cube( Type*5+i ,16+j*7+1);
                            p2 =   Cube( Type*5+i ,16+j*7+2);
                            p3 =   Cube( Type*5+i ,16+j*7+3);
                            p4 =   Cube( Type*5+i ,16+j*7+4);
                            p5 =   Cube( Type*5+i ,16+j*7+5);
                            P5 = [obj.P(p1,:)',obj.P(p2,:)',obj.P(p3,:)',obj.P(p4,:)',obj.P(p5,:)',obj.P(p1,:)'];
                            plot3(P5(1,:),P5(2,:),P5(3,:),'Color',color);
                            
                        case 6
                            p1 =   Cube( Type*5+i ,16+j*7+1);
                            p2 =   Cube( Type*5+i ,16+j*7+2);
                            p3 =   Cube( Type*5+i ,16+j*7+3);
                            p4 =   Cube( Type*5+i ,16+j*7+4);
                            p5 =   Cube( Type*5+i ,16+j*7+5);
                            p6 =   Cube( Type*5+i ,16+j*7+6);
                            
                            P6 = [obj.P(p1,:)',obj.P(p2,:)',obj.P(p3,:)',obj.P(p4,:)',obj.P(p5,:)',obj.P(p6,:)',obj.P(p1,:)'];
                            
                            plot3(P6(1,:),P6(2,:),P6(3,:),'Color',color);
                    end
                end
            end
            
            hold off;
        end
        
        function Draw_Points(obj, color)
            if(obj.m_PointA(1)+obj.m_PointG(1)<0) return;end
            if(obj.m_PointA(2)+obj.m_PointG(2)<0) return;end
            if(obj.m_PointA(3)+obj.m_PointG(3)<0) return;end

            Polygonise( obj);
            [numSolid, solid] = obj.getSolids();
            hold on
            %             Type = obj.m_cellType;
            for i = 1 : numSolid
                [V,F,nV,nF,isTrim] = obj.getVerticesFaces(solid(i));

                for j = 0 : nF-1
%                     if(isTrim(F(3*j+1))==0) continue;end
%                     if(isTrim(F(3*j+2))==0) continue;end
%                     if(isTrim(F(3*j+3))==0) continue;end
                        
                    p1 =   V(F(3*j+1),:);    p2 =   V(F(3*j+2),:);     p3 =   V(F(3*j+3),:);
                    T = [p1',p2',p3',p1'];
                    plot3(T(1,:),T(2,:),T(3,:),'Color',color);
                end
            end
            hold off;
        end
    
    
        function Draw_IntegrationPts_Points(obj, color)

            if(obj.m_PointA(1)+obj.m_PointG(1)<0) return;end
            if(obj.m_PointA(2)+obj.m_PointG(2)<0) return;end
            if(obj.m_PointA(3)+obj.m_PointG(3)<0) return;end
 
%             Polygonise( obj);
            [numSolid, solid] = obj.getSolids();
            hold on
            nPt = 5;
            for i = 1 : numSolid
                [IntegrationPts] = AllocateIntegrationPoints_Random(obj, nPt,solid(i));
%                 IntegrationPts = IntegrationPts';
                scatter3(IntegrationPts(1,:)',IntegrationPts(2,:)',IntegrationPts(3,:)');
            end
            hold off;
        end
    end
    
    
 
    methods(Static)
        %         function DrawPolygon( PolygonCornerXY, color )
        %             line( PolygonCornerXY(:,1), PolygonCornerXY(:,2), 'Color', color, 'LineWidth',2 );
        %         end
        
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
