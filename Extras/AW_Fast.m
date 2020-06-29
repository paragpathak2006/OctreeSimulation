%=======================================================================%
%                     ______________  _____          __                 %
%                    / ____/ ____/  |/  / /   ____ _/ /_                %
%                   / /_  / /   / /|_/ / /   / __ `/ __ \               %
%                  / __/ / /___/ /  / / /___/ /_/ / /_/ /               %
%                 /_/    \____/_/  /_/_____/\__,_/_.___/                %
%                                                                       %
%                                                                       %
% Copyright (c) 2012, 2013                                              %
% Computation in Engineering, Technische Universitaet Muenchen          %
%                                                                       %
% This file is part of the MATLAB toolbox FCMLab. This library is free  %
% software; you can redistribute it and/or modify it under the terms of %
% the GNU General Public License as published by the Free Software      %
% Foundation; either version 3, or (at your option) any later version.  %
%                                                                       %
% This library is distributed in the hope that it will be useful,       %
% but WITHOUT ANY WARRANTY; without even the implied warranty of        %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          %
% GNU General Public License for more details.                          %
%                                                                       %
% You should have received a copy of the GNU General Public License     %
% along with this program; see the files COPYING respectively.          %
% If not, see <http://www.gnu.org/licenses/>.                           %
%                                                                       %
% In case of a scientific publication of results obtained using FCMLab, %
% we ask the authors to cite the introductory article                   %
%                                                                       %
% N. Zander, T. Bog, M. Elhaddad, R. Espinoza, H. Hu, A. F. Joly,       %
% C. Wu, P. Zerbe, A. Duester, S. Kollmannsberger, J. Parvizian,        %
% M. Ruess, D. Schillinger, E. Rank:                                    %
% "FCMLab: A Finite Cell Research Toolbox for MATLAB."                  %
% Submitted to Advances in Engineering Software, 2013					%
%                                                                       %
%=======================================================================%
 
classdef AW_Fast < handle
    
    methods
        %% constructor
        function obj = AW_Fast(integrationOrder, embedded_domain, nAdditionalEdges)
            obj.integrationOrder = integrationOrder;
            obj.embeddedDomain = embedded_domain;
            obj.mPt = nAdditionalEdges;
        end
        
        %% getCoordinates
        function [integrationPoints, polyX, polyY, nPolyPoints, ret] = getCoordinates(obj,geometry, support, boundary, flag,cellType, CornerPoints)
            ret = true;
            
            if( isempty(geometry) || isempty(boundary) || isempty(CornerPoints) )
                ret = false;
                integrationPoints = [];
                polyX = [];
                polyY = [];
                return;
            end
                 %Boundary cells -> generate Adaptive coordinates
                 if isa(geometry,'AbsArea')
                    nPt = obj.integrationOrder;
                    [integrationPoints, polyX, polyY, nPolyPoints, ret] = AW_Fast.GetIntegrationPoints( boundary, CornerPoints, cellType, nPt, obj.mPt);
                   
                    if( ret == false || isempty(polyX) || isempty(polyY) )
                        ret = false;
                        return;
                    end
                    
                 elseif isa(geometry,'AbsVolume')
                     %ToDo
                 else
                     Logger.ThrowException(['Unknown Geometric Type: ' class(geometry) '!!']);
                 end
        end
        
     
        
        %% getWeights
        function [integrationWeights, ret] = getWeights(obj,geometry, support, xCurve, yCurve, integrationPoints, flag, cellType, CornerPoints, nPolyPoints)
            ret = true;
            if( isempty(xCurve) || isempty(yCurve) || isempty(integrationPoints) || isempty(geometry) || isempty(CornerPoints))
                ret = false;
                integrationWeights = [];
                return;
            end
            fast = AW_Fast.Fast;
            if( isa(geometry,'AbsCurve') || flag )
            %Interior cells -> generate regular weights     
            integrationWeights1D = ...
                getGaussQuadratureWeights(obj.integrationOrder);
            
            if isa(geometry,'AbsCurve')
                integrationWeights = integrationWeights1D;
            elseif isa(geometry,'AbsArea')
                integrationWeights = ...
                    kron(integrationWeights1D,integrationWeights1D);
            elseif isa(geometry,'AbsVolume')
                integrationWeights = ...
                    kron(integrationWeights1D,...
                    kron(integrationWeights1D,integrationWeights1D));
            else
                Logger.ThrowException(['Unknown Geometric Type: ' class(geometry) '!!']);
            end
            else
             %Boundary cells -> generate Adaptive weights
             if isa(geometry,'AbsArea')
                 if( AW_Fast.Debug )
                     hold on;
                     %AW_Fast.DrawPolygon( CornerPoints, 'm', 1);
                     AW_Fast.DrawPoints( integrationPoints, 'k', 1);
                     hold off;
                 end
                 %[B,~] = obj.BMatrix( integrationPoints );
                  nPt = sqrt(size(integrationPoints,1));
                 polyDegree = 2*nPt - 1;
                 dDegree =  obj.embeddedDomain.getDomainDegree();
                 totalDegree = polyDegree + dDegree;
                 nPt = ceil( (totalDegree + 1)/2);
                 [W,s] = AW_Fast.Gauss_W_X_Scaled(0,1,nPt);
                 if( isempty( nPolyPoints ) )
                    if( ~fast )
                        [LineList, velocityEdgeNumbers, diagonalLength] = AW_Fast.GetApproxPolygonEdges( xCurve, yCurve, CornerPoints, cellType );
                    else
                        [LineList, velocityEdgeNumbers, diagonalLength] = AW_Fast.GetApproxPolygonEdges_XY( xCurve, yCurve, CornerPoints, cellType );
                    end
                 else
                     if( ~fast )
                        [LineList, velocityEdgeNumbers, diagonalLength] = AW_Fast.GetApproxPolygonEdges_From_PolyPoints( xCurve, yCurve, CornerPoints, cellType, nPolyPoints );
                     else
                        [LineList, velocityEdgeNumbers, diagonalLength] = AW_Fast.GetApproxPolygonEdges_From_PolyPoints_Fast( xCurve, yCurve, CornerPoints, cellType, nPolyPoints );   
                     end
                 end
                
                 Vel = [];
                 if( fast )
                    Vel = AW_Fast.ComputeVelocities_XY( LineList, velocityEdgeNumbers, s, xCurve, yCurve, diagonalLength );
                 else
                    AW_Fast.ComputeVelocities( LineList, velocityEdgeNumbers, s, xCurve, yCurve, diagonalLength );
                 end
                    
                 
                 %AW_Fast.GetIntegrationPoints_ApproxPolygon( boundary, CornerPoints, cellType, n)
                 
                 if( AW_Fast.Debug )
%                      length(LineList);
%                      hold on;
%                      if( fast )
%                         AW_Fast.DrawArrows_XY( LineList, 'r',1);
%                      else
%                         AW_Fast.DrawArrows( LineList, 'r',1);
%                      end
%                      hold off;
                 end
                 if( ~fast )
                    [B,L] = obj.GetBMatrix_And_LHS( integrationPoints, LineList, W, s, velocityEdgeNumbers );
                 else
                    [B,L] = obj.GetBMatrix_And_LHS_Fast( integrationPoints, LineList, W, s, velocityEdgeNumbers, Vel );
                 end
                 
                  %[L] = obj.GetLVector( LineList, velocityEdgeNumber );
                 integrationWeights = B*L;
             elseif isa(geometry,'AbsVolume')
                 %ToDo
             else
                 Logger.ThrowException(['Unknown Geometric Type: ' class(geometry) '!!']);
             end
            end
        end
       
          function [B,lhs] = GetBMatrix_And_LHS_Fast( obj, XY, LineList, W, s, velocityEdgeNumber, Vel )
            degree = obj.integrationOrder;
            nPoly = (degree+1)*(degree+2)/2;
            level = degree;

            [nPt,~] = size(XY);
            A = zeros(nPoly,nPt);
            lhs = zeros(nPoly,1);
            

            i = 1;
            for l = 0 : level
                for n = 0:l
                    m = l - n;

                    for k = 1 : nPt
                        x = XY(k,1); 
                        y = XY(k,2);
                        A(i,k) = x^m * y^n;
                    end
                    %A(i,:) = (XY(:,1).^m .* XY(:,2).^n)';
                    lhs(i) = obj.L_i_Fast(LineList, m , n, W, s, velocityEdgeNumber, Vel);
                    i = i + 1;
                end
            end

            B = pinv(A);
        end
        
        
        function [B,lhs] = GetBMatrix_And_LHS( obj, XY, LineList, W, s, velocityEdgeNumber )
            degree = obj.integrationOrder;
            nPoly = (degree+1)*(degree+2)/2;
            level = degree;

            [nPt,~] = size(XY);
            A = zeros(nPoly,nPt);
            lhs = zeros(nPoly,1);
            

            i = 1;
            for l = 0 : level
                for n = 0:l
                    m = l - n;

                    for k = 1 : nPt
                        x = XY(k,1); 
                        y = XY(k,2);
                        A(i,k) = x^m * y^n;
                    end
                    lhs(i) = obj.L_i(LineList, m , n, W, s, velocityEdgeNumber);
                    i = i + 1;
                end
            end

            B = pinv(A);
        end
        
        
        function [flag,cellType, CornerPoints] = IsInteriorCell( obj, geometry, support )
            [flag1,p1] = obj.getDomainIndexForSeedPoint([-1,-1],geometry, support);
            [flag2,p2] = obj.getDomainIndexForSeedPoint([1,-1],geometry, support);
            [flag3,p3] = obj.getDomainIndexForSeedPoint([1,1],geometry, support);
            [flag4,p4] = obj.getDomainIndexForSeedPoint([-1,1],geometry, support);
            flag = flag1 && flag2 && flag3 && flag4;
            s = sprintf('%d%d%d%d',flag1,flag2,flag3,flag4);
            cellType = bin2dec(s);
            CornerPoints =[p1;p2;p3;p4];
        end
        
       
       
        
       
        
        %To be Generalized to 3D
        function Pts = GetRandomPoints( obj, minX, minY, maxX, maxY, n)
            deltaX = (maxX - minX);
            deltaY = (maxY - minY);
            
            nRandomPoints = max(round(sqrt(n))*10,100000);
            randX = rand(nRandomPoints,1) * deltaX + minX;
            randY = rand(nRandomPoints,1) * deltaY + minY;
            randZ = zeros(nRandomPoints,1);
            
            
            randXYZ = [randX randY randZ];
            
            flags = obj.embeddedDomain.getDomainIndex_vectorized(randXYZ);
            
            Pts = randXYZ(flags,:);
            %Need to add a check here - size(Pts,1) < n then ???
            Pts = Pts(1:n,:);
        end
        
        function [domainIndex, globalCoordinates] = getDomainIndexForSeedPoint(obj,seedPoint,indexGeometry,support)
            globalCoordinates = AW_Fast.getGlobalCoords(seedPoint,indexGeometry,support );
            
            if( obj.embeddedDomain.getDomainIndex(globalCoordinates)== 2 )
                domainIndex = true;
            else
                domainIndex = false;
            end
        end
        
       
        
    end
    
    methods(Static)
         function globalCoordinates = getGlobalCoords( seedPoint,indexGeometry,support )
            localCoordinates = indexGeometry.mapLocalToGlobal(seedPoint);
            globalCoordinates = support.mapLocalToGlobal(localCoordinates);
         end
         
         
          function [Pt] = GetPtOfIntersection_XY( polyX, polyY, XY)
            Pt = [];
            %hold on;
            %AW_Fast.DrawPolyLine( [x' y'],'k',1);
            %hold off;
            [xi,yi] = polyxpoly( [XY(1) XY(3)] , [XY(2) XY(4)], polyX', polyY');
            if( ~isempty(xi) && ~isempty(yi) )
                [xr,yr] = AW_Fast.Nearest( polyX, polyY, xi, yi );
                Pt = [xr yr 0];
            end
           
        end
        
          function [Vert] = GetPtOfIntersection_Basic( polyX, polyY, Line2)
           
            Vert = [];
            v = Line2.getVertices();
            x = [v(1).getX() v(2).getX()];
            y = [v(1).getY() v(2).getY()];
            %hold on;
            %AW_Fast.DrawPolyLine( [x' y'],'k',1);
            %hold off;
            [xi,yi] = polyxpoly( x , y, polyX', polyY');
            if( ~isempty(xi) && ~isempty(yi) )
                [xr,yr] = AW_Fast.Nearest( polyX, polyY, xi, yi );
                Vert = Vertex( [xr yr 0] );
            end
           
        end
        
        function [Vert] = GetPtOfIntersection( boundaryLines, Line2, minX, minY, maxX, maxY )
           
            [polyX, polyY] = AW_Fast.GetTestLines( boundaryLines, minX, minY, maxX, maxY );
            v = Line2.getVertices();
            x = [v(1).getX() v(2).getX()];
            y = [v(1).getY() v(2).getY()];
            %AW_Fast.DrawPolyLine( [x' y'],'k',1);
            [xi,yi] = intersections( x , y, polyX, polyY, false);
            Vert = Vertex( xi(1), yi(1) );
        end
        
         function [polyX, polyY] = GetTestLines_Fast( XY_Boundary, minX, minY, maxX, maxY )
             XY_Boundary = XY_Boundary((minX <= XY_Boundary(:,1) & XY_Boundary(:,1) <= maxX & minY <= XY_Boundary(:,2) & XY_Boundary(:,2) <= maxY) |  (minX <= XY_Boundary(:,3) & XY_Boundary(:,3) <= maxX & minY <= XY_Boundary(:,4) & XY_Boundary(:,4) <= maxY),:);
             polyX = [XY_Boundary(:,1);XY_Boundary(:,3)];
             polyY = [XY_Boundary(:,2);XY_Boundary(:,4)];
             XY = unique([polyX polyY], 'rows');
             polyX = XY(:,1)';
             polyY = XY(:,2)';
         end
        
         function [polyX, polyY] = GetTestLines( boundaryLines, minX, minY, maxX, maxY )
            n = length(boundaryLines);
            
            polyX = [];
            polyY = [];
            for i = 1 : n
                vertices = boundaryLines(i).getVertices();
                v1 = vertices(1);
                v2 = vertices(2);
                v1X = v1.getX();
                v2X = v2.getX();
                v1Y = v1.getY();
                v2Y = v2.getY();
                if( (v1X >= minX && v1X <= maxX && v1Y >= minY && v1Y <= maxY) || (v2X >= minX && v2X <= maxX && v2Y >= minY && v2Y <= maxY) )
                    polyX = [polyX v1X v2X];
                    polyY = [polyY v1Y v2Y];
                end
            end
            % first remove duplicates and then sort based on angles
            XY = unique([polyX' polyY'], 'rows');
            polyX = XY(:,1);
            polyY = XY(:,2);
            %Matrix = [mod(atan2(polyY-0, polyX-0),2*pi) polyX polyY];
            %Matrix = sortrows(Matrix,1);
            %polyX = Matrix(:,2)';
            %polyY = Matrix(:,3)';
         end
         
          function [ I ] = Integrand_Fast_I2( line, s, m, n)
            %UNTITLED3 Summary of this function goes here
            %   Detailed explanation goes here
            vel = line.GetVelocity();
            vert = line.getVertices();
            x1 = vert(1).getX();
            x2 = vert(2).getX();
            y1 = vert(1).getY();
            y2 = vert(2).getY();
            X_s = s.*(x2-x1) + x1;
            Y_s = s.*(y2-y1) + y1;
            l = sqrt((x2-x1)^2 + (y2-y1)^2);
            I = ((X_s.^m .* Y_s.^n).*vel).*l;
          end
          function [ I ] = Integrand_Fast_I2_New( line, s, m, n, vel)
            %UNTITLED3 Summary of this function goes here
            %   Detailed explanation goes here
            x1 = line(1);
            x2 = line(3);
            y1 = line(2);
            y2 = line(4);
            X_s = s.*(x2-x1) + x1;
            Y_s = s.*(y2-y1) + y1;
            l = sqrt((x2-x1)^2 + (y2-y1)^2);
            I = ((X_s.^m .* Y_s.^n).*vel).*l;
          end
          
          
          
          
           function [ I ] = Integrand_Fast_I1_New( line, s, m, n )
            x1 = line(1);
            x2 = line(3);
            y1 = line(2);
            y2 = line(4);
            X_s = s.*(x2-x1) + x1;
            Y_s = s.*(y2-y1) + y1;

            I = (y2-y1).*(( X_s.^(m+1) .* Y_s.^n ) / (m+1));
          end
          
           function [ I ] = Integrand_Fast_I1( line, s, m, n )
            %UNTITLED3 Summary of this function goes here
            %   Detailed explanation goes here
            vert = line.getVertices();
            x1 = vert(1).getX();
            x2 = vert(2).getX();
            y1 = vert(1).getY();
            y2 = vert(2).getY();
            X_s = s.*(x2-x1) + x1;
            Y_s = s.*(y2-y1) + y1;

            I = (y2-y1).*(( X_s.^(m+1) .* Y_s.^n ) / (m+1));
          end
        
        
         
         function [ I ] = Integrand_Fast( line, s, m, n, option )
            %UNTITLED3 Summary of this function goes here
            %   Detailed explanation goes here
            if( option == 0 )
                vel = 0;
            else
                vel = line.GetVelocity();
            end
            phi = @(x,y) ( x.^(m+1) .* y.^n ) / (m+1);
            f = @(x,y) x.^m .* y.^n;
            vert = line.getVertices();
            x1 = vert(1).getX();
            x2 = vert(2).getX();
            y1 = vert(1).getY();
            y2 = vert(2).getY();
            X_s = s.*(x2-x1) + x1;
            Y_s = s.*(y2-y1) + y1;

            
            l = sqrt((x2-x1)^2 + (y2-y1)^2);
            I1 = -(y2-y1).*(phi(X_s, Y_s));
            if( option == 1)
                I2 = (f(X_s,Y_s).*vel).*l;
            else
                I2 = 0;
            end
                
            I = I1 + I2;
          end
        
          function [ I ] = Integrand( line, s, XCurve, yCurve, diagonalLength, m, n, option )
            %UNTITLED3 Summary of this function goes here
            %   Detailed explanation goes here
            if( option == 0 )
                vel = @(s) 0;
            else
                vel = @(s) AW_Fast.Velocity(line,s, XCurve, yCurve, diagonalLength);
            end
            phi = @(x,y) ( x.^(m+1) .* y.^n ) / (m+1);
            f = @(x,y) x.^m .* y.^n;
            vert = line.getVertices();
            x1 = vert(1).getX();
            x2 = vert(2).getX();
            y1 = vert(1).getY();
            y2 = vert(2).getY();

            X_s = s.*(x2-x1) + x1;
            Y_s = s.*(y2-y1) + y1;

            
            l = sqrt((x2-x1)^2 + (y2-y1)^2);
            I1 = -(y2-y1).*(phi(X_s, Y_s));
            I2 = (f(X_s,Y_s).*vel(s)).*l;
            I = I1 + I2;
          end
          function [d] = GetDistanceToBoundary_Circle( x0, y0, n, r )
              nx = n(1);
              ny = n(2);
              
              a = nx.^2 + ny.^2;
              b = 2*(nx.*x0 + ny.*y0);
              c = x0.^2 + y0.^2 - r.^2;
              
              Disc = sqrt(b.^2 - 4.*a.*c);
              t1 = (-b + Disc)./(2*a);
              t2 = (-b - Disc)./(2*a);
              
              d = min(abs(t1),abs(t2));
          end
          
             function [ v_n ] = Velocity_Fast( line1, s, XCurve, YCurve, diagonal_length )
            %UNTITLED6 Summary of this function goes here
            %   Detailed explanation goes here
            
            [m,~] = size(s);
            v_n = zeros( m , 1 );
            %Constant Normal
            x1 = line1(1);
            y1 = line1(2);
            x2 = line1(3);
            y2 = line1(4);
            n = [-(y2-y1) (x2-x1)];
            n = -n / norm(n);
            for i = 1 : m
                x_0 = x1 + s(i)*(x2-x1);
                y_0 = y1 + s(i)*(y2-y1);
                X = [x_0 y_0];
                pt1 = X - diagonal_length.*n;
                pt2 = X + diagonal_length.*n;
                line2 = [pt1 pt2];
%                 if( AW_Fast.Debug )
%                     hold on;
%                     AW_Fast.DrawLines( line2, 'y', 1);
%                     hold off;
                end

                   d1 = AW_Fast.GetDistanceToBoundary_Circle(x_0, y_0, n, 1 );
                   v_n(i) = d1;
%                   d2 = AW_Fast.GetDistanceToBoundary_Circle(x_0, y_0, n, 0.25 );
%                   v_n(i) = min(d1,d2);

%                 v3 = AW_Fast.GetPtOfIntersection_XY( XCurve, YCurve, line2 ); 
%                 if( isempty(v3) )
%                     v_n(i) = 0;
%                     %Bug here
%                 else
%                     x3 = v3(1);
%                     y3 = v3(2);
%                     vx = [x_0 x3]';
%                     vy = [y_0 y3]';
%                     if( AW_Fast.Debug )
%                         hold on;
%                         %AW_Fast.DrawPolyLine( [vx vy], 'g', 1); 
%                         hold off;
%                     end
% 
%                     v = [x3-x_0;y3-y_0];
%                     v_n(i) = v'*n';
%                 end
       end
       
           function [ v_n ] = Velocity( line1, s, XCurve, YCurve, diagonal_length )
            %UNTITLED6 Summary of this function goes here
            %   Detailed explanation goes here
            
            [m,~] = size(s);
            v_n = zeros( m , 1 );
            %Constant Normal
            n = line1.calcNormalVector(s(1));
            for i = 1 : m
                X = line1.GetXYZ_length_scaling( s(i) );
                x_0 = X(1);
                y_0 = X(2);
                pt1 = X - diagonal_length.*[n 0];
                pt2 = X + diagonal_length.*[n 0];
                v1 = Vertex(pt1);
                v2 = Vertex( pt2 );
                line2 = Line( v1, v2 );
%                 if( AW_Fast.Debug )
%                     hold on;
%                     AW_Fast.DrawLines( line2, 'y', 1);
%                     hold off;
%                 end
                
                v3 = AW_Fast.GetPtOfIntersection_Basic( XCurve, YCurve, line2 ); 
                if( isempty(v3) )
                    v_n(i) = 0;
                    %Bug here
                else
                    coord = v3.getCoords();
                    x1 = coord(1);
                    y1 = coord(2);
                    vx = [x_0 x1]';
                    vy = [y_0 y1]';
                    if( AW_Fast.Debug )
                        hold on;
                        AW_Fast.DrawPolyLine( [vx vy], 'g', 1); 
                        hold off;
                    end

                    v = [x1-x_0;y1-y_0];
                    v_n(i) = v'*n';
                end
            end
       end
       
       function [ xMin,yMin ] = Nearest( x_0, y_0, X,Y )
            %UNTITLED3 Summary of this function goes here
            %   Detailed explanation goes here

            [m,~] = size(X);
            [m1,~] = size(Y);

            if (m == 0 || m1 ~= m)
                xMin = 0;
                yMin = 0;
                return;
            end

            d = @(x,y) sqrt((x-x_0).^2 + (y-y_0).^2);
            xMin = X(1);
            yMin = Y(1);
            dMin = d(xMin,yMin);
            for i = 2 : m
                xC = X(i);
                yC = Y(i);
                dist = d(xC,yC);
                if( dist < dMin )
                    xMin = xC;
                    yMin = yC;
                    dMin = dist;
                end
            end

       end
        
          function [ W , X ] = Gauss_W_X_Scaled( a,b, nPt )
        %UNTITLED5 Summary of this function goes here
        %   Detailed explanation goes here

            W = getGaussQuadratureWeights( nPt );
            X = getGaussQuadratureCoordinates( nPt );

            c1 = (b-a)/2;
            c2 = (b+a)/2;
            W = W*c1;
            X = X*c1 + c2;
%             for i = 1 : nPt
%                 W(i) = W(i)*c1;
%                 X(i) = X(i)*c1 + c2;
%             end
          end
          
            function [lhs_i] = L_i_Fast(LineList, m , n,W, s, velocityEdgeNumbers, vel)
            k = size(LineList,1);
            lhs_i = 0;
            for i = 1 : k
                  %Zero Velocity
                    %lhs_i = lhs_i + W0'*AW_Fast.Integrand( LineList(i), X0, xCurve, yCurve, diagonalLength, m , n, 0 );
                    lhs_i = lhs_i + W'*AW_Fast.Integrand_Fast_I1_New( LineList(i,:), s, m , n );
            end
            k = length(velocityEdgeNumbers);
            for i = 1 : k
                  %Compute Velocity
                    %lhs_i = lhs_i + W0'*AW_Fast.Integrand( LineList(i), X0, xCurve, yCurve, diagonalLength, m , n, 0 );
                    lhs_i = lhs_i + W'*AW_Fast.Integrand_Fast_I2_New( LineList(velocityEdgeNumbers(i),:), s, m , n, vel(:,i) );
            end
            
         end
         
        
         function [lhs_i] = L_i(LineList, m , n,W, s, velocityEdgeNumbers)
            k = length(LineList);
            lhs_i = 0;
            for i = 1 : k
                  %Zero Velocity
                    %lhs_i = lhs_i + W0'*AW_Fast.Integrand( LineList(i), X0, xCurve, yCurve, diagonalLength, m , n, 0 );
                    lhs_i = lhs_i + W'*AW_Fast.Integrand_Fast_I1( LineList(i), s, m , n );
            end
            k = length(velocityEdgeNumbers);
            for i = 1 : k
                  %Compute Velocity
                    %lhs_i = lhs_i + W0'*AW_Fast.Integrand( LineList(i), X0, xCurve, yCurve, diagonalLength, m , n, 0 );
                    lhs_i = lhs_i + W'*AW_Fast.Integrand_Fast_I2( LineList(velocityEdgeNumbers(i)), s, m , n );
            end
            
         end
         
         
             function [integrationPoints, polyX, polyY, nPolyPoints, ret] = GetIntegrationPoints( boundaryLines, CornerPoints, cellType, n, mPt)
             ret = true;
             p1 = CornerPoints(1,:);
             p2 = CornerPoints(2,:);
             p3 = CornerPoints(3,:);
             p4 = CornerPoints(4,:);
%              Vertex1 = Vertex( p1 );
%              Vertex2 = Vertex( p2 );
%              Vertex3 = Vertex( p3 );
%              Vertex4 = Vertex( p4 );
%              Line1 = Line( Vertex1, Vertex2 );
%              Line2 = Line( Vertex2, Vertex3 );
%              Line3 = Line( Vertex3, Vertex4 );
%              Line4 = Line( Vertex4, Vertex1 );
             Line1 = [p1(1:2) p2(1:2)];
             Line2 = [p2(1:2) p3(1:2)];
             Line3 = [p3(1:2) p4(1:2)];
             Line4 = [p4(1:2) p1(1:2)];
             minX = min(CornerPoints(:,1));
             minY = min(CornerPoints(:,2));
             maxX = max(CornerPoints(:,1));
             maxY = max(CornerPoints(:,2));
             diagonalLength = sqrt( (minX - maxX).^2 + (maxY - minY).^2 );
             nPolyPoints = [];
             
             [polyX, polyY] = AW_Fast.GetTestLines_Fast( boundaryLines, minX, minY, maxX, maxY );
             if( isempty(polyX) || isempty(polyY) )
                 ret = false;
                 integrationPoints = [];
                 return;
             end
             if( AW_Fast.Debug )
                hold on;
                AW_Fast.DrawPolyLine( [polyX' polyY'], 'r', 1);
                %AW_Fast.DrawPoints( [polyX' polyY'], 'b', 0.5);
                AW_Fast.DrawPolygon( CornerPoints, 'm', 1);
                hold off;
            end
             
              
             
              switch cellType
                case 0
                    %Exterior Cell
                    integrationPoints = [];
                case 1
                    directionVector = (p2 - p1);
                    directionVector = directionVector / norm(directionVector);
                    [integrationPoints, nPolyPoints] = AW_Fast.GetPolarCoordinates( p1, directionVector,  polyX, polyY, diagonalLength, n, mPt );
                case 2 
                    directionVector = (p3 - p2);
                    directionVector = directionVector / norm(directionVector);
                    [integrationPoints, nPolyPoints] = AW_Fast.GetPolarCoordinates( p2, directionVector, polyX, polyY, diagonalLength,  n, mPt );
                case 3 
                    directionVector =  (p4 - p1);
                    directionVector = directionVector / norm(directionVector);
                    [integrationPoints, nPolyPoints] = AW_Fast.GetCartesianCoordinates_XY( Line1, Line3, directionVector, polyX, polyY, diagonalLength, n, mPt );
                case 4 
                    directionVector = (p4 - p3);
                    directionVector = directionVector / norm(directionVector);
                    [integrationPoints, nPolyPoints] = AW_Fast.GetPolarCoordinates( p3, directionVector, polyX, polyY, diagonalLength,  n, mPt );
                case 5
                    %Ambiguous case
                    integrationPoints = [];
                case 6 
                    directionVector =  (p1 - p2);
                    directionVector = directionVector / norm(directionVector);
                    [integrationPoints, nPolyPoints] = AW_Fast.GetCartesianCoordinates_XY( Line2, Line4, directionVector,  polyX, polyY, diagonalLength,  n, mPt );
                case 7
                    directionVector =  (p1 - p2);
                    directionVector = directionVector / norm(directionVector);
                    [integrationPoints, nPolyPoints] = AW_Fast.GetCartesianCoordinates_XY( Line2, Line4, directionVector,  polyX, polyY, diagonalLength, n, mPt );
                case 8
                    directionVector = (p1 - p4);
                    directionVector = directionVector / norm(directionVector);
                    [integrationPoints, nPolyPoints] = AW_Fast.GetPolarCoordinates( p4, directionVector, polyX, polyY, diagonalLength, n, mPt );
                case 9
                    directionVector =  (p2 - p1);
                    directionVector = directionVector / norm(directionVector);
                    [integrationPoints, nPolyPoints] = AW_Fast.GetCartesianCoordinates_XY( Line4, Line2, directionVector,  polyX, polyY, diagonalLength, n, mPt );
                case 10
                    %Ambiguous case
                    integrationPoints = [];
                case 11
                    directionVector =  (p2 - p1);
                    directionVector = directionVector / norm(directionVector);
                    [integrationPoints, nPolyPoints] = AW_Fast.GetCartesianCoordinates_XY( Line4,Line2, directionVector,  polyX, polyY, diagonalLength, n, mPt );
                case 12
                    directionVector =  (p1 - p4);
                    directionVector = directionVector / norm(directionVector);
                    [integrationPoints, nPolyPoints] = AW_Fast.GetCartesianCoordinates_XY( Line3,Line1, directionVector,  polyX, polyY, diagonalLength, n, mPt );
                case 13
                    directionVector =  (p2 - p1);
                    directionVector = directionVector / norm(directionVector);
                    [integrationPoints, nPolyPoints] = AW_Fast.GetCartesianCoordinates_XY( Line4,Line2, directionVector,  polyX, polyY, diagonalLength, n, mPt);
                case 14
                    directionVector =  (p1 - p2);
                    directionVector = directionVector / norm(directionVector);
                    [integrationPoints, nPolyPoints] = AW_Fast.GetCartesianCoordinates_XY( Line2,Line4, directionVector,  polyX, polyY, diagonalLength, n, mPt );
                case 15
                    %Interior Cell
                    integrationPoints = [];
              end
             
    end
                 
 
         function [integrationPts, nPolyPoints] = GetPolarCoordinates( basePoint, directionVector, polyX, polyY, diagonalLength, nPt, mPt )
             integrationPts = zeros( nPt*nPt, 3 );
              
             deltaTheta = deg2rad(90) / (nPt + 1);
              
             k = 1;
             
             baseTheta = mod(atan2( directionVector(2), directionVector(1)), 2*pi );
             theta = baseTheta + deltaTheta;
             base = zeros(nPt,3);
             base(:,1) = basePoint(1);
             base(:,2) = basePoint(2);
             base(:,3) = basePoint(3);
             
             
             nPolyPoints = zeros( nPt, 3);
             for i = 1 : nPt
                 n = [cos(theta) sin(theta) 0];
%                  v2 = Vertex( basePoint + diagonalLength.*n );
%                  perpLine = Line( v1 , v2 );
%                  v3 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, perpLine );
                 Pt2 = basePoint + diagonalLength.*n;
                 Pt3 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, [basePoint(1:2) Pt2(1:2)] );
                 %length = norm((v3.getCoords() - v1.getCoords() ));
                 length = norm(Pt3 - basePoint);
                 dl = length / (nPt + 1) ; 
                 delta = 0+dl:dl:length-dl;
                 integrationPts( k : k+nPt-1, :) = base + delta'*n;
                 nPolyPoints(i,:) = basePoint + length*n;
                 k = k + nPt;
                 theta = theta + deltaTheta;
             end
             
             if( ~isempty(mPt) )
                  deltaTheta = deg2rad(90) / (mPt + 1);
                  nPolyPoints = zeros( mPt, 3);
                  theta = baseTheta + deltaTheta;
                  for i = 1 : mPt
                     n = [cos(theta) sin(theta) 0];
%                      v2 = Vertex( basePoint + diagonalLength.*n );
%                      perpLine = Line( v1 , v2 );
%                      v3 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, perpLine );
                     Pt2 = basePoint + diagonalLength.*n;
                     Pt3 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, [basePoint(1:2) Pt2(1:2)] );
                     %length = norm((v3.getCoords() - v1.getCoords() ));
                     length = norm(Pt3 - basePoint);
                     nPolyPoints(i,:) = basePoint + length*n;
                     theta = theta + deltaTheta;
                 end
             end
         end
        
         
          function [integrationPts, nPolyPoints] = GetCartesianCoordinates_XY( baseLineXY, oppLineXY, directionVector, polyX, polyY, diagonalLength, nPt, mPt )
             integrationPts = zeros( nPt*nPt, 3 );
             nPolyPoints = zeros( nPt, 3 );
             p1 = [baseLineXY(1) baseLineXY(2) 0];
             p2 = [baseLineXY(3) baseLineXY(4) 0];
             
             %length = norm((p2 - p1));
             %dl = length / (nPt+1);
             
             %The bug is here
             if( abs(p2(2) - p1(2)) < 1e-12 )
                min_val = p1(1);
                max_val = p2(1);
                index = 1;
             else
                min_val = p1(2);
                max_val = p2(2);
                index = 2;
             end
             dl = (max_val - min_val)/(nPt+1);
             delta = min_val + dl : dl : max_val - dl; 
             xPts = zeros( nPt, 3);
             xPts(:,1) = p1(1);
             xPts(:,2) = p1(2);
             xPts(:,3) = 0; 
             xPts(:,index) = delta';
           
             lineX = [oppLineXY(1) oppLineXY(3)];
             lineY = [oppLineXY(2) oppLineXY(4)];
          
             k = 1;
             j = 1;
             
             for i = 1 : nPt
                 flag = true;
                 perpLineXY = [xPts(i,1:2) xPts(i,1:2) + diagonalLength.*directionVector(1:2)];
                 p3 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, perpLineXY );
                 if( isempty(p3) )
                     p3 = AW_Fast.GetPtOfIntersection_XY( lineX, lineY, perpLineXY);
                     flag = false;
                 end
                 if( flag )
                    nPolyPoints(j,:) = p3;
                    j = j + 1;
                 end
                 
                 %thisLine = Line(v1,v3);
                 %AW_Fast.DrawLines( thisLine, 'g', 1);
                 length = norm(p3 - xPts(i,:));
                 dl = length /  (nPt + 1);
                 deltaY = dl : dl : length - dl;
                 base = zeros(nPt,3);
                 base(:,1) = xPts(i,1);
                 base(:,2) = xPts(i,2);
                 base(:,3) = xPts(i,3);
                 integrationPts( k : k+nPt-1, :) = base + deltaY'*directionVector;
                 k = k + nPt;
             end
             if( ~isempty(mPt) )
                 j = 1;
                 dl = (max_val - min_val)/(mPt+1);
                 delta = min_val + dl : dl : max_val - dl; 
                 xPts = zeros( mPt, 3);
                 xPts(:,1) = p1(1);
                 xPts(:,2) = p1(2);
                 xPts(:,3) = 0; 
                 xPts(:,index) = delta';
                  for i = 1 : mPt
                     perpLineXY = [xPts(i,1:2) xPts(i,1:2) + diagonalLength.*directionVector(1:2)];
                     p3 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, perpLineXY );
                     if( isempty(p3) )
                         p3 = AW_Fast.GetPtOfIntersection_XY( lineX, lineY, perpLineXY);
                     end
                     nPolyPoints(j,:) = p3;
                     j = j + 1;
                  end
             end
             nPolyPoints = nPolyPoints(1:j-1,:);
             nPolyPoints = flipud(nPolyPoints);
         end
         
         
         function [integrationPts, nPolyPoints] = GetCartesianCoordinates( baseLine, oppLine, directionVector, polyX, polyY, diagonalLength, nPt, mPt )
             integrationPts = zeros( nPt*nPt, 3 );
             nPolyPoints = zeros( nPt, 3 );
             vertices = baseLine.getVertices();
             p1 = vertices(1).getCoords();
             p2 = vertices(2).getCoords();
             
             %length = norm((p2 - p1));
             %dl = length / (nPt+1);
             
             %The bug is here
             if( abs(p2(2) - p1(2)) < 1e-12 )
                min_val = p1(1);
                max_val = p2(1);
                index = 1;
             else
                min_val = p1(2);
                max_val = p2(2);
                index = 2;
             end
             dl = (max_val - min_val)/(nPt+1);
             delta = min_val + dl : dl : max_val - dl; 
             xPts = zeros( nPt, 3);
             xPts(:,1) = p1(1);
             xPts(:,2) = p1(2);
             xPts(:,3) = 0; 
             xPts(:,index) = delta';
             

           oppVertices = oppLine.getVertices();
           oppP1 = oppVertices(1).getCoords();
           oppP2 = oppVertices(2).getCoords();
           lineX = [oppP1(1) oppP2(1)];
           lineY = [oppP1(2) oppP2(2)];
          
             k = 1;
             j = 1;
             
             for i = 1 : nPt
                 flag = true;
                 v1 = Vertex( xPts(i,:) );
                 v2 = Vertex( xPts(i,:) + diagonalLength.*directionVector);
                 perpLine = Line( v1 , v2 );
                 v3 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, perpLine );
                 if( isempty(v3) )
                     v3 = AW_Fast.GetPtOfIntersection_Basic( lineX, lineY, perpLine);
                     flag = false;
                 end
                 if( flag )
                    nPolyPoints(j,:) = v3.getCoords();
                    j = j + 1;
                 end
                 
                 %thisLine = Line(v1,v3);
                 %AW_Fast.DrawLines( thisLine, 'g', 1);
                 length = norm((v3.getCoords() - v1.getCoords() ));
                 dl = length /  (nPt + 1);
                 deltaY = dl : dl : length - dl;
                 base = zeros(nPt,3);
                 base(:,1) = xPts(i,1);
                 base(:,2) = xPts(i,2);
                 base(:,3) = xPts(i,3);
                 integrationPts( k : k+nPt-1, :) = base + deltaY'*directionVector;
                 k = k + nPt;
             end
             if( ~isempty(mPt) )
                 j = 1;
                 dl = (max_val - min_val)/(mPt+1);
                 delta = min_val + dl : dl : max_val - dl; 
                 xPts = zeros( mPt, 3);
                 xPts(:,1) = p1(1);
                 xPts(:,2) = p1(2);
                 xPts(:,3) = 0; 
                 xPts(:,index) = delta';
                  for i = 1 : mPt
                     flag = true;
                     v1 = Vertex( xPts(i,:) );
                     v2 = Vertex( xPts(i,:) + diagonalLength.*directionVector);
                     perpLine = Line( v1 , v2 );
                     v3 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, perpLine );
                     if( isempty(v3) )
                         v3 = AW_Fast.GetPtOfIntersection_Basic( lineX, lineY, perpLine);
                         flag = false;
                     end
                     if( flag )
                        nPolyPoints(j,:) = v3.getCoords();
                        j = j + 1;
                     end
                  end
             end
             nPolyPoints = nPolyPoints(1:j-1,:);
             nPolyPoints = flipud(nPolyPoints);
         end
         
         function [LineList, velocityEdgeNumber, diagonalLength] = GetApproxPolygonEdges_Fast( CornerPoints, cellType )
            Vertex1 = Vertex( CornerPoints(1,:) );
            Vertex2 = Vertex( CornerPoints(2,:) );
            Vertex3 = Vertex( CornerPoints(3,:) );
            Vertex4 = Vertex( CornerPoints(4,:) );
            Line1 = Line( Vertex1, Vertex2 );
            Line2 = Line( Vertex2, Vertex3 );
            Line3 = Line( Vertex3, Vertex4 );
            Line4 = Line( Vertex4, Vertex1 );
            minX = min(CornerPoints(:,1));
            minY = min(CornerPoints(:,2));
            maxX = max(CornerPoints(:,1));
            maxY = max(CornerPoints(:,2));
            diagonalLength = sqrt( (minX - maxX).^2 + (maxY - minY).^2 );
            LineList = [Line1; Line2; Line3; Line4];
            
            %if( AW_Fast.Debug )
                %AW_Fast.DrawPolyLine( [polyX' polyY'], 'c', 1);
                %AW_Fast.DrawPoints( [polyX' polyY'], 'b', 0.5);
            %end
            
            switch cellType
                case 0
                    %Exterior Cell
                    LineList = [];
                case 1
                    velocityEdgeNumber = 3;
                case 2 
                    velocityEdgeNumber = 4;
                case 3 
                    velocityEdgeNumber = 3;
                case 4 
                    velocityEdgeNumber = 1;
                case 5
                    %Ambiguous case
                    LineList = [];
                case 6 
                    velocityEdgeNumber = 4;
                case 7
                    velocityEdgeNumber = 3;
                case 8
                   velocityEdgeNumber = 1;
                case 9
                  velocityEdgeNumber = 2;
                case 10
                    %Ambiguous case
                    LineList = [];
                case 11
                    velocityEdgeNumber = 3;
                case 12
                    velocityEdgeNumber = 1;
                case 13
                    velocityEdgeNumber = 1;
                case 14
                    velocityEdgeNumber = 1;
                case 15
                    %Interior Cell
                    LineList = [];
            
            end
         end
        
         function [BoundaryLineList] = GetBoundaryLineList_XY(nPolyPoints)
%             [n,~] = size(nPolyPoints);
%             BoundaryLineList = zeros(n-1,4);
%             
%             for i = 1 : n-1
%                 BoundaryLineList(i,:) = [nPolyPoints(i,1:2) nPolyPoints(i+1,1:2)];
%             end
            BoundaryLineList = [nPolyPoints(1:end-1,1:2) nPolyPoints(2:end,1:2)];
         end
         
         function [BoundaryLineList] = GetBoundaryLineList(nPolyPoints)
            [n,~] = size(nPolyPoints);
            
            BoundaryLineList = Line.empty(n-1,0);
            
            for i = 1 : n-1
                v1 = Vertex( nPolyPoints(i,:) );
                v2 = Vertex( nPolyPoints(i+1,:) );
                BoundaryLineList(i) = Line(v1,v2);
            end
         end
         
         function [LineList, velocityEdgeNumbers, diagonalLength] = GetApproxPolygonEdges_From_PolyPoints_Fast( polyX, polyY, CornerPoints, cellType, nPolyPoints )
            Vertex1 = CornerPoints(1,1:2);
            Vertex2 = CornerPoints(2,1:2);
            Vertex3 = CornerPoints(3,1:2);
            Vertex4 = CornerPoints(4,1:2);
            Line1 = [Vertex1 Vertex2];
            Line2 = [Vertex2 Vertex3];
            Line3 = [Vertex3 Vertex4];
            Line4 = [Vertex4 Vertex1];
            minX = min(CornerPoints(:,1));
            minY = min(CornerPoints(:,2));
            maxX = max(CornerPoints(:,1));
            maxY = max(CornerPoints(:,2));
            diagonalLength = sqrt( (minX - maxX).^2 + (maxY - minY).^2 );
            
          
            
            %if( AW_Fast.Debug )
                %AW_Fast.DrawPolyLine( [polyX' polyY'], 'c', 1);
                %AW_Fast.DrawPoints( [polyX' polyY'], 'b', 0.5);
            %end
            
            switch cellType
                case 0
                    %Exterior Cell
                    LineList = [];
                case 1
                    v1 = Vertex1;
                    v2 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line1);
                    v3 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line4);
                    nPolyPoints = [v2;nPolyPoints;v3]; 
                    Line2 = AW_Fast.GetBoundaryLineList_XY(nPolyPoints);
                    Line3 = [v3(1:2) v1(1:2)];
                    LineList = [Line1;Line2;Line3];
                    velocityEdgeNumbers = 2:size(LineList,1)-1;
                case 2 
                    v1 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY,  Line1 );
                    v2 = Vertex2;
                    v3 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line2);
                    Line1 = [ v1(1:2)  v2(1:2)];
                    Line2 = [v2(1:2)  v3(1:2)];
                    nPolyPoints = [v3;nPolyPoints;v1]; 
                    Line3 = AW_Fast.GetBoundaryLineList_XY(nPolyPoints);
                    LineList = [Line1;Line2;Line3];
                    velocityEdgeNumbers = 3:size(LineList,1);
                case 3 
                    v1 = Vertex1;
                    v2 = Vertex2;
                    v3 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line2 );
                    v4 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line4 );
                    Line1 = [ v1(1:2)  v2(1:2)];
                    Line2 = [v2(1:2)  v3(1:2)];
                    nPolyPoints = [v3;nPolyPoints;v4]; 
                    Line3 = AW_Fast.GetBoundaryLineList_XY(nPolyPoints);
                    Line4 = [v4(1:2)  v1(1:2)];
                    LineList = [Line1;Line2; Line3; Line4];
                    velocityEdgeNumbers = 3:size(LineList,1)-1;
                case 4 
                    v1 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line2);
                    v2 = Vertex3;
                    v3 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line3);
                    Line1 = [ v1(1:2)  v2(1:2)];
                    Line2 = [v2(1:2)  v3(1:2)];
                    nPolyPoints = [v3;nPolyPoints;v1]; 
                    Line3 = AW_Fast.GetBoundaryLineList_XY(nPolyPoints);
                    LineList = [Line1; Line2; Line3];
                    velocityEdgeNumbers = 3:size(LineList,1);
                case 5
                    %Ambiguous case
                    LineList = [];
                case 6 
                    v1 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line1);
                    v2 = Vertex2;
                    v3 = Vertex3;
                    v4 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line3 );
                    Line1 = [ v1(1:2) v2(1:2)];
                    Line2 = [v2(1:2)  v3(1:2)];
                    Line3 = [v3(1:2)  v4(1:2)];
                    nPolyPoints = [v4;nPolyPoints;v1]; 
                    Line4 = AW_Fast.GetBoundaryLineList_XY(nPolyPoints);
                    LineList = [Line1; Line2; Line3; Line4];
                    velocityEdgeNumbers = 4:size(LineList,1);
                case 7
                    v1 = Vertex1;
                    v2 = Vertex2;
                    v3 = Vertex3;
                    v4 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line3);
                    v5 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line4);
                    if( isempty(v4) )
                        v4 = v3;
                    end
                    if( isempty(v5) )
                        v5 = v1;
                    end
                    Line1 = [ v1(1:2)  v2(1:2)];
                    Line2 = [v2(1:2)  v3(1:2)];
                    Line3 = [v3(1:2)  v4(1:2)];
                    nPolyPoints = [v4;nPolyPoints;v5]; 
                    Line4 = AW_Fast.GetBoundaryLineList_XY(nPolyPoints);
                    Line5 = [v5(1:2)  v1(1:2)];
                    LineList = [Line1;Line2;Line3;Line4;Line5];
                    velocityEdgeNumbers = 4:size(LineList,1)-1;
                case 8
                    v1 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line3);
                    v2 = Vertex4;
                    v3 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line4 );
                    Line1 = [ v1(1:2)  v2(1:2)];
                    Line2 = [v2(1:2)  v3(1:2)];
                    nPolyPoints = [v3;nPolyPoints;v1]; 
                    Line3 = AW_Fast.GetBoundaryLineList_XY(nPolyPoints);
                    LineList = [Line1; Line2; Line3];
                    velocityEdgeNumbers = 3:size(LineList,1);
                case 9
                    v1 = Vertex1;
                    v2 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line1);
                    v3 =AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line3);
                    v4 = Vertex4;
                    Line1 = [ v1(1:2)  v2(1:2)];
                    nPolyPoints = [v2;nPolyPoints;v3]; 
                    Line2 = AW_Fast.GetBoundaryLineList_XY(nPolyPoints);
                    Line3 = [ v3(1:2)  v4(1:2)];
                    Line4 = [ v4(1:2)  v1(1:2)];
                    LineList = [Line1; Line2; Line3; Line4];
                   velocityEdgeNumbers = 2:size(LineList,1)-2;
                case 10
                    %Ambiguous case
                    LineList = [];
                case 11
                    v1 = Vertex1;
                    v2 = Vertex2;
                    v3 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line2);
                    v4 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line3 );
                    v5 = Vertex4;
                     if( isempty(v3) )
                        v3 = v2;
                    end
                    if( isempty(v4) )
                        v4 = v5;
                    end
                    Line1 = [ v1(1:2)  v2(1:2)];
                    Line2 = [ v2(1:2)  v3(1:2)];
                    nPolyPoints = [v3;nPolyPoints;v4]; 
                    Line3 = AW_Fast.GetBoundaryLineList_XY(nPolyPoints);
                    Line4 = [ v4(1:2)  v5(1:2)];
                    Line5 = [ v5(1:2)  v1(1:2)];
                   
                    LineList = [Line1; Line2; Line3; Line4; Line5];
                    velocityEdgeNumbers = 3:size(LineList,1)-2;
                case 12
                    v1 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line4);
                    v2 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line2 );
                    v3 = Vertex3;
                    v4 = Vertex4;
                    nPolyPoints = [v1;nPolyPoints;v2]; 
                    Line1 = AW_Fast.GetBoundaryLineList_XY(nPolyPoints);
                    Line2 = [ v2(1:2)  v3(1:2)];
                    Line3 = [ v3(1:2)  v4(1:2)];
                    Line4 = [ v4(1:2)  v1(1:2)];
                    LineList = [Line1; Line2; Line3;Line4];
                    velocityEdgeNumbers = 1:size(LineList,1)-3;
                case 13
                    v1 = Vertex1;
                    v2 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line1 );
                    v3 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line2 );
                    v4 = Vertex3;
                    v5 = Vertex4;
                    if( isempty(v2) )
                        v2 = v1;
                    end
                    if( isempty( v3) )
                        v3 = v4;
                    end
                    Line1 = [ v1(1:2)  v2(1:2)];
                    nPolyPoints = [v2;nPolyPoints;v3]; 
                    Line2 = AW_Fast.GetBoundaryLineList_XY(nPolyPoints);
                    Line3 = [ v3(1:2)  v4(1:2)];
                    Line4 = [ v4(1:2)  v5(1:2)];
                    Line5 = [ v5(1:2)  v1(1:2)];
                    LineList = [Line1; Line2; Line3; Line4; Line5];
                    velocityEdgeNumbers = 2:size(LineList,1)-3;
                case 14
                    v1 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line1);
                    v2 = Vertex2;
                    v3 = Vertex3;
                    v4 = Vertex4;
                    v5 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line4);
                    if( isempty(v1) )
                        v1 = v2;
                    end
                    if( isempty(v5) )
                        v5 = v4;
                    end
                    Line1 = [ v1(1:2)  v2(1:2)];
                    Line2 = [ v2(1:2)  v3(1:2)];
                    Line3 = [ v3(1:2)  v4(1:2)];
                    Line4 = [ v4(1:2)  v5(1:2)];
                    nPolyPoints = [v5;nPolyPoints;v1]; 
                    Line5 = AW_Fast.GetBoundaryLineList_XY(nPolyPoints);
                    LineList = [Line1; Line2 ;Line3 ;Line4; Line5];
                    velocityEdgeNumbers = 5:size(LineList,1);
                case 15
                    %Interior Cell
                    LineList = [];
            
            end
        end
         
           function [LineList, velocityEdgeNumbers, diagonalLength] = GetApproxPolygonEdges_From_PolyPoints( polyX, polyY, CornerPoints, cellType, nPolyPoints )
            Vertex1 = Vertex( CornerPoints(1,:) );
            Vertex2 = Vertex( CornerPoints(2,:) );
            Vertex3 = Vertex( CornerPoints(3,:) );
            Vertex4 = Vertex( CornerPoints(4,:) );
            Line1 = Line( Vertex1, Vertex2 );
            Line2 = Line( Vertex2, Vertex3 );
            Line3 = Line( Vertex3, Vertex4 );
            Line4 = Line( Vertex4, Vertex1 );
            minX = min(CornerPoints(:,1));
            minY = min(CornerPoints(:,2));
            maxX = max(CornerPoints(:,1));
            maxY = max(CornerPoints(:,2));
            diagonalLength = sqrt( (minX - maxX).^2 + (maxY - minY).^2 );
            
          
            
            %if( AW_Fast.Debug )
                %AW_Fast.DrawPolyLine( [polyX' polyY'], 'c', 1);
                %AW_Fast.DrawPoints( [polyX' polyY'], 'b', 0.5);
            %end
            
            switch cellType
                case 0
                    %Exterior Cell
                    LineList = [];
                case 1
                    v1 = Vertex1;
                    v2 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line1);
                    v3 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line4);
                    Line1 = Line( v1, v2 );
                    nPolyPoints = [v2.getCoords();nPolyPoints;v3.getCoords()]; 
                    Line2 = AW_Fast.GetBoundaryLineList(nPolyPoints);
                    Line3 = Line( v3, v1 );
                    LineList = [Line1 Line2 Line3];
                    velocityEdgeNumbers = 2:length(LineList)-1;
                case 2 
                    v1 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY,  Line1 );
                    v2 = Vertex2;
                    v3 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line2);
                    Line1 = Line( v1, v2 );
                    Line2 = Line( v2, v3 );
                    nPolyPoints = [v3.getCoords();nPolyPoints;v1.getCoords()]; 
                    Line3 = AW_Fast.GetBoundaryLineList(nPolyPoints);
                    LineList = [Line1 Line2 Line3];
                    velocityEdgeNumbers = 3:length(LineList);
                case 3 
                    v1 = Vertex1;
                    v2 = Vertex2;
                    v3 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line2 );
                    v4 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line4 );
                    Line1 = Line( v1, v2 );
                    Line2 = Line( v2, v3 );
                    nPolyPoints = [v3.getCoords();nPolyPoints;v4.getCoords()]; 
                    Line3 = AW_Fast.GetBoundaryLineList(nPolyPoints);
                    Line4 = Line( v4, v1 );
                    LineList = [Line1 Line2 Line3 Line4];
                    velocityEdgeNumbers = 3:length(LineList)-1;
                case 4 
                    v1 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line2);
                    v2 = Vertex3;
                    v3 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line3);
                    Line1 = Line( v1, v2 );
                    Line2 = Line( v2, v3 );
                    nPolyPoints = [v3.getCoords();nPolyPoints;v1.getCoords()]; 
                    Line3 = AW_Fast.GetBoundaryLineList(nPolyPoints);
                    LineList = [Line1 Line2 Line3];
                    velocityEdgeNumbers = 3:length(LineList);
                case 5
                    %Ambiguous case
                    LineList = [];
                case 6 
                    v1 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line1);
                    v2 = Vertex2;
                    v3 = Vertex3;
                    v4 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line3 );
                    Line1 = Line( v1, v2 );
                    Line2 = Line( v2, v3 );
                    Line3 = Line( v3, v4 );
                    nPolyPoints = [v4.getCoords();nPolyPoints;v1.getCoords()]; 
                    Line4 = AW_Fast.GetBoundaryLineList(nPolyPoints);
                    LineList = [Line1 Line2 Line3 Line4];
                    velocityEdgeNumbers = 4:length(LineList);
                case 7
                    v1 = Vertex1;
                    v2 = Vertex2;
                    v3 = Vertex3;
                    v4 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line3);
                    v5 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line4);
                    if( isempty(v4) )
                        v4 = v3;
                    end
                    if( isempty(v5) )
                        v5 = v1;
                    end
                    Line1 = Line( v1, v2 );
                    Line2 = Line( v2, v3 );
                    Line3 = Line( v3, v4 );
                    nPolyPoints = [v4.getCoords();nPolyPoints;v5.getCoords()]; 
                    Line4 = AW_Fast.GetBoundaryLineList(nPolyPoints);
                    Line5 = Line( v5, v1 );
                    LineList = [Line1 Line2 Line3 Line4 Line5];
                    velocityEdgeNumbers = 4:length(LineList)-1;
                case 8
                    v1 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line3);
                    v2 = Vertex4;
                    v3 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line4 );
                    Line1 = Line( v1, v2 );
                    Line2 = Line( v2, v3 );
                    nPolyPoints = [v3.getCoords();nPolyPoints;v1.getCoords()]; 
                    Line3 = AW_Fast.GetBoundaryLineList(nPolyPoints);
                    LineList = [Line1 Line2 Line3];
                    velocityEdgeNumbers = 3:length(LineList);
                case 9
                    v1 = Vertex1;
                    v2 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line1);
                    v3 =AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line3);
                    v4 = Vertex4;
                    Line1 = Line( v1, v2 );
                    nPolyPoints = [v2.getCoords();nPolyPoints;v3.getCoords()]; 
                    Line2 = AW_Fast.GetBoundaryLineList(nPolyPoints);
                    Line3 = Line( v3, v4 );
                    Line4 = Line( v4, v1 );
                    LineList = [Line1 Line2 Line3 Line4];
                   velocityEdgeNumbers = 2:length(LineList)-2;
                case 10
                    %Ambiguous case
                    LineList = [];
                case 11
                    v1 = Vertex1;
                    v2 = Vertex2;
                    v3 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line2);
                    v4 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line3 );
                    v5 = Vertex4;
                     if( isempty(v3) )
                        v3 = v2;
                    end
                    if( isempty(v4) )
                        v4 = v5;
                    end
                    Line1 = Line( v1, v2 );
                    Line2 = Line( v2, v3 );
                    nPolyPoints = [v3.getCoords();nPolyPoints;v4.getCoords()]; 
                    Line3 = AW_Fast.GetBoundaryLineList(nPolyPoints);
                    Line4 = Line( v4, v5 );
                    Line5 = Line( v5, v1 );
                   
                    LineList = [Line1 Line2 Line3 Line4 Line5];
                    velocityEdgeNumbers = 3:length(LineList)-2;
                case 12
                    v1 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line4);
                    v2 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line2 );
                    v3 = Vertex3;
                    v4 = Vertex4;
                    nPolyPoints = [v1.getCoords();nPolyPoints;v2.getCoords()]; 
                    Line1 = AW_Fast.GetBoundaryLineList(nPolyPoints);
                    Line2 = Line( v2, v3 );
                    Line3 = Line( v3, v4 );
                    Line4 = Line( v4, v1 );
                    LineList = [Line1 Line2 Line3 Line4];
                    velocityEdgeNumbers = 1:length(LineList)-3;
                case 13
                    v1 = Vertex1;
                    v2 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line1 );
                    v3 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line2 );
                    v4 = Vertex3;
                    v5 = Vertex4;
                    if( isempty(v2) )
                        v2 = v1;
                    end
                    if( isempty( v3) )
                        v3 = v4;
                    end
                    Line1 = Line( v1, v2 );
                    nPolyPoints = [v2.getCoords();nPolyPoints;v3.getCoords()]; 
                    Line2 = AW_Fast.GetBoundaryLineList(nPolyPoints);
                    Line3 = Line( v3, v4 );
                    Line4 = Line( v4, v5 );
                    Line5 = Line( v5, v1 );
                    LineList = [Line1 Line2 Line3 Line4 Line5];
                    velocityEdgeNumbers = 2:length(LineList)-3;
                case 14
                    v1 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line1);
                    v2 = Vertex2;
                    v3 = Vertex3;
                    v4 = Vertex4;
                    v5 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line4);
                    if( isempty(v1) )
                        v1 = v2;
                    end
                    if( isempty(v5) )
                        v5 = v4;
                    end
                    Line1 = Line( v1, v2 );
                    Line2 = Line( v2, v3 );
                    Line3 = Line( v3, v4 );
                    Line4 = Line( v4, v5 );
                    nPolyPoints = [v5.getCoords();nPolyPoints;v1.getCoords()]; 
                    Line5 = AW_Fast.GetBoundaryLineList(nPolyPoints);
                    LineList = [Line1 Line2 Line3 Line4 Line5];
                    velocityEdgeNumbers = 5:length(LineList);
                case 15
                    %Interior Cell
                    LineList = [];
            
            end
           end
        
             function [LineList, velocityEdgeNumber, diagonalLength] = GetApproxPolygonEdges_XY( polyX, polyY, CornerPoints, cellType )
            Vertex1 = CornerPoints(1,1:2);
            Vertex2 = CornerPoints(2,1:2);
            Vertex3 = CornerPoints(3,1:2);
            Vertex4 = CornerPoints(4,1:2);
            Line1 = [Vertex1, Vertex2];
            Line2 = [Vertex2, Vertex3];
            Line3 =[ Vertex3, Vertex4];
            Line4 = [Vertex4, Vertex1];
            minX = min(CornerPoints(:,1));
            minY = min(CornerPoints(:,2));
            maxX = max(CornerPoints(:,1));
            maxY = max(CornerPoints(:,2));
            diagonalLength = sqrt( (minX - maxX).^2 + (maxY - minY).^2 );
            
            %if( AW_Fast.Debug )
                %AW_Fast.DrawPolyLine( [polyX' polyY'], 'c', 1);
                %AW_Fast.DrawPoints( [polyX' polyY'], 'b', 0.5);
            %end
            
            switch cellType
                case 0
                    %Exterior Cell
                    LineList = [];
                case 1
                    v1 = Vertex1;
                    v2 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line1);
                    v3 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line4);
                    Line1 = [v1(1:2), v2(1:2)];
                    Line2 = [v2(1:2), v3(1:2)];
                    Line3 = [ v3(1:2), v1(1:2)];
                    LineList = [Line1; Line2; Line3];
                    velocityEdgeNumber = 2;
                case 2 
                    v1 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY,  Line1 );
                    v2 = Vertex2;
                    v3 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line2);
                    Line1 = [v1(1:2), v2(1:2)];
                    Line2 = [v2(1:2), v3(1:2)];
                    Line3 = [ v3(1:2), v1(1:2)];
                    LineList = [Line1; Line2; Line3];
                    velocityEdgeNumber = 3;
                case 3 
                    v1 = Vertex1;
                    v2 = Vertex2;
                    v3 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line2 );
                    v4 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line4 );
                    Line1 = [v1(1:2), v2(1:2)];
                    Line2 = [v2(1:2), v3(1:2)];
                    Line3 = [ v3(1:2), v4(1:2)];
                    Line4 = [v4(1:2), v1(1:2)];
                    LineList = [Line1 ;Line2 ;Line3; Line4];
                    velocityEdgeNumber = 3;
                case 4 
                    v1 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line2);
                    v2 = Vertex3;
                    v3 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line3);
                    Line1 = [v1(1:2), v2(1:2)];
                    Line2 = [v2(1:2), v3(1:2)];
                    Line3 = [ v3(1:2), v1(1:2)];
                    LineList = [Line1; Line2; Line3];
                    velocityEdgeNumber = 3;
                case 5
                    %Ambiguous case
                    LineList = [];
                case 6 
                    v1 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line1);
                    v2 = Vertex2;
                    v3 = Vertex3;
                    v4 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line3 );
                    Line1 = [v1(1:2), v2(1:2)];
                    Line2 = [v2(1:2), v3(1:2)];
                    Line3 = [ v3(1:2), v4(1:2)];
                    Line4 = [v4(1:2), v1(1:2)];
                    LineList = [Line1 ;Line2 ;Line3; Line4];
                    velocityEdgeNumber = 4;
                case 7
                    v1 = Vertex1;
                    v2 = Vertex2;
                    v3 = Vertex3;
                    v4 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line3);
                    v5 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line4);
                    if( isempty(v4) )
                        v4 = v3;
                    end
                    if( isempty(v5) )
                        v5 = v1;
                    end
                    Line1 = [v1(1:2), v2(1:2)];
                    Line2 = [v2(1:2), v3(1:2)];
                    Line3 = [ v3(1:2), v4(1:2)];
                    Line4 = [v4(1:2), v5(1:2)];
                    Line5 = [v5(1:2), v1(1:2)];
                    LineList = [Line1; Line2; Line3 ;Line4; Line5];
                    velocityEdgeNumber = 4;
                case 8
                    v1 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line3);
                    v2 = Vertex4;
                    v3 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line4 );
                    Line1 = [v1(1:2), v2(1:2)];
                    Line2 = [v2(1:2), v3(1:2)];
                    Line3 = [ v3(1:2), v1(1:2)];
                    LineList = [Line1; Line2; Line3];
                    velocityEdgeNumber = 3;
                case 9
                    v1 = Vertex1;
                    v2 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line1);
                    v3 =AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line3);
                    v4 = Vertex4;
                    Line1 = [v1(1:2), v2(1:2)];
                    Line2 = [v2(1:2), v3(1:2)];
                    Line3 = [ v3(1:2), v4(1:2)];
                    Line4 = [v4(1:2), v1(1:2)];
                    LineList = [Line1 ;Line2 ;Line3; Line4];
                    velocityEdgeNumber = 2;
                case 10
                    %Ambiguous case
                    LineList = [];
                case 11
                    v1 = Vertex1;
                    v2 = Vertex2;
                    v3 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line2);
                    v4 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line3 );
                    v5 = Vertex4;
                     if( isempty(v3) )
                        v3 = v2;
                    end
                    if( isempty(v4) )
                        v4 = v5;
                    end
                    Line1 = [v1(1:2), v2(1:2)];
                    Line2 = [v2(1:2), v3(1:2)];
                    Line3 = [ v3(1:2), v4(1:2)];
                    Line4 = [v4(1:2), v5(1:2)];
                    Line5 = [v5(1:2), v1(1:2)];
                    LineList = [Line1; Line2; Line3 ;Line4; Line5];
                    velocityEdgeNumber = 3;
                case 12
                    v1 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line4);
                    v2 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line2 );
                    v3 = Vertex3;
                    v4 = Vertex4;
                    Line1 = [v1(1:2), v2(1:2)];
                    Line2 = [v2(1:2), v3(1:2)];
                    Line3 = [ v3(1:2), v4(1:2)];
                    Line4 = [v4(1:2), v1(1:2)];
                    LineList = [Line1 ;Line2 ;Line3; Line4];
                    velocityEdgeNumber = 1;
                case 13
                    v1 = Vertex1;
                    v2 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line1 );
                    v3 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line2 );
                    v4 = Vertex3;
                    v5 = Vertex4;
                    if( isempty(v2) )
                        v2 = v1;
                    end
                    if( isempty( v3) )
                        v3 = v4;
                    end
                    Line1 = [v1(1:2), v2(1:2)];
                    Line2 = [v2(1:2), v3(1:2)];
                    Line3 = [ v3(1:2), v4(1:2)];
                    Line4 = [v4(1:2), v5(1:2)];
                    Line5 = [v5(1:2), v1(1:2)];
                    LineList = [Line1; Line2; Line3 ;Line4; Line5];
                    velocityEdgeNumber = 2;
                case 14
                    v1 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line1);
                    v2 = Vertex2;
                    v3 = Vertex3;
                    v4 = Vertex4;
                    v5 = AW_Fast.GetPtOfIntersection_XY( polyX, polyY, Line4);
                    if( isempty(v1) )
                        v1 = v2;
                    end
                    if( isempty(v5) )
                        v5 = v4;
                    end
                    Line1 = [v1(1:2), v2(1:2)];
                    Line2 = [v2(1:2), v3(1:2)];
                    Line3 = [ v3(1:2), v4(1:2)];
                    Line4 = [v4(1:2), v5(1:2)];
                    Line5 = [v5(1:2), v1(1:2)];
                    LineList = [Line1; Line2; Line3 ;Line4; Line5];
                    velocityEdgeNumber = 5;
                case 15
                    %Interior Cell
                    LineList = [];
            
            end
         end
        
        
         function [LineList, velocityEdgeNumber, diagonalLength] = GetApproxPolygonEdges( polyX, polyY, CornerPoints, cellType )
            Vertex1 = Vertex( CornerPoints(1,:) );
            Vertex2 = Vertex( CornerPoints(2,:) );
            Vertex3 = Vertex( CornerPoints(3,:) );
            Vertex4 = Vertex( CornerPoints(4,:) );
            Line1 = Line( Vertex1, Vertex2 );
            Line2 = Line( Vertex2, Vertex3 );
            Line3 = Line( Vertex3, Vertex4 );
            Line4 = Line( Vertex4, Vertex1 );
            minX = min(CornerPoints(:,1));
            minY = min(CornerPoints(:,2));
            maxX = max(CornerPoints(:,1));
            maxY = max(CornerPoints(:,2));
            diagonalLength = sqrt( (minX - maxX).^2 + (maxY - minY).^2 );
            
            %if( AW_Fast.Debug )
                %AW_Fast.DrawPolyLine( [polyX' polyY'], 'c', 1);
                %AW_Fast.DrawPoints( [polyX' polyY'], 'b', 0.5);
            %end
            
            switch cellType
                case 0
                    %Exterior Cell
                    LineList = [];
                case 1
                    v1 = Vertex1;
                    v2 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line1);
                    v3 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line4);
                    Line1 = Line( v1, v2 );
                    Line2 = Line( v2, v3 );
                    Line3 = Line( v3, v1 );
                    LineList = [Line1 Line2 Line3];
                    velocityEdgeNumber = 2;
                case 2 
                    v1 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY,  Line1 );
                    v2 = Vertex2;
                    v3 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line2);
                    Line1 = Line( v1, v2 );
                    Line2 = Line( v2, v3 );
                    Line3 = Line( v3, v1 );
                    LineList = [Line1 Line2 Line3];
                    velocityEdgeNumber = 3;
                case 3 
                    v1 = Vertex1;
                    v2 = Vertex2;
                    v3 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line2 );
                    v4 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line4 );
                    Line1 = Line( v1, v2 );
                    Line2 = Line( v2, v3 );
                    Line3 = Line( v3, v4 );
                    Line4 = Line( v4, v1 );
                    LineList = [Line1 Line2 Line3 Line4];
                    velocityEdgeNumber = 3;
                case 4 
                    v1 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line2);
                    v2 = Vertex3;
                    v3 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line3);
                    Line1 = Line( v1, v2 );
                    Line2 = Line( v2, v3 );
                    Line3 = Line( v3, v1 );
                    LineList = [Line1 Line2 Line3];
                    velocityEdgeNumber = 3;
                case 5
                    %Ambiguous case
                    LineList = [];
                case 6 
                    v1 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line1);
                    v2 = Vertex2;
                    v3 = Vertex3;
                    v4 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line3 );
                    Line1 = Line( v1, v2 );
                    Line2 = Line( v2, v3 );
                    Line3 = Line( v3, v4 );
                    Line4 = Line( v4, v1 );
                    LineList = [Line1 Line2 Line3 Line4];
                    velocityEdgeNumber = 4;
                case 7
                    v1 = Vertex1;
                    v2 = Vertex2;
                    v3 = Vertex3;
                    v4 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line3);
                    v5 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line4);
                    if( isempty(v4) )
                        v4 = v3;
                    end
                    if( isempty(v5) )
                        v5 = v1;
                    end
                    Line1 = Line( v1, v2 );
                    Line2 = Line( v2, v3 );
                    Line3 = Line( v3, v4 );
                    Line4 = Line( v4, v5 );
                    Line5 = Line( v5, v1 );
                    LineList = [Line1 Line2 Line3 Line4 Line5];
                    velocityEdgeNumber = 4;
                case 8
                    v1 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line3);
                    v2 = Vertex4;
                    v3 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line4 );
                    Line1 = Line( v1, v2 );
                    Line2 = Line( v2, v3 );
                    Line3 = Line( v3, v1 );
                    LineList = [Line1 Line2 Line3];
                    velocityEdgeNumber = 3;
                case 9
                    v1 = Vertex1;
                    v2 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line1);
                    v3 =AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line3);
                    v4 = Vertex4;
                    Line1 = Line( v1, v2 );
                    Line2 = Line( v2, v3 );
                    Line3 = Line( v3, v4 );
                    Line4 = Line( v4, v1 );
                    LineList = [Line1 Line2 Line3 Line4];
                    velocityEdgeNumber = 2;
                case 10
                    %Ambiguous case
                    LineList = [];
                case 11
                    v1 = Vertex1;
                    v2 = Vertex2;
                    v3 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line2);
                    v4 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line3 );
                    v5 = Vertex4;
                     if( isempty(v3) )
                        v3 = v2;
                    end
                    if( isempty(v4) )
                        v4 = v5;
                    end
                    Line1 = Line( v1, v2 );
                    Line2 = Line( v2, v3 );
                    Line3 = Line( v3, v4 );
                    Line4 = Line( v4, v5 );
                    Line5 = Line( v5, v1 );
                   
                    LineList = [Line1 Line2 Line3 Line4 Line5];
                    velocityEdgeNumber = 3;
                case 12
                    v1 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line4);
                    v2 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line2 );
                    v3 = Vertex3;
                    v4 = Vertex4;
                    Line1 = Line( v1, v2 );
                    Line2 = Line( v2, v3 );
                    Line3 = Line( v3, v4 );
                    Line4 = Line( v4, v1 );
                    LineList = [Line1 Line2 Line3 Line4];
                    velocityEdgeNumber = 1;
                case 13
                    v1 = Vertex1;
                    v2 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line1 );
                    v3 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line2 );
                    v4 = Vertex3;
                    v5 = Vertex4;
                    if( isempty(v2) )
                        v2 = v1;
                    end
                    if( isempty( v3) )
                        v3 = v4;
                    end
                    Line1 = Line( v1, v2 );
                    Line2 = Line( v2, v3 );
                    Line3 = Line( v3, v4 );
                    Line4 = Line( v4, v5 );
                    Line5 = Line( v5, v1 );
                    LineList = [Line1 Line2 Line3 Line4 Line5];
                    velocityEdgeNumber = 2;
                case 14
                    v1 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line1);
                    v2 = Vertex2;
                    v3 = Vertex3;
                    v4 = Vertex4;
                    v5 = AW_Fast.GetPtOfIntersection_Basic( polyX, polyY, Line4);
                    if( isempty(v1) )
                        v1 = v2;
                    end
                    if( isempty(v5) )
                        v5 = v4;
                    end
                    Line1 = Line( v1, v2 );
                    Line2 = Line( v2, v3 );
                    Line3 = Line( v3, v4 );
                    Line4 = Line( v4, v5 );
                    Line5 = Line( v5, v1 );
                    LineList = [Line1 Line2 Line3 Line4 Line5];
                    velocityEdgeNumber = 5;
                case 15
                    %Interior Cell
                    LineList = [];
            
            end
         end
        
          function [velocity] = ComputeVelocities_XY( LineList_XY, velocityEdgeNumbers, s, xCurve, yCurve, diagonalLength )
               n = size( velocityEdgeNumbers,2 );
               m = length(s);
               velocity = zeros(m,n);
               for i = 1 : n
                thisLine = LineList_XY(velocityEdgeNumbers(i),:);
                velocity(:,i) = AW_Fast.Velocity_Fast( thisLine, s, xCurve, yCurve, diagonalLength );
               end
        end
        
             
        function [] = ComputeVelocities( LineList, velocityEdgeNumbers, s, xCurve, yCurve, diagonalLength )
               n = length( velocityEdgeNumbers );
               for i = 1 : n
                thisLine = LineList(velocityEdgeNumbers(i));
                v = AW_Fast.Velocity( thisLine, s, xCurve, yCurve, diagonalLength );
                thisLine.SetVelocity(v);
               end
        end
        
        
        function [] = DrawPolygon( Points, color, lineWidth )
          X = [Points(:,1); Points(1,1)];
          Y = [Points(:,2); Points(1,2)];
          plot(X,Y,color,'LineWidth',lineWidth);
        end
        
        function [] = DrawPolyLine( Points, color, lineWidth )
          X = Points(:,1);
          Y = Points(:,2);
          plot(X,Y,color,'LineWidth',lineWidth);
        end
        
        function [] = DrawLines_XY( LineList, color, lineWidth )
            m = length( LineList );
            XY = zeros( m*2, 2 );
            k = 1;
            for i = 1 : m
                XY(k,1) = LineList(i,1);
                XY(k,2) = LineList(i,2);
                k = k + 1;
                XY(k,1) = LineList(i,3);
                XY(k,2) = LineList(i,4);
                k = k + 1;
            end
            AW_Fast.DrawPolyLine( XY, color, lineWidth);
        end
        
        function [] = DrawLines( LineList, color, lineWidth )
            m = length( LineList );
            XY = zeros( m*2, 2 );
            k = 1;
            for i = 1 : m
                vertices = LineList(i).getVertices();
                XY(k,1) = vertices(1).getX();
                XY(k,2) = vertices(1).getY();
                k = k + 1;
                XY(k,1) = vertices(2).getX();
                XY(k,2) = vertices(2).getY();
                k = k + 1;
            end
            AW_Fast.DrawPolyLine( XY, color, lineWidth);
        end
        
        function [] = DrawPoints( Points, color, pointDia )
          X = Points(:,1);
          Y = Points(:,2);
          scatter(X,Y,pointDia,color,'filled');
        end
             
        function [] = DrawArrows( LineList, color, lineWidth)
           m = length( LineList );
             
             X = [];
             Y = [];
             for i = 1 : m
                vertices = LineList(i).getVertices();
                x1 = vertices(1).getX();
                y1 = vertices(1).getY();
                x2 = vertices(2).getX();
                y2 = vertices(2).getY();
                uv = [(x1-x2);(y1-y2)];
                l = norm(uv);
                uv = uv/l;
                n = [-uv(2);uv(1)];
                A = ([x2;y2]+[x1;y1])/2;
                %A = A - n*0.05;
                arrowHeadSize = 0.04*l;
                B = A + uv*arrowHeadSize;
                F = A + uv*(arrowHeadSize)*3;
                C = B + arrowHeadSize*n;
                D = B - arrowHeadSize*n;
                Ax = A(1);
                Ay = A(2);
                Cx = C(1);
                Cy = C(2);
                Dx = D(1);
                Dy = D(2);
                Fx = F(1);
                Fy = F(2);
                X = [X;x1 x2; Ax Cx; Ax Dx; Ax Fx];
                Y = [Y;y1 y2; Ay Cy; Ay Dy; Ay Fy];
             end
             line(X',Y','Color', color, 'LineWidth',lineWidth);
        end
        
        
          function [] = DrawArrows_XY( LineList, color, lineWidth)
           m = length( LineList );
             
             X = [];
             Y = [];
             for i = 1 : m
                x1 = LineList(i,1);
                y1 = LineList(i,2);
                x2 = LineList(i,3);
                y2 = LineList(i,4);
                uv = [(x1-x2);(y1-y2)];
                l = norm(uv);
                uv = uv/l;
                n = [-uv(2);uv(1)];
                A = ([x2;y2]+[x1;y1])/2;
                %A = A - n*0.05;
                arrowHeadSize = 0.04*l;
                B = A + uv*arrowHeadSize;
                F = A + uv*(arrowHeadSize)*3;
                C = B + arrowHeadSize*n;
                D = B - arrowHeadSize*n;
                Ax = A(1);
                Ay = A(2);
                Cx = C(1);
                Cy = C(2);
                Dx = D(1);
                Dy = D(2);
                Fx = F(1);
                Fy = F(2);
                X = [X;x1 x2; Ax Cx; Ax Dx; Ax Fx];
                Y = [Y;y1 y2; Ay Cy; Ay Dy; Ay Fy];
             end
             line(X',Y','Color', color, 'LineWidth',lineWidth);
        end
    end
    
    properties (Constant)
        Debug=false; 
        Fast = true;
    end
    properties (Access = private)
        integrationOrder
        embeddedDomain
        mPt;
    end
    
end

