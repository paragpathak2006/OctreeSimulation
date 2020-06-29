classdef CircularHoleFeature < SmallFeatures
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        m_X_c;
        m_Y_c;
        m_Radius;
        m_Polygon;
        m_ComputeSSACorrection;
        m_SensitivityType;
    end
    
    methods
        function obj = CircularHoleFeature( xc, yc, radius )
            obj.m_X_c = xc;
            obj.m_Y_c = yc;
            obj.m_Radius = radius;
        end
        
        function SetSensitivityType( obj, sType )
            obj.m_SensitivityType = sType;
            
            switch( sType )
                case 0
                    obj.SetSSACorrection(0);
                case 1
                    obj.SetSSACorrection(1);
            end
        end
        
        function sType = GetSensitivityType( obj )
            sType = obj.m_SensitivityType;
        end
        function flag = ImplicitFunction ( obj, X, Y )
            flag = -(X-ones(size(X))*obj.m_X_c).^2 - (Y-ones(size(Y))*obj.m_Y_c).^2 + ones(size(X))*obj.m_Radius.^2;
        end
        function flag = PMC( obj, X , Y )
            flag = 1*(obj.ImplicitFunction(X,Y) >= zeros(size(X)));
        end
        
        function XY = GetCenter( obj )
            XY = [ obj.m_X_c obj.m_Y_c];
        end
        
        function meas = GetMeasure( obj )
            meas = pi*obj.m_Radius.^2;
        end
           
        
        function [X Y] = Polygonize(obj, n)
            t = 0:2*pi/n:2*pi;
            X = (obj.m_Radius.*cos(t))';
            Y = (obj.m_Radius.*sin(t))';
            X = X + ones(size(X)).*obj.m_X_c;
            Y = Y + ones(size(Y)).*obj.m_Y_c;
            obj.m_Polygon = [X Y];
            obj.m_ComputeSSACorrection = zeros(size(X));
        end
        
        function ComputeSSACorrection = GetSSACorrection( obj )
            ComputeSSACorrection = obj.m_ComputeSSACorrection;
        end
        function [] = Display(obj, color)
            k = 1000;
            t = 0:2*pi/k:2*pi;
            X = obj.m_Radius.*cos(t);
            Y = obj.m_Radius.*sin(t);
            X = X + ones(size(X)).*obj.m_X_c;
            Y = Y + ones(size(Y)).*obj.m_Y_c;
            hold on;
            plot(X,Y,color);
            hold off;
        end
        
        function [] = SetSSACorrection( obj, n )
            
            s = size(obj.m_Polygon,1);
            if( n == 0 )
                obj.m_ComputeSSACorrection = zeros(s,1);
            else
                obj.m_ComputeSSACorrection = ones(s,1);
            end
        end
        
        function [] = DisplayNormal( obj, color )
            
            n = size(obj.m_Polygon(:,1),1);
            
            for i = 1 : n-1
                StartPt = obj.m_Polygon(i,:);
                EndPt = obj.m_Polygon(i+1,:);
                n = Cell_2D.GetEdgeNormal(StartPt,EndPt);
                p_c_1 = (StartPt + EndPt)/2;
                p_c_2 = p_c_1' + n*1;
                XY = [p_c_1;p_c_2'];
                hold on;
                plot(XY(:,1),XY(:,2),color);
                hold off;
            end
        end
        
        function [P] = GetPolygon( obj )
            P = obj.m_Polygon;
        end
        function [] = DisplayPolygon(obj, color )
            if(~isempty( obj.m_Polygon ) )
                hold on;
                plot(obj.m_Polygon(:,1),obj.m_Polygon(:,2), color );
                hold off;
            end
        end
        
        function deg = GetDomainDegree(obj)
            deg = 4;
        end
        
        function I = GetAnalyticalIntegral( obj, degree )
               y1 = @(x) obj.m_Y_c+sqrt(obj.m_Radius.^2 - (x-obj.m_X_c).^2);
               y2 = @(x) obj.m_Y_c-sqrt(obj.m_Radius.^2 - (x-obj.m_X_c).^2);
               StartX = obj.m_X_c - obj.m_Radius;
               EndX = obj.m_X_c + obj.m_Radius;
               n = Cell_2D.GetNumberOfBasisFunctions(degree);
               I = zeros(n,1);
               for z = 1 : n
                    A1 = double(Analytical_Integral_NonConvex_Numerical( z, y1, StartX , EndX));
                    A2 = double(Analytical_Integral_NonConvex_Numerical( z, y2, EndX , StartX));
                    I(z) = A1 + A2;
               end

        end
        
        function d = DistanceToBoundary( obj, QP, direction )
            FBeta = @(beta,px,py,nx,ny) (px + beta.*nx - ones(size(px)).*obj.m_X_c).^2 + (py + beta.*ny - ones(size(py)).*obj.m_Y_c).^2 - ones(size(px)).*obj.m_Radius.^2;
            d = DistanceFunction_Circle_General(QP, FBeta, direction,obj.m_X_c,obj.m_Y_c, obj.m_Radius);
        end

    end
    
end

