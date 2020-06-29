classdef Cell_2D < handle
    %UNTITLED9 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        m_ApproxPolygonVertices;
        m_DistanceToBoundaryFun;
        m_DomainDegree;
        m_nEdges;
        m_Adaptive1DQuadrature;
        m_Integral;
        m_Error;
    end
    
    methods
        function obj = Cell_2D( ApproxPolygon, DistanceToBoundary, DomainDegree, AdaptiveQuadrature )
            obj.m_ApproxPolygonVertices = ApproxPolygon;
            obj.m_DistanceToBoundaryFun = DistanceToBoundary;
            obj.m_DomainDegree = DomainDegree;
            obj.m_nEdges = size(obj.m_ApproxPolygonVertices,1)-1;
            obj.m_Adaptive1DQuadrature = AdaptiveQuadrature;
        end
        
        function A = GetAMatrix_t( obj, t, StartPt, EndPt, degree )
           [n1,m1] = size(t);
           X_0 = ones(n1,m1).*StartPt(1);
           Y_0 = ones(n1,m1).*StartPt(2);
           E = Cell_2D.GetEdgeDirection(StartPt,EndPt);
           X = X_0 + t.*E(1);
           Y = Y_0 + t.*E(2);
           XY = [X' Y'];
           hold on;
           %scatter( X , Y );
           hold off;
           A = Cell_2D.GetAMatrix(XY, degree);
        end
        
        function Ax = GetA_x_t( obj, t, StartPt, EndPt, degree )
           [n1,m1] = size(t);
           X_0 = ones(n1,m1).*StartPt(1);
           Y_0 = ones(n1,m1).*StartPt(2);
           E = Cell_2D.GetEdgeDirection(StartPt,EndPt);
           X = X_0 + t.*E(1);
           Y = Y_0 + t.*E(2);
           XY = [X' Y'];
           hold on;
           %scatter( X , Y );
           hold off;
           Ax = Cell_2D.GetA_Grad_X(XY, degree);
        end
        
        function Ay = GetA_y_t( obj, t, StartPt, EndPt, degree )
           [n1,m1] = size(t);
           X_0 = ones(n1,m1).*StartPt(1);
           Y_0 = ones(n1,m1).*StartPt(2);
           E = Cell_2D.GetEdgeDirection(StartPt,EndPt);
           X = X_0 + t.*E(1);
           Y = Y_0 + t.*E(2);
           XY = [X' Y'];
           hold on;
           %scatter( X , Y );
           hold off;
           Ay = Cell_2D.GetA_Grad_Y(XY, degree);
        end
        
        
        
        function d = GetDistanceToBoundary( obj, StartPt, EndPt, t)
            n = length(t);
            X_0 = StartPt(1);
            Y_0 = StartPt(2);
            d = zeros(n,1);
            EdgeDirection = Cell_2D.GetEdgeDirection(StartPt,EndPt);
            EdgeNormal = Cell_2D.GetEdgeNormal( StartPt, EndPt );    
            for i = 1 : n
                X = X_0 + t(i).*EdgeDirection(1);
                Y = Y_0 + t(i).*EdgeDirection(2);
                %d(i) = (obj.m_DistanceToBoundaryFunction( [X Y], EdgeNormal )*[0;1])'*EdgeNormal; 
                d(i) = obj.m_DistanceToBoundaryFun( [X Y], EdgeNormal ); 
                X2 = X + d(i).*EdgeNormal(1);
                Y2 = Y + d(i).*EdgeNormal(2);
                hold on;
                line([X X2],[Y Y2],'Color','r');
                Cell_2D.PlotLine(X,Y,X2,Y2);
                hold off
            end
        end
        
        function I_i = GetSSATerm( obj, i, A, d, t )
            AMatrix = A(t);
            n = size(AMatrix,1);
            dist = d(t)';
            D = repmat(dist,n,1);
            AD = AMatrix.*D;
            I_i = AD(i,:);
            
        end
        
        function C1 = GetC1( obj, StartPt, EndPt, degree )
            AMatrix = @(t) obj.GetAMatrix_t( t, StartPt, EndPt, degree );
            Distance = @(t) obj.GetDistanceToBoundary( StartPt, EndPt, t);
            Length = Cell_2D.GetEdgeLength(StartPt,EndPt);
            nTerms = Cell_2D.GetNumberOfBasisFunctions( degree );
            C1 = zeros(nTerms,1);
            for i = 1 : nTerms
                Integrand = @(t) obj.GetSSATerm( i, AMatrix, Distance, t );
                C1(i) = integral(Integrand,0,Length);
            end
        end
        
        function I_i = GetSSATerm_2( obj, i, Ax, Ay, d, t, n )
            A_x = Ax(t);
            A_y = Ay(t);
            AMatrix = A_x*n(1) + A_y*n(2);
            n = size(AMatrix,1);
            dist = d(t).^2';
            D = repmat(dist,n,1);
            AD = 0.5*AMatrix.*D;
            I_i = AD(i,:);
        end
        
        function C1 = GetC2( obj, StartPt, EndPt, degree )
            [Ax] = @(t) obj.GetA_x_t( t, StartPt, EndPt, degree );
            [Ay] = @(t) obj.GetA_y_t( t, StartPt, EndPt, degree );
            Distance = @(t) obj.GetDistanceToBoundary( StartPt, EndPt, t);
            Length = Cell_2D.GetEdgeLength(StartPt,EndPt);
            nTerms = Cell_2D.GetNumberOfBasisFunctions( degree );
            C1 = zeros(nTerms,1);
            n = Cell_2D.GetEdgeNormal(StartPt,EndPt);
            for i = 1 : nTerms
                Integrand = @(t) obj.GetSSATerm_2( i, Ax, Ay, Distance, t, n );
                C1(i) = integral(Integrand,0,Length);
            end
        end
        
        
%         function [PolygonalArea,C1] = GetC1_Actual( obj, StartPt, EndPt )
%             base = StartPt(1) - EndPt(1);
%                 h1 = StartPt(2);
%                 h2 = EndPt(2);
%                 x1 = EndPt(1);
%                 x2 = StartPt(1);
%                 PolygonalArea = 0.5*base*(h1 + h2);
%                 TotalArea = integral(obj.m_yFun,x1,x2);
%                 C1 = TotalArea - PolygonalArea;
%         end
%         
       function [lhs_1,Moment_Error] = GetLHS( obj, degree, ComputeSSACorrection, SmallFeaturesList )
           nPt = Cell_2D.GetNumOfIntegrationPoints_1D( degree + 1 ) ;
           [ W , t ] = Cell_2D.Gauss_W_X_Scaled( 0,1, nPt );
           nBasisFunctions = Cell_2D.GetNumberOfBasisFunctions( degree );
           lhs_1 = zeros(nBasisFunctions,1);
           Moment_Error = zeros(nBasisFunctions, 1);
           PArea_Computed = 0;
           PArea_Actual = 0;
           for i = 1 : obj.m_nEdges 
               StartPt = obj.m_ApproxPolygonVertices(i,:)';
               EndPt = obj.m_ApproxPolygonVertices(i+1,:)';
               L = Cell_2D.GetEdgeLength(StartPt,EndPt);
               if( L <= eps )
                   continue;
               end
               X = StartPt(1) + t.*(EndPt(1) - StartPt(1));
               Y = StartPt(2) + t.*(EndPt(2) - StartPt(2));
               XY = [X Y];
               %Cell_2D.PlotPoints( X, Y );
               dy = EndPt(2) - StartPt(2);
               Beta = Cell_2D.GetBetaMatrix(XY, degree);
               PolygonalIntegral = dy*Beta*W;
               lhs_1 = lhs_1 + PolygonalIntegral;
               PArea_Computed = PArea_Computed + PolygonalIntegral;
               
               if( ComputeSSACorrection(i) )
                %Compute 1st Order SSA Correction
                if(~obj.m_Adaptive1DQuadrature)
                nPt1 = Cell_2D.GetNumOfIntegrationPoints_1D( degree + obj.m_DomainDegree ) ;
                
                [ W1 , t1 ] = Cell_2D.Gauss_W_X_Scaled( 0,1, nPt1 );
                X = StartPt(1) + t1.*(EndPt(1) - StartPt(1));
                Y = StartPt(2) + t1.*(EndPt(2) - StartPt(2));
                XY = [X Y];
                EdgeNormal = Cell_2D.GetEdgeNormal( StartPt, EndPt );    
                A = Cell_2D.GetAMatrix(XY, degree);
                d = zeros(size(W1,1),1);
                for j = 1 : nPt1
                    d(j) = obj.m_DistanceToBoundaryFun( XY(j,:), EdgeNormal ); 
                    pointOnCurve = XY(j,:)' + d(j)*EdgeNormal;
                    %pointOnCurve = XY(j,:)' + d(j)*[0;1];
                    Cell_2D.PlotLine(XY(j,1), XY(j,2),pointOnCurve(1), pointOnCurve(2));
                end
                C1 = A*diag(d)*W1*L;
                %[P_Actual,C1_Actual]= obj.GetC1_Actual( StartPt, EndPt );
                %Error(i) = abs(C1_Actual - C1);
                lhs_1 = lhs_1 + C1;
                %PArea_Actual = PArea_Actual + P_Actual;
                else
                  C1 = obj.GetC1( StartPt, EndPt, degree );
                  %[P_Actual,C1_Actual] = obj.GetC1_Actual( StartPt, EndPt );
                  %Error(i) = abs(C1_Actual - C1);
                  lhs_1 = lhs_1 + C1;
                  Moment_Error = Moment_Error + obj.GetC2( StartPt, EndPt, degree);
                  %PArea_Actual = PArea_Actual + P_Actual;
                end
               end
           end
           if(~isempty(SmallFeaturesList) )
               [C_S,e_S] = obj.GetC_S( degree, SmallFeaturesList );
               lhs_1 = lhs_1 - C_S;
               Moment_Error = Moment_Error + e_S;
           end
           %Error
           PolygonalArea_Error = abs(PArea_Actual - PArea_Computed);
       end
       
       function [lhs_1, e_1] = GetLHS_TS( obj, SmallFeature, polynomialDegree )
           XY = SmallFeature.GetCenter();
           Area = SmallFeature.GetMeasure();
           A = Cell_2D.GetAMatrix(XY, polynomialDegree);
           lhs_1 = A*Area;
           e_1 = zeros(size(lhs_1));
       end
       
       function [C_s,e_s] = GetC_S( obj, polynomialDegree, SmallFeaturesList )
           if( isempty(SmallFeaturesList) )
               C_s = [];
               return;
           end
               
           nSmallFeatures = length(SmallFeaturesList);

                for i = 1 : nSmallFeatures
                  SensitivityType = SmallFeaturesList{i}.GetSensitivityType();
                  if( SensitivityType == 1 ||  SensitivityType == 0)
                    ApproxPolygon = SmallFeaturesList{i}.GetPolygon();
                    ComputeSSACorrection = SmallFeaturesList{i}.GetSSACorrection();
                    dist = @(QP,Direction) SmallFeaturesList{i}.DistanceToBoundary(QP,Direction);
                    AdaptiveQuadrature = false;
                    domainDegree = SmallFeaturesList{i}.GetDomainDegree();
                    thisCell = Cell_2D(ApproxPolygon, dist, domainDegree, AdaptiveQuadrature);
                    [lhs_1, e_1] = thisCell.GetLHS( polynomialDegree, ComputeSSACorrection, [] );
                  else
                      [lhs_1,e_1] = obj.GetLHS_TS( SmallFeaturesList{i},polynomialDegree );
                  end
                    if( i == 1 )
                        C_s = lhs_1;
                        e_s = e_1;
                    else
                        C_s = C_s + lhs_1;
                        e_s = e_s + e_1;
                    end
                end
     
    end
       
       
       
       function [I, e] = GetIntegral( obj, IntegrationPoints_XY, PolynomialDegree, ComputeSSACorrection, integrand, SmallFeaturesList )
           A = Cell_2D.GetAMatrix( IntegrationPoints_XY, PolynomialDegree );
           [L,M_e] = obj.GetLHS( PolynomialDegree, ComputeSSACorrection, SmallFeaturesList  )  ;
           W = A \ L;
           f = integrand(IntegrationPoints_XY(:,1),IntegrationPoints_XY(:,2));
           I = W'*f;
           z = A'\f;
           e = M_e' * z;
           obj.m_Integral = I;
           obj.m_Error = e;
       end
       
       function [W] = GetAdaptiveWeights( obj, IntegrationPoints_XY, PolynomialDegree, ComputeSSACorrection )
           A = Cell_2D.GetAMatrix( IntegrationPoints_XY, PolynomialDegree );
           [L,~] = obj.GetLHS( PolynomialDegree, ComputeSSACorrection, []  )  ;
           W = A \ L;
       end
       
          
    end
    
    methods(Static)
        function nBasis = GetNumberOfBasisFunctions( PolynomialDegree )
            level = PolynomialDegree + 1;
            nBasis = level*(level+1)/2;
        end
        
        function PlotPoints( X , Y )
         hold on;
           scatter(X,Y);
         hold off;
        end
        
        function PlotLine( x1 , y1, x2, y2 )
            X = [x1 x2];
            Y = [y1 y2];
         hold on;
            line(X,Y);
         hold off;
        end
        
        function nPoints = GetNumberOfPoints_2D( PolynomialDegree )
            nPoints = Cell_2D.GetNumberOfBasisFunctions( PolynomialDegree );
        end
        
        function Length = GetEdgeLength( StartPt, EndPt )
            EdgeDirection = (EndPt - StartPt);
            Length = norm(EdgeDirection);
        end
        
        function [EdgeDirection] = GetEdgeDirection( StartPt, EndPt )
               EdgeDirection = (EndPt - StartPt);
               Length = norm(EdgeDirection);
               if( Length <= eps )
                   EdgeDirection = [ 0 0];
                   return;
               end
               EdgeDirection = EdgeDirection / Length;
        end
        
        function [EdgeNormal] = GetEdgeNormal( StartPt, EndPt )
            EdgeDirection = Cell_2D.GetEdgeDirection( StartPt, EndPt );
            EdgeNormal = [EdgeDirection(2); -EdgeDirection(1)];
        end
        
        function nDegree = GetNumOfIntegrationPoints_1D( polynomialDegree )
            nDegree = ceil((polynomialDegree + 1)/2);
        end
        
           function [ W , X ] = Gauss_W_X_Scaled( a,b, nPt )
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
          
          function [C] = GetA_Grad_Y( XY, polynomialDegree )
            
            C = [];
            nPt = size(XY,1);
            zer = zeros(1,nPt);
            one = ones(1,nPt);
            X = XY(:,1)';
            Y = XY(:,2)';
            if( polynomialDegree >= 0 )
                %A = [A; Ones];
                C= [C; zer];
            end
            
            if( polynomialDegree >= 1 )
                %A = [A; X; Y;];
                C = [C; zer;one];
            end
            
            if( polynomialDegree >= 2 )
                %A = [A; X.^2; X.*Y; Y.^2;];
                C = [C; zer; X; 2.*Y];
            end
            
            if( polynomialDegree >= 3 )
                %A = [A; X.^3; X.^2.*Y; Y.^2.*X; Y.^3;];
                C = [C; zer; X.^2; 2*Y.*X; 3*Y.^2];
            end
            
            if( polynomialDegree >= 4 )
                %A = [A; X.^4; X.^3.*Y; X.^2.*Y.^2; Y.^4; Y.^3.*X;];
                C = [C; zer; X.^3; 2*X.^2.*Y; 4*Y.^3; 3*Y.^2.*X];
            end
            
            if( polynomialDegree >= 5 )
                %A = [A; X.^5; X.^4.*Y; X.^3.*Y.^2; Y.^5; Y.^4.*X; Y.^3.*X.^2;];
                C = [C; zer; X.^4; 2*X.^3.*Y; 5*Y.^4; 4*Y.^3.*X; 3*Y.^2.*X.^2;];
            end
            
            if( polynomialDegree >= 6 )
                %A = [A;  X.^6; X.^5.*Y; X.^4.*Y.^2; X.^3.*Y.^3; Y.^6;  Y.^5.*X; Y.^4.*X.^2;];
                C = [C;  zer; X.^5; 2*X.^4.*Y; 3*X.^3.*Y.^2; 6*Y.^5;  5*Y.^4.*X; 4*Y.^3.*X.^2;];
            end
            
            if( polynomialDegree >=7 )
                disp('Error: A Matrix not implemented for level >= 7' );
                return;
            end
          end

           
          function [B] = GetA_Grad_X( XY, polynomialDegree )
            
            B = [];
            nPt = size(XY,1);
            zer = zeros(1,nPt);
            one = ones(1,nPt);
            X = XY(:,1)';
            Y = XY(:,2)';
            if( polynomialDegree >= 0 )
                %A = [A; Ones];
                B = [B; zer];
            end
            
            if( polynomialDegree >= 1 )
                %A = [A; X; Y;];
                B = [B; one; zer];
            end
            
            if( polynomialDegree >= 2 )
                %A = [A; X.^2; X.*Y; Y.^2;];
                B = [B; 2.*X; Y; zer];
            end
            
            if( polynomialDegree >= 3 )
                %A = [A; X.^3; X.^2.*Y; Y.^2.*X; Y.^3;];
                B = [B; 3*X.^2; 2*X.*Y; Y.^2; zer];
            end
            
            if( polynomialDegree >= 4 )
                %A = [A; X.^4; X.^3.*Y; X.^2.*Y.^2; Y.^4; Y.^3.*X;];
                B = [B; 4*X.^3; 3*X.^2.*Y; 2*X.*Y.^2; zer; Y.^3];
            end
            
            if( polynomialDegree >= 5 )
                %A = [A; X.^5; X.^4.*Y; X.^3.*Y.^2; Y.^5; Y.^4.*X; Y.^3.*X.^2;];
                B = [B;5*X.^4; 4*X.^3.*Y; 3*X.^2.*Y.^2; zer; Y.^4; 2*Y.^3.*X;];
            end
            
            if( polynomialDegree >= 6 )
                %A = [A;  X.^6; X.^5.*Y; X.^4.*Y.^2; X.^3.*Y.^3; Y.^6;  Y.^5.*X; Y.^4.*X.^2;];
                B = [B;  6*X.^5; 5*X.^4.*Y; 4*X.^3.*Y.^2; 3*X.^2.*Y.^3; zer;  Y.^5; 2*Y.^4.*X;];
            end
            
            if( polynomialDegree >=7 )
                disp('Error: A Matrix not implemented for level >= 7' );
                return;
            end
          end

           
           
           
          function A = GetAMatrix( XY, polynomialDegree )
            
            A = [];
            nPt = size(XY,1);
            one = ones(1,nPt);
            X = XY(:,1)';
            Y = XY(:,2)';
            if( polynomialDegree >= 0 )
                A = [A; one];
            end
            
            if( polynomialDegree >= 1 )
                A = [A; X; Y;];
            end
            
            if( polynomialDegree >= 2 )
                A = [A; X.^2; X.*Y; Y.^2;];
            end
            
            if( polynomialDegree >= 3 )
                A = [A; X.^3; X.^2.*Y; Y.^2.*X; Y.^3;];
            end
            
            if( polynomialDegree >= 4 )
                A = [A; X.^4; X.^3.*Y; X.^2.*Y.^2; Y.^4; Y.^3.*X;];
            end
            
            if( polynomialDegree >= 5 )
                A = [A; X.^5; X.^4.*Y; X.^3.*Y.^2; Y.^5; Y.^4.*X; Y.^3.*X.^2;];
            end
            
            if( polynomialDegree >= 6 )
                A = [A; 
                    X.^6; 
                    X.^5.*Y; 
                    X.^4.*Y.^2; 
                    X.^3.*Y.^3; 
                    Y.^6;
                    Y.^5.*X;
                    Y.^4.*X.^2;
                    ];
                    
            end
            
            if( polynomialDegree >=7 )
                disp('Error: A Matrix not implemented for level >= 7' );
                return;
            end
          end
        function A = GetBetaMatrix( XY, polynomialDegree )
            
            A = [];
            X = XY(:,1)';
            Y = XY(:,2)';
            if( polynomialDegree >= 0 )
                A = [A; X];
            end
            
            if( polynomialDegree >= 1 )
                A = [A; X.^2 ./2; X.*Y;];
            end

            if( polynomialDegree >= 2 )
                A = [A; X.^3 ./3; (X.^2 ./ 2).*Y;  X.*Y.^2; ];
            end
              %A = [A; X.^3; X.^2.*Y; X.^2.*Z; X.*Y.*Z; Y.^3; Y.^2.*X; Y.^2.*Z; Z.^3; Z.^2.*X; Z.^2.*Y];
              % A = [A; X.^3; X.^2.*Y; Y.^3; Y.^2.*X; ];

            if( polynomialDegree >= 3 )
                A = [A; X.^4 ./ 4; (X.^3 ./3).*Y; Y.^2.*(X.^2 ./2);X.*Y.^3;];
            end
            
            if( polynomialDegree >= 4 )
                A = [A; X.^5 ./5; (X.^4 ./ 4).*Y; (X.^3 ./3).*Y.^2; X.*Y.^4; Y.^3.*(X.^2 /2);];
            end
            
            if( polynomialDegree >= 5 )
                A = [A; X.^6./6; (X.^5./5).*Y; (X.^4./4).*Y.^2; X.*Y.^5; Y.^4.*(X.^2./2); Y.^3.*(X.^3./3);];
            end
            
            if( polynomialDegree >= 6 )
               A = [A; X.^7/7; (X.^6.*Y)/6; (X.^5.*Y.^2)/5; (X.^4.*Y.^3)/4; X.*Y.^6; (X.^2.*Y.^5)/2; (X.^3.*Y.^4)/3; ];
            end
            
            
            if( polynomialDegree >=7 )
                disp('Error: Beta Matrix not implemented for level >= 5' );
                return;
            end
        end
        
    end
    
end

