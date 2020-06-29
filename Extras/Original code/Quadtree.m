classdef Quadtree
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        m_PMC;
        m_Bbox;
        m_Depth;
        m_RootNode;
        m_ImplicitFn;
    end
    
    methods
         function obj = Quadtree( PMC,implicitFun, bbox )
             obj.m_PMC = PMC;
             obj.m_Bbox = bbox;
             obj.m_ImplicitFn = implicitFun;
             PointA = obj.m_Bbox.left;
             PointC = obj.m_Bbox.right;
             PointB = [PointC(1); PointA(2)];
             PointD = [PointA(1); PointC(2)];
            
             obj.m_RootNode = Node_2D( PMC,implicitFun, PointA, PointB, PointC, PointD);
         end
        
         function obj = Partition( obj, depth )
             obj.m_RootNode.Partition( depth );
         end
         
         function nodes = GetLeafNodes( obj )
             nodeList = NodeList_2D();
             obj.m_RootNode.GetLeafNodes( nodeList );
             nodes = nodeList.m_NodeList;
         end
         
         function [interiorNodes, boundaryNodes] = GetInterior_BoundaryNodes(obj)
             nodes = obj.GetLeafNodes();
             interiorNodes = Quadtree.GetInteriorNodes(nodes);
             boundaryNodes = Quadtree.GetBoundaryNodes(nodes);
         end
         
%          function PlotError( obj )
%              [interiorNodes, boundaryNodes] = GetInterior_BoundaryNodes(obj);
%              nodes = boundaryNodes;
%              [Integral,Error] = Quadtree.GetError_Integral_Vector( nodes );
%              
%              AbsError = abs(Error);
%              TotalError = sum(AbsError);
%              Error_Percentage = 100*AbsError/TotalError;
%              
%              [Sorted_Error_Percent,I] = sort(Error_Percentage);
%              nodes = nodes(I);
%              
%              Sorted_Error = AbsError(I);
%              
%              %Percent = [eps;1;10;20;40;60;80;100];
%              %Color   = [0; 1;10;20 ; 40; 60; 80; 100];
%              Percent = [0 1 5 6:100];
%              Color = Percent;
%              
%              col = zeros( size(nodes) );
%              faces = [];
%              CornerPoints=[];
%              for i = 1 : length(nodes)
%                 thisCoordinates = [nodes(i).m_PointA';nodes(i).m_PointB'; nodes(i).m_PointC';nodes(i).m_PointD'];
%                 CornerPoints = [CornerPoints;thisCoordinates];
%                 index_i = (4*i)-3:4*i;
%                 faces = [faces; index_i];
%                 e_Val = Sorted_Error_Percent(i);
%                 
%                 for j = 1 : length(Percent)
%                     if( e_Val <= Percent(j) )
%                         col(i) = Color(j);
%                         break;
%                     end
%                 end
%                 
%                 str1 = sprintf('%0.2f',e_Val);
%                 str1 = strcat(str1,' %');
%                 str2 = sprintf(' (%0.2e)',Sorted_Error(i));
%                 str = strcat(str1,str2);
%                 MidPt = mean(thisCoordinates);
%                 text(MidPt(1),MidPt(2),str);
%              end
%              
%             hold on;
%             patch('Faces',faces,'Vertices',CornerPoints,'FaceVertexCData',col,'FaceColor','flat','FaceAlpha',0.3);
%             colorbar;
%             hold off;
%          end

         
         
         function DrawTree( obj, color )
             obj.m_RootNode.Draw(color);
         end
    end
    
    methods(Static)
%          function [Integral,Error] = GetError_Integral_Vector( nodes )
%              n = length(nodes);
%              Error = zeros(n,1);
%              Integral = zeros(n,1);
%              for i = 1 : n
%                  Integral(i) = nodes(i).m_Integral;
%                  Error(i) = nodes(i).m_Error;
%              end
%          end
% 
         function [boundaryNodes] = GetBoundaryNodes(nodes)
             n = length(nodes);
             boundaryNodes = [];
             for i = 1 : n
                 if( nodes(i).IsBoundary() )
                    boundaryNodes = [boundaryNodes;nodes(i)];
                 end
             end
         end
         
         function [interiorNodes] = GetInteriorNodes(nodes)
             n = length(nodes);
             interiorNodes = [];
             for i = 1 : n
                 if( nodes(i).IsInterior() )
                    interiorNodes = [interiorNodes;nodes(i)];
                 end
             end
         end
    end
    
end

