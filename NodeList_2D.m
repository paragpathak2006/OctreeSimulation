classdef NodeList_2D < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        m_NodeList;
    end
    
    methods
            function obj = NodeList_2D()
            end
            
            function Add( obj, thisNode )
                obj.m_NodeList = [obj.m_NodeList;thisNode];
            end
            
            function nodeList = GetNodeList( obj )
                nodeList = obj.m_NodeList;
            end
    end
    
end

