classdef SmallFeatures < handle
   methods(Abstract)
      flag = PMC( obj, X , Y );
      d = DistanceToBoundary( obj, QP, direction );
      [X Y] = Polygonize(obj, n);
      Display(obj, color);
      DisplayPolygon(obj, color );
   end
end