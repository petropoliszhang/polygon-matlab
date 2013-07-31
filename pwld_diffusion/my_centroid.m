function [Cx,Cy] = my_centroid( x , y )
% Computes the centroid of a polygon

% aux var
aux = x.*y([2:end 1]) - y.*x([2:end 1]);

% polygon's signed area
signed_area=0.5* sum(aux);

% centroid coordinates
Cx = sum( ( x + x([2:end 1]) ).*aux ) / (6*signed_area);
Cy = sum( ( y + y([2:end 1]) ).*aux ) / (6*signed_area);

end
