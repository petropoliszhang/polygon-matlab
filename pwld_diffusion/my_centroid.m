function [C] = my_centroid( vert, centroid_opt )
% Computes the centroid of a polygon
% global centroid_opt

switch centroid_opt
    case{1}
        C=mean(vert);
    case{2}
        % aux var
        x=vert(:,1);
        y=vert(:,2);
        aux = x.*y([2:end 1]) - y.*x([2:end 1]);
        
        % polygon's signed area
        signed_area=0.5* sum(aux);
        
        % centroid coordinates
        C(1) = sum( ( x + x([2:end 1]) ).*aux ) / (6*signed_area);
        C(2) = sum( ( y + y([2:end 1]) ).*aux ) / (6*signed_area);
    otherwise
        error('centroid option unknown');
end

return
end
