function isOnSegment = check_point_on_segment(A, B, P)
    % Check if point P lies on the line segment AB in 2D
    % A, B, and P are 1x2 vectors: [x, y]

    % Compute direction vectors
    AB = B - A;
    AP = P - A;
    
    % Check collinearity using the cross product (should be close to zero)
    crossProd = AB(1) * AP(2) - AB(2) * AP(1);
    
    if abs(crossProd) < 1e-8  % Small tolerance for numerical errors
        % Compute dot product to ensure P is between A and B
        dotProd = dot(AB, AP);
        if dotProd >= 0 && dotProd <= dot(AB, AB)
            isOnSegment = true;
            return;
        end
    end
    
    isOnSegment = false;
end