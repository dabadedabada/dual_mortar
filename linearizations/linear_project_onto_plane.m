function delt_proj_point = linear_project_onto_plane(point, x0, n0, delt_point, delt_x0, delt_n0)
    % Popp 3d (A17)
    % point and delt_point (nodal coordinate linearizations) are 1x3
    % vectors
    
    delt_proj_point = delt_point - (dot(delt_point-delt_x0,n0) + ...
      dot(point-x0, delt_n0))*n0 - dot(point-x0,n0)*delt_n0;
end