function projPoints = project_onto_plane(points, x0, n0)
    % points are Nx3
    n0 = n0 / norm(n0);  % Ensure normal is unit length
    dists = (points - x0) * n0';  % Compute distances along n0
    projPoints = points - dists .* n0;  % Project points
end