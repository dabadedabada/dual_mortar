function Py = project_onto_plane(y, x0, n0)
    % points are Nx3
    nN_ele = size(y,2);
    Py = zeros(size(y));
    for k=1:nN_ele
      alpha = (y(:,k)- x0)' * n0;  % Compute distances along n0
      Py(:,k) = y(:,k) - alpha * n0;  % Project points
    end
end