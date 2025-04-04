function sum_normals = get_sum_normals(mesh, face, normals_list)
% averages normals from different elements on each node
% returns (nN*nele)x3 avg_normals - averaged node normals sorted by element
% (e.g. first nN rows will be normals in nodes (in counterclockwise order) of first element)
    node_list = reshape(mesh.bou{face}.nod', [], 1);
    unique_nodes = unique(node_list); % sorted unique nodes
    sum_normals = zeros(size(normals_list));
    normals_list = normalize(normals_list, 2);
    
    % Average normals per unique node
    for k = 1:length(unique_nodes)
      node_id = unique_nodes(k);
      indices = find(node_list == node_id); % Find all occurrences
      same_node_normals = normals_list(indices, :); % normals on sorted nodes
      avg_normal = mean(same_node_normals, 1); % Average normals
      for j=1:size(indices)
        avg_normals(indices(j), :) = avg_normal / norm(avg_normal); % Normalize
      end
    end