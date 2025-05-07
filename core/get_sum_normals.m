function sum_normals = get_sum_normals(cont_face, ele_normals)
% sums normals from different elements on each node
% returns (nN*nele)x3 sum_normals - summed node normals sorted by element
% (e.g. first nN rows will be normals in nodes (in counterclockwise order) of first element)
nN = cont_face.info.nN;         % Number of Nodes on the face

sum_normals = zeros(nN, 3);

for i=1:nN
  [ele_id, xi_id] = find(cont_face.nod==i);
  for e = 1:length(ele_id)
    curr_normal = ele_normals{ele_id(e)}(xi_id(e),:);
    sum_normals(i,:) = sum_normals(i,:) + curr_normal/norm(curr_normal);
  end
end


%{
node_list = reshape(cont_face.nod', [], 1);
unique_nodes = unique(node_list); % sorted unique nodes
sum_normals = zeros(nN, 3);
unit_normals_list = zeros(nele*nN_ele, 3);
for i=1:nele
  for j=1:nN_ele
    unit_normals_list((i-1)*nN_ele+j,:) = ele_normals{i}(j,:)/norm(ele_normals{i}(j,:));
  end
end

% Average normals per unique node
for k = 1:nN
  node_id = unique_nodes(k);
  indices = node_list == node_id; % Find all occurrences
  same_node_normals = unit_normals_list(indices, :); % normals on sorted nodes
  % returns list of normals
  sum_normals(k, :) = sum(same_node_normals, 1); 
end
%}