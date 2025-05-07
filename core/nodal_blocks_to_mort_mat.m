function mort_mat = nodal_blocks_to_mort_mat(mat_nodal)

if ismatrix(mat_nodal)
  rows = size(mat_nodal, 1);
  cols = size(mat_nodal, 2);
  mort_mat = zeros(3*rows, 3*cols);
  I = eye(3);
  
  for i = 1:rows
      for j = 1:cols
          row_idx = 3*(i-1) + (1:3); 
          col_idx = 3*(j-1) + (1:3);
          mort_mat(row_idx, col_idx) = mat_nodal(i, j) * I;
      end
  end
  return;
elseif ndims(mat_nodal) == 3
  rows = size(mat_nodal, 1);
  cols = size(mat_nodal, 2);
  third_dim = size(mat_nodal, 3);
  mort_mat = zeros(3*rows, 3*cols, third_dim);
  for i = 1:rows
    for j = 1:cols
        row_idx = 3*(i-1) + (1:3); 
        col_idx = 3*(j-1) + (1:3);
        temp = zeros(3,3,third_dim);
        for k=1:3
          temp(k,k,:) = mat_nodal(i, j,:);
        end
        mort_mat(row_idx, col_idx,:) = temp;
    end
  end
  return;
end