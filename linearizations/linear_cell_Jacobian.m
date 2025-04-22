function jcell = linear_cell_Jacobian(v_coo, V1, V2, V3)
v1_coo = v_coo(:,1);
v2_coo = v_coo(:,2);
v3_coo = v_coo(:,3);
jcell = norm(cross(v2_coo-v1_coo,v3_coo-v1_coo))^(-1)*cross(v2_coo-v1_coo,v3_coo-v1_coo)'*...
  (skew_mat(v3_coo-v1_coo)'*(V2-V1)+skew_mat(v2_coo-v1_coo)*(V3-V1));