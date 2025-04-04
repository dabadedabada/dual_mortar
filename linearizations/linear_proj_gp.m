function delt_proj_gp = linear_proj_gp(ele_nod, delt_ele_nod, ele_fe, proj_gp_coo, alphas, n0, delt_n0, delt_cell_gp)

n_gp = size(proj_gp_coo, 1);
delt_proj_gp = zeros(n_gp, 3);

for i=1:n_gp
  Ncell_in_gp = ele_fe.N(proj_gp_coo(i,:)');
  dN_in_gp = ele_fe.dNdxi(proj_gp_coo(i,:)');
  l1 = dN_in_gp(:,1)'*ele_nod;
  l2 = dN_in_gp(:,2)'*ele_nod;
  L = [l1', l2', -n0'];
  term = alphas(i)*delt_n0+delt_cell_gp(i,:)-Ncell_in_gp'*delt_ele_nod;
  result = L^(-1)*term';
  delt_proj_gp(i,:)=result';
end
