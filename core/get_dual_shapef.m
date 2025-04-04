function [Ae, Me] = get_dual_shapef(ele_nod, ele_fe)
  fe = ele_fe;
  gp = fe.QP;
  weights = fe.QW;
  n_gp = size(gp, 2);  % number of gp
  dN_in_gp = fe.dNdxi(gp);
  N_in_gp = fe.N(gp);
  n_N = size(N_in_gp, 1);
  De = zeros(n_N,1);  % vzorec 25
  Me = zeros(n_N);


  for g=1:n_gp
    v1 = ele_nod'*dN_in_gp(:,g);
    v2 = ele_nod'*dN_in_gp(:,n_gp+g);
    J_in_gp = norm(cross(v1, v2));
    De = De + N_in_gp(:,g)*weights(g)*J_in_gp;
    Me = Me + N_in_gp(:,g)*N_in_gp(:,g)'*weights(g)*J_in_gp;
  end
  Ae = diag(De)*Me^(-1); % vzorec (25)
  
end