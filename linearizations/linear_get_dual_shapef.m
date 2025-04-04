function delt_A = linear_get_dual_shapef(ele_nod, delt_ele_nod, ele_fe, Ae, Me)

  fe = ele_fe;
  gp = fe.QP;
  weights = fe.QW;
  n_gp = size(gp, 2);  % number of gp
  dN_in_gp = fe.dNdxi(gp);
  N_in_gp = fe.N(gp);
  n_N = size(N_in_gp, 1);
  delt_De = zeros(n_N,1);  % (A24)
  delt_Me = zeros(n_N);


  for g=1:n_gp
    v1 = ele_nod'*dN_in_gp(:,g); % deriv in the xi direction Popp 3D (A26)
    delt_v1 = delt_ele_nod'*dN_in_gp(:,g);
    v2 = ele_nod'*dN_in_gp(:,n_gp+g); % deriv in the eta direction
    delt_v2 = delt_ele_nod'*dN_in_gp(:,n_gp+g);

    delt_J_in_gp = dot(cross(v1, v2)/norm(cross(v1, v2)),cross(delt_v1, v2)+cross(v1, delt_v2));
    delt_De = delt_De + N_in_gp(:,g)*weights(g)*delt_J_in_gp;
    delt_Me = delt_Me + N_in_gp(:,g)*N_in_gp(:,g)'*weights(g)*delt_J_in_gp;
  end
  delt_A = diag(delt_De)*Me^(-1)-Ae*delt_Me*Me^(-1); % Popp 3D (A23)
  
end