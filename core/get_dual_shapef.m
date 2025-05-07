function [Ae, Me] = get_dual_shapef(coo, fe)
  gp_coo = fe.QP;
  weights = fe.QW;
  n_gp = size(gp_coo, 2);  % number of gp
  n_N = size(fe.xi, 1);
  De = zeros(n_N,1);
  Me = zeros(n_N);

  for gp=1:n_gp
    dN_in_gp = fe.dNdxi(gp_coo(:,gp));
    N_in_gp = fe.N(gp_coo(:,gp));
    v1 = coo*dN_in_gp(:,1);
    v2 = coo*dN_in_gp(:,2);
    J_in_gp = norm(cross(v1, v2));
    De = De + N_in_gp*weights(gp)*J_in_gp;    % vector which is made into diagonal later
    Me = Me + N_in_gp*N_in_gp'*weights(gp)*J_in_gp;
  end
  Ae = diag(De)*Me^(-1); 
  
end