function Aalg = linear_get_dual_shapef(coo, nod, fe, Ae, Me, nN, nN_ele)

gauss_points = fe.QP;
weights = fe.QW;
n_gp = size(gauss_points, 2);  % number of gp
Dalg = zeros(nN_ele,nN_ele,3*nN);
Malg = zeros(nN_ele,nN_ele,3*nN);
for gp=1:n_gp
  dN_in_gp = fe.dNdxi(gauss_points(:,gp));
  Ttilde_xi = zeros(3,3*nN);
  Ttilde_eta = zeros(3,3*nN);
  % mapping local node id k to global node id nod(k) to contact id
  for k=1:nN_ele
    con_id =  nod(k)*3-2;
    Ttilde_xi(:,con_id:con_id+2) = eye(3)*dN_in_gp(k,1);
    Ttilde_eta(:,con_id:con_id+2) = eye(3)*dN_in_gp(k,2);
  end
  t_xi = coo*dN_in_gp(:,1);
  t_eta = coo*dN_in_gp(:,2);
  N_in_gp = fe.N(gauss_points(:,gp));
  t_gp = cross(t_eta,t_xi)'*norm(cross(t_eta,t_xi))^(-1)*(skew_mat(t_xi)'*Ttilde_eta + skew_mat(t_eta)*Ttilde_xi);
  for k=1:nN_ele
    Dalg(k,k,:) = Dalg(k,k,:) + reshape(weights(gp)*N_in_gp(k)*t_gp, [1,1,3*nN]);
    for l=1:nN_ele
      Malg(k,l,:) = Malg(k,l,:) + reshape(weights(gp)*N_in_gp(k)*N_in_gp(l)*t_gp, [1,1,3*nN]);
    end
  end
end
Aalg = pagemtimes(Dalg,Me^(-1)) - pagemtimes(Ae,pagemtimes(Malg,Me^(-1)));