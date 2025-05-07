function G = linear_proj_gp(coo, nod, fe, proj_gp_coo, alpha, n0, N0, Ghat, sl_or_m, nN_m)

nN_ele = length(nod);

Gtilde = zeros(size(N0));
N_in_gp = fe.N(proj_gp_coo);
% mapping local node id k to global node id nod(k) to contact id
for k=1:nN_ele
  if (sl_or_m == 1)
    con_id = nN_m*3 + nod(k)*3-2;
  elseif (sl_or_m == 2)
    con_id = nod(k)*3-2;
  end
  Gtilde(:,con_id:con_id+2) = eye(3)*N_in_gp(k);
end

J = newton_it_J(coo, n0, fe, proj_gp_coo);

G = -J^(-1)*(Gtilde - alpha*N0-Ghat);


