function Dalg_loc = linear_local_D_mortmat(s_fe, w, jalg, Jcell, sl_proj_gp, G)

N_in_gp = s_fe.N(sl_proj_gp);
nN_ele = length(N_in_gp);
size_con = size(G,2);
dN_in_gp = s_fe.dNdxi(sl_proj_gp);

Dalg_loc = zeros(nN_ele,nN_ele,size_con);
for k=1:nN_ele
  Dalg_loc(k,k,:) = w*N_in_gp(k)*jalg + w*Jcell*dN_in_gp(k,:)*G(1:2,:);
end

