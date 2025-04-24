function ghat = linear_disc_gapf(Ngp, normal_gp, G_m, G_s, s_proj_gp, m_proj_gp, s_coo, m_coo, s_nod, m_nod, s_fe, m_fe, nN_m)

C_s = zeros(size(G_s));
C_m = zeros(size(G_m));

nN_ele_s = length(s_nod);
nN_ele_m = length(m_nod);

dNs_in_sgp = s_fe.dNdxi(s_proj_gp);
dNm_in_mgp = m_fe.dNdxi(m_proj_gp);

Ns_in_sgp = s_fe.N(s_proj_gp);
Nm_in_mgp = m_fe.N(m_proj_gp);

for k=1:nN_ele_s
  con_id = nN_m*3 + s_nod(k)*3-2;
  C_s(:,con_id:con_id+2) = eye(3)*Ns_in_sgp(k);
end
for l=1:nN_ele_m
  con_id = m_nod(l)*3-2;
  C_m(:,con_id:con_id+2) = eye(3)*Nm_in_mgp(l);
end

ghat = (m_coo*Nm_in_mgp-s_coo*Ns_in_sgp)'*Ngp + normal_gp'*(m_coo*dNm_in_mgp)*G_m(1:2,:)...
  -normal_gp'*(s_coo*dNs_in_sgp)*G_s(1:2,:)+normal_gp'*C_m-normal_gp'*C_s;