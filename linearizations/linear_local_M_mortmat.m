function Malg_loc = linear_local_M_mortmat(s_fe, m_fe, w, jalg, Jcell, s_proj_gp, m_proj_gp, G_m, Ae, F)

s_N_in_gp = s_fe.N(s_proj_gp);
m_N_in_gp = m_fe.N(m_proj_gp);
s_nN_ele = length(s_N_in_gp);
m_nN_ele = length(m_N_in_gp);
size_con = size(G_m,2);
m_dN_in_gp = m_fe.dNdxi(m_proj_gp);
dshpf_in_gp = Ae*s_N_in_gp;

Malg_loc = zeros(s_nN_ele,m_nN_ele, size_con);
for k=1:s_nN_ele
  for l=1:m_nN_ele
    Malg_loc(k,l,:) = w*dshpf_in_gp(k)*m_N_in_gp(l)*jalg + ...
      w*Jcell*dshpf_in_gp(k)*m_dN_in_gp(l,:)*G_m(1:2,:)+...
      w*Jcell*m_N_in_gp(l)*F(k,:);
  end
end
