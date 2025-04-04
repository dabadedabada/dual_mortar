function delt_disc_gapf = linear_disc_gapf(s_nod, m_nod, delt_s_nod, delt_m_nod, s_fe, m_fe, delt_avg_normals, avg_normals, s_delt_proj_gp, m_delt_proj_gp, s_proj_gp, m_proj_gp)

%Ns_in_sgp = s_fe.N(s_proj_gp');
%Nm_in_mgp = m_fe.N(m_proj_gp');
%disc_gapf = dot(Ns_in_sgp'*s.normals, Nm_in_mgp'*m.nod-Ns_in_sgp'*s.nod, 2); % (A9) Popp 3D 

n_gp = size(s_delt_proj_gp,1);
delt_disc_gapf = zeros(n_gp, 1);

for i=1:n_gp
  Ns_in_sgp = s_fe.N(s_proj_gp(i,:)');
  Nm_in_mgp = m_fe.N(m_proj_gp(i,:)');
  dNs_in_sgp = s_fe.dNdxi(s_proj_gp(i,:)');
  dNm_in_mgp = m_fe.dNdxi(m_proj_gp(i,:)');
  term1 = dot(s_delt_proj_gp(i,:)*dNs_in_sgp'*avg_normals+Ns_in_sgp'*delt_avg_normals, Nm_in_mgp'*m_nod-Ns_in_sgp'*s_nod);
  term2 = dot(Ns_in_sgp'*avg_normals, m_delt_proj_gp(i,:)*dNm_in_mgp'*m_nod+Nm_in_mgp'*delt_m_nod);
  term3 = dot(Ns_in_sgp'*avg_normals, s_delt_proj_gp(i,:)*dNs_in_sgp'*s_nod+Ns_in_sgp'*delt_s_nod);
  delt_disc_gapf(i)=term1+term2-term3; % Popp 3D (A10)
end