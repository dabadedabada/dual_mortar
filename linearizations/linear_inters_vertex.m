function T = linear_inters_vertex(s_s, s_e, m_s, m_e, T_s_s, T_s_e, T_m_s, T_m_e, n0, N0)
% returns derivative of vertex originated as intersection of line sections
% of slave and master start and end nodes, Popp 3d (A19)

sl_edge = s_e-s_s;
m_edge = m_e-m_s;
weird_term = s_s-m_s;

f = dot(cross(weird_term, m_edge), n0);
g = dot(cross(sl_edge, m_edge), n0);

delt_sl_edge = T_s_e-T_s_s;
delt_m_edge = T_m_e-T_m_s;
delt_weird_term = T_s_s-T_m_s;

delt_f = cross(weird_term, m_edge)*N0 + n0*skew_mat(m_edge)'*delt_weird_term + n0*skew_mat(weird_term)*delt_m_edge;
delt_g = cross(sl_edge, m_edge)*N0 + n0*skew_mat(m_edge)'*delt_sl_edge + n0*skew_mat(sl_edge)*delt_m_edge;

T = T_s_s - 1/(g*g)*sl_edge'*(delt_f*g - f*delt_g) - f/g*delt_sl_edge;