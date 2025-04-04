function delt_inter_vert = linear_inters_vertex(s_s, s_e, m_s, m_e, delt_s_s, delt_s_e, delt_m_s, delt_m_e, n0, delt_n0)
% returns derivative of vertex originated as intersection of line sections
% of slave and master start and end nodes, Popp 3d (A19)

term2 = cross(s_s-m_s, m_e-m_s);
delt_term2 = cross(delt_s_s-delt_m_s, m_e-m_s) + cross(s_s-m_s, delt_m_e-delt_m_s);
term3 = cross(s_e-s_s, m_e-m_s);
delt_term3 = cross(delt_s_e-delt_s_s, m_e-m_s) + cross(s_e-s_s, delt_m_e-delt_m_s);
numerator = dot(term2, n0);
delt_numerator = dot(delt_term2, n0) + dot(term2, delt_n0);
denominator = dot(term3, n0);
delt_denominator = dot(delt_term3, n0) + dot(term3, delt_n0);
frac = numerator/denominator;
delt_frac = (delt_numerator*denominator - numerator*delt_denominator)/(denominator)^2;
term4 = s_e-s_s;
delt_term4 = delt_s_e-delt_s_s;
delt_inter_vert = delt_s_s - (delt_frac*term4 + frac*delt_term4);