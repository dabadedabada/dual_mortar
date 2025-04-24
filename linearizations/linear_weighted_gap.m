function gtilde = linear_weighted_gap(disc_gapf, w, jalg, Jcell, Ae, s_proj_gp, s_fe, F, ghat)

Ns_in_sgp = s_fe.N(s_proj_gp);
dshpf_in_sgp = Ae*Ns_in_sgp;

gtilde = w*(disc_gapf*Jcell*F + dshpf_in_sgp*Jcell*ghat + dshpf_in_sgp*disc_gapf*jalg);