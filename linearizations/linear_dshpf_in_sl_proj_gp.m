function F = linear_dshpf_in_sl_proj_gp(Ae, Aalg, fe, proj_gp, G, nN_m)

N_in_proj_gp = fe.N(proj_gp);
dN_in_proj_gp = fe.dNdxi(proj_gp);
ddshpf_in_proj_gp = Ae*dN_in_proj_gp;

Fs = [zeros(size(N_in_proj_gp,1),3*nN_m),squeeze(pagemtimes(Aalg,N_in_proj_gp))];

F = Fs + ddshpf_in_proj_gp*G(1:2,:);
