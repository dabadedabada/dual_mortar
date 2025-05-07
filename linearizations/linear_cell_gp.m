function Ghat = linear_cell_gp(V1, V2, V3, fe_cell, gp_id)

gauss_points = fe_cell.QP;
N_in_gp = fe_cell.N(gauss_points(:,gp_id));
Ghat = N_in_gp(1)*V1 + N_in_gp(2)*V2 + N_in_gp(3)*V3;

