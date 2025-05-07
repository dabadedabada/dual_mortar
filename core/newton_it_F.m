function F = newton_it_F(coo, n0, x_gp, fe, xi, alpha)
N = fe.N(xi);
F = coo*N - alpha*n0 - x_gp;