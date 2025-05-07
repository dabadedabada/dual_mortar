function J = newton_it_J(coo, n0, fe, xi)
dN = fe.dNdxi(xi);
J = [coo*dN(:,1), coo*dN(:,2), -n0];
