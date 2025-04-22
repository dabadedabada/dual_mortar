function proj_gp = newton_it_for_projgp(coo, fe, n0, x_gp)
% returns projections of Gauss points x_g (global coordinates) onto the
% isoparametric element (slave or master) along the slave element normal n0
% and coeficient alpha for scaling n0

n_gp = size(x_gp, 2);
proj_gp = zeros(3, n_gp);
tol = 1e-6;
max_iter = 10;

for gp=1:n_gp
  % starting point n0+x_g, alpha*n0+x_g is the projection from aux plane
  % to the element in current configuration
  xi_1 = [0;0;0]; 
  xi_2 = [1;1;1]; 
  iter = 0;

  while (norm(xi_2-xi_1) > tol) && (iter <= max_iter)
    xi_1 = xi_2;
    J = newton_it_J(coo, n0, fe, xi_1(1:2));
    F = newton_it_F(coo, n0, x_gp(:,gp), fe, xi_1(1:2), xi_1(3));
    xi_2 = xi_1 - J^-1 *F;
    iter = iter +1;
  end
  if iter == max_iter
      warning('Newton method did not converge within %d iterations on the %d . point', max_iter, gp);
  end
  proj_gp(1,gp) = xi_2(1);
  proj_gp(2,gp) = xi_2(2);
  proj_gp(3,gp) = xi_2(3); 
end