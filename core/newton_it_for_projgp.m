function proj_gp = newton_it_for_projgp(element_nod, element_fe, n0, x_g)
% returns projections of Gauss points x_g (global coordinates) onto the
% isoparametric element (slave or master) along the slave element normal n0
% and coeficient alpha for scaling n0
  n_gp = size(x_g, 1);
  proj_gp = zeros(n_gp, 3);
  fe = element_fe;
  dN = fe.dNdxi([0,0,0]'); % derivatives are constant
  tol = 1e-6;
  max_iter = 15;

  % assemble Jacobian matrix
  J = zeros(3);
  for i=1:3
    J(i,1) = dN(:,1)'*element_nod(:,i);
    J(i,2) = dN(:,2)'*element_nod(:,i);
    J(i,3) = -n0(i);
  end

  for i=1:n_gp
    % starting point n0+x_g, alpha*n0+x_g is the projection from aux plane
    % to the element in current configuration
    xi = n0(1)+x_g(i,1); 
    eta = n0(2)+x_g(i,2); 
    alpha = 1;
    F = ones(3,1);
    iter = 0;

    while (norm(F) > tol) && (iter < max_iter)
      % assemble the equations
      N = fe.N([xi, eta]');
      for j=1:3
        F(j) = N'*element_nod(:,j) - alpha*n0(j)-x_g(i, j);
      end
      
      c = -J^(-1)*F;
      xi = xi + c(1); eta = eta + c(2); alpha = alpha + c(3);
      iter = iter+1;
    end

    if norm(F) > tol
        warning('Newton method did not converge within %d iterations on the %d . point', max_iter, i);
    end
    proj_gp(i, 1) = xi;
    proj_gp(i, 2) = eta;
    proj_gp(i, 3) = alpha;
    
  end

end