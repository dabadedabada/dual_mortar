%% test Popp Gitterle Gee Wall (A3) for tria3
nk  = 4;
nxi = 2;
nj  = nk;
nd  = 3;
ne  = 1;

dN = sym('dN_%d%d%d',[nk,nxi,nj],'real'); % k xi j
x  = sym( 'x_%d%d'  ,[nk,nd]    ,'real'); % k d
Dx = sym('Dx_%d%d'  ,[nk,nd]    ,'real'); % k d

n  = sym(zeros(      [ne,nj,nd]      ));  % e j d  
Dn = sym(zeros(      [ne,nj,nd]      ));  % e j d  
M_ = sym(zeros(      [ne,nj,nd,nk,nd]));  % e j d  k d,   Dn = M(e,j,d,:,:)*Dx
M  = M_;
Nhat_je = M;

for e = 1:ne
  for j = 1:nj
    tmpl = sym(zeros(3,1));   tmpDl = sym(zeros(3,1)); 
    tmpr = sym(zeros(3,1));   tmpDr = sym(zeros(3,1));
    for k = 1:nk
      tmpl  = tmpl  + dN(k,1,j)* x(k,:)';   tmpDl = tmpDl + dN(k,1,j)*Dx(k,:)';
      tmpr  = tmpr  + dN(k,2,j)* x(k,:)';   tmpDr = tmpDr + dN(k,2,j)*Dx(k,:)';
    end
    n( e,j,:) = n( e,j,:) + reshape( cross( tmpl , tmpr)                      , [1 1 nd]);
    Dn(e,j,:) = Dn(e,j,:) + reshape( cross( tmpDl, tmpr) + cross( tmpl, tmpDr), [1 1 nd]);
  end
end
subs_for = Dx;
switch nk
  case 3
    subs_what_ = {      0,      0,      0;       0,      0,      0;       0,      0,      0};
  case 4
    subs_what_ = {      0,      0,      0;       0,      0,      0;       0,      0,      0;       0,      0,      0};
end
for e = 1:ne
  for j = 1:nj
    for d1 = 1:nd
      for k = 1:nk
        for d2 = 1:nd
          subs_what = subs_what_;   subs_what{k,d2} = 1;
          M_(e,j,d1,k,d2) = simplify(subs( Dn(e,j,d1),subs_for,subs_what));
        end
      end
    end
  end
end

for e=1:ne
  for j=1:nj
    tmpl = dN(:,1,j)'*x; tmpr = dN(:,2,j)'*x;
    for k=1:nk
      Nhat_je(e,j,:,k,:) = skew_mat(tmpl)'*(eye(3)*dN(k,2,j)) + skew_mat(tmpr)*(eye(3)*dN(k,1,j));
    end
  end
end


for e = 1:ne
  for j = 1:nj
    for k = 1:nk
      tmp = 1:nk; tmp(k) = [];
      switch nk
        case 3
          M(e,j,1,k,2) =  (dN(k,1,j)*dN(tmp(1),2,j)-dN(k,2,j)*dN(tmp(1),1,j))*x(tmp(1),3) + (dN(k,1,j)*dN(tmp(2),2,j)-dN(k,2,j)*dN(tmp(2),1,j))*x(tmp(2),3);
          M(e,j,1,k,3) = -(dN(k,1,j)*dN(tmp(1),2,j)-dN(k,2,j)*dN(tmp(1),1,j))*x(tmp(1),2) - (dN(k,1,j)*dN(tmp(2),2,j)-dN(k,2,j)*dN(tmp(2),1,j))*x(tmp(2),2);
          M(e,j,2,k,3) =  (dN(k,1,j)*dN(tmp(1),2,j)-dN(k,2,j)*dN(tmp(1),1,j))*x(tmp(1),1) + (dN(k,1,j)*dN(tmp(2),2,j)-dN(k,2,j)*dN(tmp(2),1,j))*x(tmp(2),1);
        case 4
          M(e,j,1,k,2) =  (dN(k,1,j)*dN(tmp(1),2,j)-dN(k,2,j)*dN(tmp(1),1,j))*x(tmp(1),3) + (dN(k,1,j)*dN(tmp(2),2,j)-dN(k,2,j)*dN(tmp(2),1,j))*x(tmp(2),3) + (dN(k,1,j)*dN(tmp(3),2,j)-dN(k,2,j)*dN(tmp(3),1,j))*x(tmp(3),3);
          M(e,j,1,k,3) = -(dN(k,1,j)*dN(tmp(1),2,j)-dN(k,2,j)*dN(tmp(1),1,j))*x(tmp(1),2) - (dN(k,1,j)*dN(tmp(2),2,j)-dN(k,2,j)*dN(tmp(2),1,j))*x(tmp(2),2) - (dN(k,1,j)*dN(tmp(3),2,j)-dN(k,2,j)*dN(tmp(3),1,j))*x(tmp(3),2);
          M(e,j,2,k,3) =  (dN(k,1,j)*dN(tmp(1),2,j)-dN(k,2,j)*dN(tmp(1),1,j))*x(tmp(1),1) + (dN(k,1,j)*dN(tmp(2),2,j)-dN(k,2,j)*dN(tmp(2),1,j))*x(tmp(2),1) + (dN(k,1,j)*dN(tmp(3),2,j)-dN(k,2,j)*dN(tmp(3),1,j))*x(tmp(3),1);
      end
      M(e,j,2,k,1) = -M(e,j,1,k,2);
      M(e,j,3,k,1) = -M(e,j,1,k,3);
      M(e,j,3,k,2) = -M(e,j,2,k,3);
    end
  end
end

% --------------------------------------- %
% ----- Numerical testing --------------- %
mesh = mesh_generator([3,3,2], [1,1,1], 'hexa8');
cont_face = init_cont_face(mesh, 5);
normals_storage = get_normals_centroids(cont_face);

dN_num = zeros(nk,nxi,nj);
for j=1:nj
  dN_num(:,:,j) = cont_face.fe.dNdxi(cont_face.fe.xi(j,:)');
end

x_num = cont_face.coo(cont_face.nod(1,:),:);

Dx_num = zeros(nk*nd,1);
[subs_vars, subs_vals] = build_substitution_lists(x, dN, x_num, dN_num, Dx, reshape(Dx_num,nd,nk)');
n0 = reshape(double(subs(n, subs_vars, subs_vals)),nk,nd)'; 
n0 = reshape(n0, nk*nd,1);

err = 10^-6;
epsilon = 10^-8;
for kd = 1:nk*nd
  Dx_num = zeros(nk*nd,1);
  Dx_num(kd) = 1; 

  x_num = x_num + epsilon*reshape(Dx_num, nd, nk)';

  [subs_vars, subs_vals] = build_substitution_lists(x, dN, x_num, dN_num, Dx, reshape(Dx_num,nd,nk)');
  
  M_num = double(subs(M, subs_vars, subs_vals));
  M_num = reshape(permute(M_num,[1 3 2 5 4]),nj*nd,nk*nd);    % e j.d  k.d
  Dn_num_alg = M_num*Dx_num;
  
  n_curr = reshape(double(subs(n, subs_vars, subs_vals)),nk,nd)';
  n_curr = reshape(n_curr, nk*nd,1);  
  Dn_num_lim = (n_curr-n0)/epsilon; 
  if (norm(Dn_num_lim-Dn_num_alg) < err)
     fprintf('For epsilon = %g: Dn  DoF: %d         - Success - with difference of dir derivatives %g \n', epsilon, kd, norm(Dn_num_lim-Dn_num_alg));
  else
     fprintf('For epsilon = %g: Dn  DoF: %d         - Fail    - with difference of dir derivatives %g \n', epsilon, kd, norm(Dn_num_lim-Dn_num_alg));
  end
end


