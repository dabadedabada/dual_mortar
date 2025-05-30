%% test Popp Gitterle Gee Wall (A3) for quad4

addpath("core");
addpath("core_matfem");
addpath(genpath("linearizations"));

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
    n( e,j,:) = n( e,j,:) + reshape( cross( tmpr , tmpl)                      , [1 1 nd]);
    Dn(e,j,:) = Dn(e,j,:) + reshape( cross( tmpDr, tmpl) + cross( tmpr, tmpDl), [1 1 nd]);
  end
end
subs_for   = {Dx(1,1),Dx(1,2),Dx(1,3); Dx(2,1),Dx(2,2),Dx(2,3); Dx(3,1),Dx(3,2),Dx(3,3); Dx(4,1),Dx(4,2),Dx(4,3)};
subs_what_ = {      0,      0,      0;       0,      0,      0;       0,      0,      0;       0,      0,      0};
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

switch_to_interlaced_notation = true;

if switch_to_interlaced_notation
  for e=1:ne
    for j=1:nj
      tmpl = dN(:,1,j)'*x; tmpr = dN(:,2,j)'*x;
      for k=1:nk
        Nhat_je(e,j,:,k,:) = skew_mat(tmpl)'*(eye(3)*dN(k,2,j)) + skew_mat(tmpr)*(eye(3)*dN(k,1,j));
      end
    end
  end
end

for e = 1:ne
  for j = 1:nj
    for k = 1:nk
      tmp = 1:4; tmp(k) = [];
      M(e,j,1,k,2) =  (dN(k,1,j)*dN(tmp(1),2,j)-dN(k,2,j)*dN(tmp(1),1,j))*x(tmp(1),3) + (dN(k,1,j)*dN(tmp(2),2,j)-dN(k,2,j)*dN(tmp(2),1,j))*x(tmp(2),3) + (dN(k,1,j)*dN(tmp(3),2,j)-dN(k,2,j)*dN(tmp(3),1,j))*x(tmp(3),3);   M(e,j,2,k,1) = -M(e,j,1,k,2);
      M(e,j,1,k,3) = -(dN(k,1,j)*dN(tmp(1),2,j)-dN(k,2,j)*dN(tmp(1),1,j))*x(tmp(1),2) - (dN(k,1,j)*dN(tmp(2),2,j)-dN(k,2,j)*dN(tmp(2),1,j))*x(tmp(2),2) - (dN(k,1,j)*dN(tmp(3),2,j)-dN(k,2,j)*dN(tmp(3),1,j))*x(tmp(3),2);   M(e,j,3,k,1) = -M(e,j,1,k,3);
      M(e,j,2,k,3) =  (dN(k,1,j)*dN(tmp(1),2,j)-dN(k,2,j)*dN(tmp(1),1,j))*x(tmp(1),1) + (dN(k,1,j)*dN(tmp(2),2,j)-dN(k,2,j)*dN(tmp(2),1,j))*x(tmp(2),1) + (dN(k,1,j)*dN(tmp(3),2,j)-dN(k,2,j)*dN(tmp(3),1,j))*x(tmp(3),1);   M(e,j,3,k,2) = -M(e,j,2,k,3);
    end
  end
end

% --------------------------------------- %
% ----- Numerical testing --------------- %
mesh = mesh_generator([3,3,2], [1,1,1],"hexa8");
cont_face = init_cont_face(mesh, 5);
normals_storage = get_normals_centroids(cont_face);

x_num = cont_face.coo;
dN_num = zeros(nk,nxi,nj);
for j=1:nj
  dN_num(:,:,j) = cont_face.fe.dNdxi(cont_face.fe.xi(j,:)');
end

subs_vars = {};
subs_vals = {};

% x substitution
for k = 1:nk
  for d = 1:nd
    subs_vars{end+1} = x(k,d);
    subs_vals{end+1} = x_num(k,d);
  end
end

% dN substitution
for k = 1:nk
  for xi = 1:nxi
    for j = 1:nj
      subs_vars{end+1} = dN(k,xi,j);
      subs_vals{end+1} = dN_num(k,xi,j);
    end
  end
end

Mnum = double(subs(M, subs_vars, subs_vals));
M_num = double(subs(M_, subs_vars, subs_vals));

nnum = double(subs(n, subs_vars, subs_vals));

