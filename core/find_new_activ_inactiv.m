function [new_activ_loc,new_inact_loc] = find_new_activ_inactiv(weighted_gap, compl_param, z, normals_storage)

n = size(z,1);
z_nu = zeros(n/3,1);
%avg_normals = reshape(normals_storage.averaged_normals',n,1);
%tangents1 = reshape(normals_storage.t1',n,1);
%tangents2 = reshape(normals_storage.t2',n,1);
%{
for i=1:n/3
  id = 3*i-2;
  components = linsolve([avg_normals(id:id+2),tangents1(id:id+2),tangents2(id:id+2)],z(id:id+2));
  z_nu(i)=components(1);
end
%}
for i=1:n/3
  z_nu(i) = normals_storage.averaged_normals(i,:)*z(3*i-2:3*i);
end
penet_cond = z_nu - compl_param * weighted_gap;  % vector of conditions
penet_cond(abs(weighted_gap) < 1e-4) = 1;

new_activ_loc = find(penet_cond > 0);              % indices where condition is >0
new_inact_loc = find(penet_cond <= 0);             % indices where condition is <=0


