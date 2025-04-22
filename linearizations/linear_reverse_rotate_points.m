function V = linear_reverse_rotate_points(R, Ralg, x0, C0, rot_vert_coo, Vtilde)
% rot_vert_coo is 3x1 coo of rotated polygon vertice
% Vtilde is 3x3(nNm + nNs) representing dir der of rot_vert_coo
% x0 is 3x1 centroid coo
% C0 is 3x3(nNm + nNs) representing dir der of x0
% R is rotation matrix
% Ralg is 3D matrix representing dir der of R

V = squeeze(pagemtimes(permute(Ralg, [2, 1, 3]),rot_vert_coo-x0))+R'*(Vtilde-C0) + C0;