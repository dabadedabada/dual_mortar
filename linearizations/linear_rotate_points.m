function Ptilde = linear_rotate_points(R, Ralg, proj_coo, P, x0, C0)
% proj_coo is 3x1 coo of projected node
% P is 3x3(nNm + nNs) representing dir der of proj_coo
% x0 is 3x1 centroid coo
% C0 is 3x3(nNm + nNs) representing dir der of x0
% R is rotation matrix
% Ralg is 3D matrix representing dir der of R

Ptilde = squeeze(pagemtimes(Ralg,proj_coo-x0))+R*(P-C0) + C0;