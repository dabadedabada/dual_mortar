function V = linear_rot_polyg_vert(s_coo, m_coo, x0, n0, N0, C0, ...
  clip_origin, proj_s_coo, proj_m_coo, rot_s_coo, rot_m_coo, nN_m, R, Ralg, rot_clip)

n_vert = size(clip_origin, 1);
clip_Vtilde = cell(n_vert,1);
for v=1:n_vert
  % case 1 - origin projected slave node
  if clip_origin(v, 1) == 1
    P = linear_project_onto_plane(s_coo(:,clip_origin(v,2)), clip_origin(v,2), x0, n0, N0, C0, 1, nN_m);  
    clip_Vtilde{v} = linear_rotate_points(R, Ralg, proj_s_coo(:,clip_origin(v,2)), P, x0, C0);
  % case 1 - origin projected master node
  elseif clip_origin(v, 1) == 2
    P = linear_project_onto_plane(m_coo(:,clip_origin(v,2)), clip_origin(v,2), x0, n0, N0, C0, 2, nN_m);  
    clip_Vtilde{v} = linear_rotate_points(R, Ralg, proj_m_coo(:,clip_origin(v,2)), P, x0, C0);
  % case 2 - line clipping
  elseif clip_origin(v, 1) == 3
    A = rot_s_coo(:,clip_origin(v,2)); B = rot_s_coo(:,clip_origin(v,3));   
    C = rot_m_coo(:,clip_origin(v,4)); D = rot_m_coo(:,clip_origin(v,5));
    P_A = linear_project_onto_plane(s_coo(:,clip_origin(v,2)), clip_origin(v,2), x0, n0, N0, C0, 1, nN_m);
    P_B = linear_project_onto_plane(s_coo(:,clip_origin(v,3)), clip_origin(v,3), x0, n0, N0, C0, 1, nN_m);
    P_C = linear_project_onto_plane(s_coo(:,clip_origin(v,4)), clip_origin(v,4), x0, n0, N0, C0, 2, nN_m);
    P_D = linear_project_onto_plane(s_coo(:,clip_origin(v,5)), clip_origin(v,5), x0, n0, N0, C0, 2, nN_m);
    Ptilde_A = linear_rotate_points(R, Ralg, proj_s_coo(:,clip_origin(v,2)), P_A, x0, C0);
    Ptilde_B = linear_rotate_points(R, Ralg, proj_s_coo(:,clip_origin(v,3)), P_B, x0, C0);
    Ptilde_C = linear_rotate_points(R, Ralg, proj_m_coo(:,clip_origin(v,4)), P_C, x0, C0);
    Ptilde_D = linear_rotate_points(R, Ralg, proj_m_coo(:,clip_origin(v,5)), P_D, x0, C0);
    
    clip_Vtilde{v} = linear_inters_vertex(A,B,C,D,Ptilde_A,Ptilde_B,Ptilde_C,Ptilde_D);
  end
end

% Centroid of the clip polygon (by averaging)
V0tilde = zeros(size(clip_Vtilde{1}));
for v=1:n_vert
  V0tilde = V0tilde + clip_Vtilde{v};
end
V0tilde = V0tilde/n_vert;
rot_centr = mean(rot_clip, 2);

% create integration cells
% matrix has 3 matrices of vertices on cell
V = cell(n_vert,3);
for v=1:n_vert
  V{v,1} = linear_reverse_rotate_points(R, Ralg, x0, C0, rot_centr, V0tilde);
  V{v,2} = linear_reverse_rotate_points(R, Ralg, x0, C0, rot_clip(:,v), clip_Vtilde{v});
  V{v,3} = linear_reverse_rotate_points(R, Ralg, x0, C0, rot_clip(:,mod(v, n_vert)+1), clip_Vtilde{mod(v, n_vert)+1});
end



