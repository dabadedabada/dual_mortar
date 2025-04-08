function T_cell = linear_integr_cell(s_coo, m_coo, x0, n0, N0, C0, clip_origin, proj_s, proj_m, nN_m)
% compute linearizations of clip vertices and its center, return matrices
% for delta x^cell_v = T_cell{cell,v} * delta d

% center is computed by simple averaging, as in Puso, Laursen
% vertices originated from node projections as Popp 3d (A17)
% intersection vertices as Popp 3d (A19) (propably should be different as our clipping aproach is also different)

% s_coo, m_coo are element node coordinates in anti-clockwise direction
% proj_s, proj_m are the same nodes projected onto aux plane
% clip_origin(:,1) are cases of origin of vertices as in Popp 3d Sec A.5
% if clip_origin(:,1)= 1 or 2, clip_origin(:,2) are indices of rows of
% s_coo, m_coo (or consequently proj_s, proj_m)
% if clip_origin(:,1)= 3, clip_origin(:,2:5) are indices of slave start
% and end and master start and end point in this order in proj_s, proj_m
ncells = size(clip_origin, 1);
clip_T = cell(ncells,1);
for i=1:ncells
  % case 1 - origin projected slave node
  if clip_origin(i, 1) == 1
    clip_T{i} = linear_project_onto_plane(s_coo(clip_origin(i,2),:), clip_origin(i,2), x0, n0, N0, C0, 1, nN_m);
  end
  % case 2 - origin projected master node
  if clip_origin(i, 1) == 2
    clip_T{i} = linear_project_onto_plane(m_coo(clip_origin(i,2),:), clip_origin(i,2), x0, n0, N0, C0, 2, nN_m);
  end
  % case 3 - line clipping
  if clip_origin(i, 1) == 3
    s_s = proj_s(clip_origin(i,2),:); s_e = proj_s(clip_origin(i,3),:);   
    m_s = proj_m(clip_origin(i,4),:); m_e = proj_m(clip_origin(i,5),:);
    T_s_s = linear_project_onto_plane(s_coo(clip_origin(i,2),:), clip_origin(i,2), x0, n0, N0, C0, 1, nN_m);
    T_s_e = linear_project_onto_plane(s_coo(clip_origin(i,3),:), clip_origin(i,2), x0, n0, N0, C0, 1, nN_m);
    T_m_s = linear_project_onto_plane(m_coo(clip_origin(i,4),:), clip_origin(i,2), x0, n0, N0, C0, 2, nN_m);
    T_m_e = linear_project_onto_plane(m_coo(clip_origin(i,5),:), clip_origin(i,2), x0, n0, N0, C0, 2, nN_m);
    
    clip_T{i} = linear_inters_vertex(s_s, s_e, m_s, m_e,...
      T_s_s,T_s_e,T_m_s,T_m_e,n0,N0);
  end
end

% case 4 - center of the clip polygon (by averaging)
T_centr = zeros(size(clip_T{1}));
for i=1:ncells
  T_centr = T_centr + clip_T{i};
end
T_centr = T_centr/ncells;

% create integration cells same as in Algorithm 1, now each row of cell
% matrix has 3 matrices of vertices on cell
T_cell = cell(ncells,3);
for j=1:ncells
  T_cell{j,1} = T_centr;
  T_cell{j,2} = clip_T{j};
  T_cell{j,3} = clip_T{mod(j, ncells)+1};
end





