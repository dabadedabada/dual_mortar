function delt_integr_cells = linear_integr_cell(s_nod, m_nod, delt_s_nod, delt_m_nod, x0, n0, delt_n0, delt_x0, clip, clip_origin, proj_s, proj_m)
% compute linearizations of clip vertices and its center, return derivatives
% of triangle integration cells as cell(ncells,1)
% center is computed by simple averaging, as in Puso, Laursen
% vertices originated from node projections as Popp 3d (A17)
% intersection vertices as Popp 3d (A19) (propably should be different as our clipping aproach is also different)

% s_nod, m_nod are element nodes in anti-clockwise direction
% proj_s, proj_m are the same nodes projected onto aux plane
% clip are vertices in anti-clockwise direction
% clip_origin(:,1) are cases of origin of vertices as in Popp 3d Sec A.5
% if clip_origin(:,1)= 1 or 2, clip_origin(:,2) are indices of rows of
% s_nod, m_nod (or consequently proj_s, proj_m)
% if clip_origin(:,1)= 3, clip_origin(:,2:5) are indices of slave start
% and end and master start and end point in this order in proj_s, proj_m (or consequently the unprojected s_nod, m_nod)


delt_clip = zeros(size(clip));
for i=1:size(clip_origin, 1)
  % case 1 - origin projected slave node
  if clip_origin(i, 1) == 1
    delt_clip(i,:) = linear_project_onto_plane(s_nod(clip_origin(i,2),:), x0, n0, delt_s_nod(clip_origin(i,2),:), delt_x0, delt_n0);
  end
  % case 2 - origin projected master node
  if clip_origin(i, 1) == 2
    delt_clip(i,:) = linear_project_onto_plane(m_nod(clip_origin(i,2),:), x0, n0, delt_m_nod(clip_origin(i,2),:), delt_x0, delt_n0);
  end
  % case 3 - line clipping
  if clip_origin(i, 1) == 3
    s_s = proj_s(clip_origin(i,2),:); s_e = proj_s(clip_origin(i,3),:);   
    m_s = proj_m(clip_origin(i,4),:); m_e = proj_m(clip_origin(i,5),:);
    delt_s_s = linear_project_onto_plane(s_nod(clip_origin(i,2),:), x0, n0, delt_s_nod(clip_origin(i,2),:), delt_x0, delt_n0);
    delt_s_e = linear_project_onto_plane(s_nod(clip_origin(i,3),:), x0, n0, delt_s_nod(clip_origin(i,3),:), delt_x0, delt_n0);
    delt_m_s = linear_project_onto_plane(m_nod(clip_origin(i,4),:), x0, n0, delt_m_nod(clip_origin(i,4),:), delt_x0, delt_n0);
    delt_m_e = linear_project_onto_plane(m_nod(clip_origin(i,5),:), x0, n0, delt_m_nod(clip_origin(i,5),:), delt_x0, delt_n0);
    
    delt_clip(i,:) = linear_inters_vertex(s_s,s_e,m_s,m_e,delt_s_s,delt_s_e,delt_m_s,delt_m_e,n0,delt_n0);
  end
end

% case 4 - center of the clip polygon (by averaging)
delt_clip_centr = mean(delt_clip,1);

% create integration cells same as in Algorithm 1
ncells = size(delt_clip, 1);
delt_integr_cells = cell(ncells,1);
for j=1:ncells
  delt_integr_cells{j} = [delt_clip_centr; delt_clip(j,:); delt_clip(mod(j, ncells)+1,:)];
end





