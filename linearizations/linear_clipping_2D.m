function clip_Vtilde = linear_clipping_2D(clip_origin, sl_rotP_alg, mast_rotP_alg, rot_s_coo, rot_m_coo)

n_vert = size(clip_origin, 1);
clip_Vtilde = cell(n_vert,1);

D_new_z = zeros(1,size(sl_rotP_alg{1},2));
for k=1:length(sl_rotP_alg)
  D_new_z = D_new_z + sl_rotP_alg{k}(3,:);
end
for l=1:length(mast_rotP_alg)
  D_new_z = D_new_z + mast_rotP_alg{l}(3,:);
end
D_new_z = D_new_z/(length(sl_rotP_alg) + length(mast_rotP_alg));

for v=1:n_vert
  % case 1: slave - origin projected slave node
  if clip_origin(v, 1) == 1
    clip_Vtilde{v} = sl_rotP_alg{clip_origin(v,2)};
    clip_Vtilde{v}(3,:) = D_new_z;
    % case 1: master - origin projected master node
  elseif clip_origin(v, 1) == 2
    clip_Vtilde{v} = mast_rotP_alg{clip_origin(v,2)};
    clip_Vtilde{v}(3,:) = D_new_z;
    % case 2 - line clipping
  elseif clip_origin(v, 1) == 3
    A = rot_s_coo(:,clip_origin(v,2)); B = rot_s_coo(:,clip_origin(v,3));
    C = rot_m_coo(:,clip_origin(v,4)); D = rot_m_coo(:,clip_origin(v,5));
    Ptilde_A = sl_rotP_alg{clip_origin(v,2)};
    Ptilde_B = sl_rotP_alg{clip_origin(v,3)};
    Ptilde_C = mast_rotP_alg{clip_origin(v,4)};
    Ptilde_D = mast_rotP_alg{clip_origin(v,5)};

    clip_Vtilde{v} = linear_inters_vertex(A,B,C,D,Ptilde_A,Ptilde_B,Ptilde_C,Ptilde_D);
    clip_Vtilde{v}(3,:) = D_new_z;
  end
end


