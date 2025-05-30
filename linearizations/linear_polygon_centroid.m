function V0tilde = linear_polygon_centroid(clip_Vtilde)

n_vert = length(clip_Vtilde);
V0tilde = zeros(size(clip_Vtilde{1}));
for v=1:n_vert
  V0tilde = V0tilde + clip_Vtilde{v};
end
V0tilde = V0tilde/n_vert;

