function [rot_clip, clip_origin] = clipping_2D(rot_s, rot_m)

rot_s = rot_s';
rot_m = rot_m';
nN_s = size(rot_s,1);
nN_m = size(rot_m,1);

%plot_3d_polygon(rot_m, rot_s);

[X, Y] = polyclip(rot_s(:,1), rot_s(:,2), rot_m(:,1), rot_m(:,2), 1); %clip

if isempty(X) || isempty(Y)
  rot_clip = [];  % Return empty if no intersection is found
  clip_origin = [];
  return;
end

% take the clipped intersection and give it 3rd dim
rot_clip = zeros(size(X{1},1), 3);
clip_origin = zeros(size(X{1},1), 5);
rot_clip(:, 1) = X{1}; rot_clip(:, 2) = Y{1}; rot_clip(:, 3) = rot_s(1,3);


% need origin of vertices for linearization, look Popp 3d Figure A1
% we dont really need to distinguish between case 1 and 2 (slave and master)
tol = 1e-8; % Small tolerance for numerical precision
for i = 1:size(rot_clip,1)
  for j=1:nN_s
    if (norm(rot_clip(i,1:2)-rot_s(j,1:2)) < tol)
      clip_origin(i,1) = 1; % From projected slave
      clip_origin(i,2) = j; % From projected slave
    end
  end
  for j=1:nN_m
    if (norm(rot_clip(i,1:2)-rot_m(j,1:2)) < tol)
      clip_origin(i,1) = 2; % From projected master
      clip_origin(i,2) = j; % From projected master
    end
  end
end

for i=1:size(clip_origin, 1)
  if clip_origin(i,1) == 0
    clip_origin(i,1) = 3; % Indicates intersection between edges

    % Check which segment of projected slave contains the intersection
    for j=1:nN_s
      j_next = mod(j, nN_s) + 1; % Next index with wrap-around
      if check_point_on_segment(rot_s(j,1:2), rot_s(j_next, 1:2), rot_clip(i,1:2))
        clip_origin(i,2) = j;      % Start index of segment
        clip_origin(i,3) = j_next; % End index of segment
        
      end
    end

    % Check which segment of projected master contains the intersection
    for j=1:nN_m
      j_next = mod(j, nN_m) + 1;
      if check_point_on_segment(rot_m(j,1:2), rot_m(j_next, 1:2), rot_clip(i,1:2))
        clip_origin(i,4) = j;      % Start index of segment
        clip_origin(i,5) = j_next; % End index of segment
        
      end
    end
  end
end
rot_clip = rot_clip';



