function S = skew_mat(v)
% Returns the skew-symmetric matrix of a 3D vector v

S = [   0   -v(3)  v(2);
       v(3)   0   -v(1);
      -v(2)  v(1)   0 ];
end