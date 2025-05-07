function Ngp = linear_normal_gp(normal_dir_gp, Nhat)

Ngp = (eye(3)/norm(normal_dir_gp)-normal_dir_gp*normal_dir_gp'*norm(normal_dir_gp)^(-3))*Nhat;