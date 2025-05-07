function [Aext,Bext] = extend_points_to_r(A,B,r)

A_2D = [A(1);A(3)];
B_2D = [B(1);B(3)];

center = (A_2D+B_2D)/2;

vA = r*(A_2D-center)/norm(A_2D-center);
vB = r*(B_2D-center)/norm(B_2D-center);

A_2D_ext = center+vA;
B_2D_ext = center+vB;

Aext = [A_2D_ext(1),A(2),A_2D_ext(2)];
Bext = [B_2D_ext(1),B(2),B_2D_ext(2)];