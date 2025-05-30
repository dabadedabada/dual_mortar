function Vtilde = linear_inters_vertex(A, B, C, D, Ptilde_A, Ptilde_B, Ptilde_C, Ptilde_D)

if length(A) == 3
  P = [0,1,0;-1,0,0;0,0,0];
elseif length(A) == 2
  P = [0,1;-1,0];
end

Vtilde = Ptilde_D + (B-D)'*P*(A-B)*((C-D)'*P*(A-B))^(-1)*(Ptilde_C-Ptilde_D)+...
  (C-D)*(A-B)'*P'*((C-D)'*P*(A-B))^(-2)*((C-D)*(A-B)'*P'*Ptilde_B +...
  (D-B)*(A-B)'*P'*Ptilde_C + (B-C)*(A-B)'*P'*Ptilde_D +...
  (B*(D-C)'+C*(B-D)'+D*(C-B)')*P*(Ptilde_A-Ptilde_B));