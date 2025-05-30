%% test Popp Gitterle Gee Wall (A3) for quad4
nj  = 4; % global indices
nk  = 4; % local  indices
nxi = 2;

nd  = 3;
ne  = 1;

% think of a table e*k, same as .nod in cont_face or mesh
%   k  1  2  3  4
% e    -----------
% 1 |  
% 2 |
% 3 |
% 4 |

dN = sym('dN_%d%d%d%d',[ne, nk, nj, nxi]   ,'real');    % e k j [xi,eta]
x  = sym('x_%d%d%d' ,  [ne, nk, nd]        ,'real');    % e k d
Dx = sym('Dx_%d%d%d',  [ne, nk, nd]        ,'real');    % e k d
% each reshape(Dx(e,:,:), )

n  = sym(zeros(        [ne, nj, nd]      ));  % e j d  
Dn = sym(zeros(        [ne, nj, nd]      ));  % e j d  
M_ = sym(zeros(        [ne, nj, nd,nk,nd]));  % e j d k d,   Dn = M(e,j,d,:,:)*Dx


%{

for e = 1:ne
  for j = 1:nj
    tmpl = sym(zeros(3,1));   tmpDl = sym(zeros(3,1)); 
    tmpr = sym(zeros(3,1));   tmpDr = sym(zeros(3,1));
    for k = 1:nk
      tmpl  = tmpl  + dN(k,1,j)* x(k,:)';   tmpDl = tmpDl + dN(k,1,j)*Dx(k,:)';
      tmpr  = tmpr  + dN(k,2,j)* x(k,:)';   tmpDr = tmpDr + dN(k,2,j)*Dx(k,:)';
    end
    n( e,j,:) = n( e,j,:) + reshape( cross( tmpl , tmpr)                      , [1 1 nd]);
    Dn(e,j,:) = Dn(e,j,:) + reshape( cross( tmpDl, tmpr) + cross( tmpl, tmpDr), [1 1 nd]);
  end
end
subs_for   = {Dx(1,1),Dx(1,2),Dx(1,3), Dx(2,1),Dx(2,2),Dx(2,3), Dx(3,1),Dx(3,2),Dx(3,3), Dx(4,1),Dx(4,2),Dx(4,3)};
subs_what_ = {      0,      0,      0,       0,      0,      0,       0,      0,      0,       0,      0,      0};
for e = 1:ne
  for j = 1:nj
    for d1 = 1:nd
      for k = 1:nk
        for d2 = 1:nd
          subs_what = subs_what_;   subs_what{nd*(k-1)+d2} = 1;
          M_(1,j,d1,k,d2) = simplify(subs( Dn(e,j,d1),subs_for,subs_what));
        end
      end
    end
  end
end
%}
%A = sym('A','real');   B = sym('B','real');   C = sym('C','real');   ABC = [A; B; C];   D = sym('D','real');   E = sym('E','real');   F = sym('F','real');   DEF = [A; B; C]; 
%a = sym('a','real');   b = sym('b','real');   c = sym('c','real');   abc = [a; b; c];   d = sym('d','real');   e = sym('e','real');   f = sym('f','real');   def = [d; e; f];
%disp([0 -C  B;  C 0 -A; -B  A 0]*abc - cross(ABC,abc));
%disp([0  C -B; -C 0  A;  B -A 0]*abc - cross(abc,ABC));