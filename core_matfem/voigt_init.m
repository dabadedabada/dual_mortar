function V = voigt_init(dim, version)

%% map_vA4_A4
  function A = Voigt_map_vA2_A2(vA)
    sA = size(vA);
    A = zeros(V.dim,V.dim,sA(2));
    if isa(vA,'sym'),   A = sym(A);   end
    for i = 1:sA(1)
      A(V.o2_mapsequence(1,i),V.o2_mapsequence(2,i),:) = vA(i,:);
      A(V.o2_mapsequence(2,i),V.o2_mapsequence(1,i),:) = vA(i,:);
    end
  end

  function vA = Voigt_map_A2_vA2(A)
    vA = zeros(V.o2_voidim,size(A,3));
    if isa(A,'sym'),    vA = sym(vA);   end
    for i = 1:V.o2_voidim
      vA(i,:) = A(V.o2_mapsequence(1,i),V.o2_mapsequence(2,i),:);
    end
  end

  function vC = Voigt_ddot_vA2_vB2(vA,vB)
    vC = V.o2_multiplicity*(vA.*vB);
  end

%% map_vA4m_A4
  function A = Voigt_map_vA4m_A4(vA)
    A = zeros(V.dim,V.dim,V.dim,V.dim,size(vA,3));
    if isa(vA,'sym'),   A = sym(A);   end
    for k_ = 1:V.o4_vdim
      % ijkl = ijlk = jikl = jilk = klij = klji = lkij = lkji
      tmp = vA(V.o4_mapsequence2(1,k_),V.o4_mapsequence2(2,k_),:);
      A(V.o4_mapsequence4(1,k_),V.o4_mapsequence4(2,k_),V.o4_mapsequence4(3,k_),V.o4_mapsequence4(4,k_),:) = tmp;
      A(V.o4_mapsequence4(1,k_),V.o4_mapsequence4(2,k_),V.o4_mapsequence4(4,k_),V.o4_mapsequence4(3,k_),:) = tmp;
      A(V.o4_mapsequence4(2,k_),V.o4_mapsequence4(1,k_),V.o4_mapsequence4(3,k_),V.o4_mapsequence4(4,k_),:) = tmp;
      A(V.o4_mapsequence4(2,k_),V.o4_mapsequence4(1,k_),V.o4_mapsequence4(4,k_),V.o4_mapsequence4(3,k_),:) = tmp;
      A(V.o4_mapsequence4(3,k_),V.o4_mapsequence4(4,k_),V.o4_mapsequence4(1,k_),V.o4_mapsequence4(2,k_),:) = tmp;
      A(V.o4_mapsequence4(3,k_),V.o4_mapsequence4(4,k_),V.o4_mapsequence4(2,k_),V.o4_mapsequence4(1,k_),:) = tmp;
      A(V.o4_mapsequence4(4,k_),V.o4_mapsequence4(3,k_),V.o4_mapsequence4(1,k_),V.o4_mapsequence4(2,k_),:) = tmp;
      A(V.o4_mapsequence4(4,k_),V.o4_mapsequence4(3,k_),V.o4_mapsequence4(2,k_),V.o4_mapsequence4(1,k_),:) = tmp;
    end
  end

  function vAmat = Voigt_map_A4_vA4m(A)
    vAmat = zeros(V.o2_voidim,V.o2_voidim,size(A,5));
    if isa(A,'sym'),  vAmat = sym(vAmat);   end
    for k_ = 1:V.o4_vdim
      tmp = A(V.o4_mapsequence4(1,k_),V.o4_mapsequence4(2,k_),V.o4_mapsequence4(3,k_),V.o4_mapsequence4(4,k_),:);
      vAmat(V.o4_mapsequence2(1,k_),V.o4_mapsequence2(2,k_),:) = tmp;
      vAmat(V.o4_mapsequence2(2,k_),V.o4_mapsequence2(1,k_),:) = tmp;
    end
  end

  function vAvec = Voigt_map_A4_vA4(A)
    vAvec = zeros(V.o4_vdim,size(A,5));
    if isa(A,'sym'),  vAvec = sym(vAvec);   end
    for i = 1:V.o4_vdim
      vAvec(i,:) = A(V.o4_mapsequence4(1,i),V.o4_mapsequence4(2,i),V.o4_mapsequence4(3,i),V.o4_mapsequence4(4,i),:);
    end
  end

  function A = Voigt_map_vA4_A4(vAvec)
    A = zeros(V.dim,V.dim,V.dim,V.dim,size(vAvec,2));
    if isa(vAvec,'sym'),  A = sym(A);   end
    for k_ = 1:V.o4_vdim
      % ijkl = ijlk = jikl = jilk = klij = klji = lkij = lkji
      tmp = vAvec(k_,:);
      A(V.o4_mapsequence4(1,k_),V.o4_mapsequence4(2,k_),V.o4_mapsequence4(3,k_),V.o4_mapsequence4(4,k_),:) = tmp;
      A(V.o4_mapsequence4(1,k_),V.o4_mapsequence4(2,k_),V.o4_mapsequence4(4,k_),V.o4_mapsequence4(3,k_),:) = tmp;
      A(V.o4_mapsequence4(2,k_),V.o4_mapsequence4(1,k_),V.o4_mapsequence4(3,k_),V.o4_mapsequence4(4,k_),:) = tmp;
      A(V.o4_mapsequence4(2,k_),V.o4_mapsequence4(1,k_),V.o4_mapsequence4(4,k_),V.o4_mapsequence4(3,k_),:) = tmp;
      A(V.o4_mapsequence4(3,k_),V.o4_mapsequence4(4,k_),V.o4_mapsequence4(1,k_),V.o4_mapsequence4(2,k_),:) = tmp;
      A(V.o4_mapsequence4(3,k_),V.o4_mapsequence4(4,k_),V.o4_mapsequence4(2,k_),V.o4_mapsequence4(1,k_),:) = tmp;
      A(V.o4_mapsequence4(4,k_),V.o4_mapsequence4(3,k_),V.o4_mapsequence4(1,k_),V.o4_mapsequence4(2,k_),:) = tmp;
      A(V.o4_mapsequence4(4,k_),V.o4_mapsequence4(3,k_),V.o4_mapsequence4(2,k_),V.o4_mapsequence4(1,k_),:) = tmp;
    end
  end

  function vC = Voigt_ddot_vA4m_vB2(vAm,vB)
    if isa(vB,'sym') || isa(vAm,'sym')
      vC = (V.o4_multiplicity_ddot.*vAm)*vB;
    else
      vC = reshape(pagemtimes(V.o4_multiplicity_ddot.*vAm,vB),V.o2_voidim,size(vB,2));
    end
  end

  function vC = Voigt_ddot_vA2_vB4m(vA,vBm)
    vC = reshape(pagemtimes(vA,'transpose',V.o4_multiplicity_ddot'.*vBm,'none'),V.o2_voidim,size(vA,2));
  end

  function D = Voigt_ddot_vA2_vB4m_vC2(vA,vBm,vC)
    D = reshape(pagemtimes(vA,'transpose',pagemtimes(V.o4_multiplicity_ddddot'.*vBm,vC),'none'),1,1,size(vA,2));
  end

  function C = Voigt_symotimes_A_B(A,B)
    C = zeros(3,3,size(A,2));
    if nargin == 1
      for i_ = 1:3, for j_ = 1:3   C(i_,j_,:) = A(i_,:)*A(j_,:);                 end, end
    else
      for i_ = 1:3, for j_ = 1:3   C(i_,j_,:) = A(i_,:)*B(j_,:)+B(i_,:)*A(j_,:); end, end
    end
  end

  function C = Voigt_symotimes_A2_B2(A,B)
    C = zeros(3,3,3,3,size(A,3));
    if nargin == 1
      for i_ = 1:3, for j_ = 1:3, for k_ = 1:3, for l_ = 1:3,   C(i_,j_,k_,l_,:) = A(i_,j_,:)*A(k_,l_,:);                       end, end, end, end
    else
      for i_ = 1:3, for j_ = 1:3, for k_ = 1:3, for l_ = 1:3,   C(i_,j_,k_,l_,:) = A(i_,j_,:)*B(k_,l_,:)+B(i_,j_,:)*A(k_,l_,:); end, end, end, end
    end
  end

  function vCm = Voigt_symotimes_vA2_vB2(vA,vB)
    % returns C = A x A
    %         C = A x B + B x A
    nq = size(vA,2);
    if nargin == 1
      vCm = pagemtimes(reshape(vA,V.o2_voidim,1,nq),reshape(vA,1,V.o2_voidim,nq));
    else
      vCm = pagemtimes(reshape(vA,V.o2_voidim,1,nq),reshape(vB,1,V.o2_voidim,nq)) + pagemtimes(reshape(vB,V.o2_voidim,1,nq),reshape(vA,1,V.o2_voidim,nq));
    end
  end

  function C = Voigt_fourotimes_N(N,i,j)
    % returns C = (Ni x Nj x Ni x Nj) + (Nj x Ni x Nj x Ni) + (Ni x Nj x Nj x Ni) + (Nj x Ni x Ni x Nj)
    nq = size(N,3);
    C = zeros(V.dim,V.dim,V.dim,V.dim,nq);
    if isa(N,'sym'), C = sym(C); end
    for i_ = 1:V.dim
      for j_ = 1:V.dim
        for k_ = 1:V.dim
          for l_ = 1:V.dim
            C(i_,j_,k_,l_,:) = N(i_,i,:)*N(j_,j,:)*N(k_,i,:)*N(l_,j,:) + N(i_,j,:)*N(j_,i,:)*N(k_,j,:)*N(l_,i,:) + N(i_,i,:)*N(j_,j,:)*N(k_,j,:)*N(l_,i,:) + N(i_,j,:)*N(j_,i,:)*N(k_,i,:)*N(l_,j,:);
          end
        end
      end
    end
  end

  function vCm = Voigt_symfourotimes_N(N,i,j)
    % returns C = (Ni x Nj x Ni x Nj) + (Nj x Ni x Nj x Ni) + (Ni x Nj x Nj x Ni) + (Nj x Ni x Ni x Nj)
    nq = size(N,3);
    vCm = zeros(V.o2_voidim,V.o2_voidim,nq);
    if isa(N,'sym'), vCm = sym(vCm); end
    switch V.dim
      case 2

      case 3
        vCm(1,1,:) = 2*N(1,i,:)*N(1,i,:)*N(1,j,:)*N(1,j,:);
        vCm(2,2,:) = 2*N(2,i,:)*N(2,i,:)*N(2,j,:)*N(2,j,:);
        vCm(3,3,:) = 2*N(3,i,:)*N(3,i,:)*N(3,j,:)*N(3,j,:);
        vCm(4,4,:) =   N(2,i,:)*N(2,i,:)*N(3,j,:)*N(3,j,:) + N(2,j,:)*N(3,i,:)*N(2,i,:)*N(3,j,:);
        vCm(5,5,:) =   N(2,i,:)*N(2,i,:)*N(3,j,:)*N(3,j,:) + N(2,j,:)*N(3,i,:)*N(2,i,:)*N(3,j,:);
    end
  end

  function C = Voigt_symodot_A2(A)
    C = zeros(3,3,3,3,size(A,3));
    if isa(A,'sym'), C = sym(C); end
    for i_ = 1:3
      for j_ = 1:3
        for k_ = 1:3
          for l_ = 1:3
            C(i_,j_,k_,l_,:) = (A(i_,k_)*A(j_,l_,:)+A(i_,l_)*A(j_,k_,:))/2;
          end
        end
      end
    end
  end

  function vCm = Voigt_symodot_vA2(vA)
    % returns C_ijkl = (A_ik * A_jl)/2
    nq = size(vA,2);
    vCm = zeros(V.o2_voidim,V.o2_voidim,nq);
    if isa(vA,'sym'), vCm = sym(vCm); end
    switch V.dim
      case 2
        vCm(1,1,:) =  vA(1,:).*vA(1,:);
        vCm(2,2,:) =  vA(2,:).*vA(2,:);
        vCm(1,2,:) = (vA(1,:).*vA(2,:)+vA(3,:).*vA(3,:))/2;   vCm(2,1,:) = vCm(1,2,:);
      case 3
        switch V.version
          case 1
            vCm(1,1,:) =  vA(1,:).*vA(1,:);
            vCm(1,2,:) =  vA(6,:).*vA(6,:);                       vCm(2,1,:) = vCm(1,2,:);
            vCm(1,3,:) =  vA(5,:).*vA(5,:);                       vCm(3,1,:) = vCm(1,3,:);
            vCm(1,4,:) =  vA(6,:).*vA(5,:);                       vCm(4,1,:) = vCm(1,4,:);
            vCm(1,5,:) =  vA(1,:).*vA(5,:);                       vCm(5,1,:) = vCm(1,5,:);
            vCm(1,6,:) =  vA(1,:).*vA(6,:);                       vCm(6,1,:) = vCm(1,6,:);
            vCm(2,2,:) =  vA(2,:).*vA(2,:);
            vCm(2,3,:) =  vA(4,:).*vA(4,:);                       vCm(3,2,:) = vCm(2,3,:);
            vCm(2,4,:) =  vA(2,:).*vA(4,:);                       vCm(4,2,:) = vCm(2,4,:);
            vCm(2,5,:) =  vA(6,:).*vA(4,:);                       vCm(5,2,:) = vCm(2,5,:);
            vCm(2,6,:) =  vA(6,:).*vA(2,:);                       vCm(6,2,:) = vCm(2,6,:);
            vCm(3,3,:) =  vA(3,:).*vA(3,:);
            vCm(3,4,:) =  vA(4,:).*vA(3,:);                       vCm(4,3,:) = vCm(3,4,:);
            vCm(3,5,:) =  vA(5,:).*vA(3,:);                       vCm(5,3,:) = vCm(3,5,:);
            vCm(3,6,:) =  vA(5,:).*vA(4,:);                       vCm(6,3,:) = vCm(3,6,:);
            vCm(4,4,:) = (vA(2,:).*vA(3,:)+vA(4,:).*vA(4,:))/2;
            vCm(4,5,:) = (vA(6,:).*vA(3,:)+vA(4,:).*vA(5,:))/2;   vCm(5,4,:) = vCm(4,5,:);
            vCm(4,6,:) = (vA(6,:).*vA(4,:)+vA(2,:).*vA(5,:))/2;   vCm(6,4,:) = vCm(4,6,:);
            vCm(5,5,:) = (vA(1,:).*vA(3,:)+vA(5,:).*vA(5,:))/2;
            vCm(5,6,:) = (vA(1,:).*vA(4,:)+vA(6,:).*vA(5,:))/2;   vCm(6,5,:) = vCm(5,6,:);
            vCm(6,6,:) = (vA(1,:).*vA(2,:)+vA(6,:).*vA(6,:))/2;
          case 2
            vCm(1,1,:) =  vA(1,:).*vA(1,:);
            vCm(1,2,:) =  vA(4,:).*vA(4,:);                       vCm(2,1,:) = vCm(1,2,:);
            vCm(1,3,:) =  vA(6,:).*vA(6,:);                       vCm(3,1,:) = vCm(1,3,:);
            vCm(1,5,:) =  vA(4,:).*vA(6,:);                       vCm(5,1,:) = vCm(1,5,:);
            vCm(1,6,:) =  vA(1,:).*vA(6,:);                       vCm(6,1,:) = vCm(1,6,:);
            vCm(1,4,:) =  vA(1,:).*vA(4,:);                       vCm(4,1,:) = vCm(1,4,:);
            vCm(2,2,:) =  vA(2,:).*vA(2,:);
            vCm(2,3,:) =  vA(5,:).*vA(5,:);                       vCm(3,2,:) = vCm(2,3,:);
            vCm(2,5,:) =  vA(2,:).*vA(5,:);                       vCm(5,2,:) = vCm(2,5,:);
            vCm(2,6,:) =  vA(4,:).*vA(5,:);                       vCm(6,2,:) = vCm(2,6,:);
            vCm(2,4,:) =  vA(4,:).*vA(2,:);                       vCm(4,2,:) = vCm(2,4,:);
            vCm(3,3,:) =  vA(3,:).*vA(3,:);
            vCm(3,5,:) =  vA(5,:).*vA(3,:);                       vCm(5,3,:) = vCm(3,5,:);
            vCm(3,6,:) =  vA(6,:).*vA(3,:);                       vCm(6,3,:) = vCm(3,6,:);
            vCm(3,4,:) =  vA(6,:).*vA(5,:);                       vCm(4,3,:) = vCm(3,4,:);
            vCm(5,5,:) = (vA(2,:).*vA(3,:)+vA(5,:).*vA(5,:))/2;
            vCm(5,6,:) = (vA(4,:).*vA(3,:)+vA(5,:).*vA(6,:))/2;   vCm(6,5,:) = vCm(5,6,:);
            vCm(5,4,:) = (vA(4,:).*vA(5,:)+vA(2,:).*vA(6,:))/2;   vCm(4,5,:) = vCm(5,4,:);
            vCm(6,6,:) = (vA(1,:).*vA(3,:)+vA(6,:).*vA(6,:))/2;
            vCm(6,4,:) = (vA(1,:).*vA(5,:)+vA(4,:).*vA(6,:))/2;   vCm(4,6,:) = vCm(6,4,:);
            vCm(4,4,:) = (vA(1,:).*vA(2,:)+vA(4,:).*vA(4,:))/2;
          case 3
            vCm(1,1,:) =  vA(1,:).*vA(1,:);
            vCm(1,2,:) =  vA(4,:).*vA(4,:);                       vCm(2,1,:) = vCm(1,2,:);
            vCm(1,3,:) =  vA(5,:).*vA(5,:);                       vCm(3,1,:) = vCm(1,3,:);
            vCm(1,6,:) =  vA(4,:).*vA(5,:);                       vCm(6,1,:) = vCm(1,6,:);
            vCm(1,5,:) =  vA(1,:).*vA(5,:);                       vCm(5,1,:) = vCm(1,5,:);
            vCm(1,4,:) =  vA(1,:).*vA(4,:);                       vCm(4,1,:) = vCm(1,4,:);
            vCm(2,2,:) =  vA(2,:).*vA(2,:);
            vCm(2,3,:) =  vA(6,:).*vA(6,:);                       vCm(3,2,:) = vCm(2,3,:);
            vCm(2,6,:) =  vA(2,:).*vA(6,:);                       vCm(6,2,:) = vCm(2,6,:);
            vCm(2,5,:) =  vA(4,:).*vA(6,:);                       vCm(5,2,:) = vCm(2,5,:);
            vCm(2,4,:) =  vA(4,:).*vA(2,:);                       vCm(4,2,:) = vCm(2,4,:);
            vCm(3,3,:) =  vA(3,:).*vA(3,:);
            vCm(3,6,:) =  vA(6,:).*vA(3,:);                       vCm(6,3,:) = vCm(3,6,:);
            vCm(3,5,:) =  vA(5,:).*vA(3,:);                       vCm(5,3,:) = vCm(3,5,:);
            vCm(3,4,:) =  vA(5,:).*vA(6,:);                       vCm(4,3,:) = vCm(3,4,:);
            vCm(6,6,:) = (vA(2,:).*vA(3,:)+vA(6,:).*vA(6,:))/2;
            vCm(6,5,:) = (vA(4,:).*vA(3,:)+vA(6,:).*vA(5,:))/2;   vCm(5,6,:) = vCm(6,5,:);
            vCm(6,4,:) = (vA(4,:).*vA(6,:)+vA(2,:).*vA(5,:))/2;   vCm(4,6,:) = vCm(6,4,:);
            vCm(5,5,:) = (vA(1,:).*vA(3,:)+vA(5,:).*vA(5,:))/2;
            vCm(5,4,:) = (vA(1,:).*vA(6,:)+vA(4,:).*vA(5,:))/2;   vCm(4,5,:) = vCm(5,4,:);
            vCm(4,4,:) = (vA(1,:).*vA(2,:)+vA(4,:).*vA(4,:))/2;
        end
    end
  end



if nargin ==1,   version = 1; end
V = struct(...
  'dim',dim, 'version',version, ...
  'map_vA2_A2',[], 'map_A2_vA2',[], 'o2_mapsequence',[], 'o2_voidim',[], 'map_vA2_A2_inversemapping',[], 'o2_multiplicity',[], ...
  'map_vA4m_A4',[], 'map_A4_vA4m',[], 'map_vA4m_A4_mapping',[], 'o4_vdim',[], 'map2to3_inversemapping',[], 'map_vA4m_A4_multiplicity',[], ...
  'map_vA4_A4',[],...
  'ddot_vA2_vB2',[]...
  );

%% map_vA2_A2
switch dim
  case 2
    V.o2_voidim = 3;
    V.o2_mapsequence = [...
      1 2 1;...
      1 2 2];
    switch version
      case 1,      o2_versionpermute = [1 2 3];
      otherwise,   fprintf('voigt_init(): dim==2, version can only be 1\n');   return;
    end
  case 3
    V.o2_voidim = 6;
    V.o2_mapsequence = [...
      1 2 3 2 1 1;...
      1 2 3 3 3 2];
    switch version
      case 1,      o2_versionpermute = [1 2 3 4 5 6];   o2_versionipermute = [1 2 3 4 5 6];
      case 2,      o2_versionpermute = [1 2 3 6 4 5];   o2_versionipermute = [1 2 3 5 6 4];
      case 3,      o2_versionpermute = [1 2 3 6 5 4];   o2_versionipermute = [1 2 3 6 5 4];
      otherwise,   fprintf('voigt_init(): dim==3, version can only be 1 or 2 or 3\n');  return;
    end
    V.o2_mapsequence = V.o2_mapsequence(:,o2_versionpermute);
  otherwise
    fprintf('voigt_init(): dim can only be 2 or 3\n');
    return;
end
V.map_vA2_A2    = @Voigt_map_vA2_A2;
V.map_vA2_A2_inversemapping = V.map_vA2_A2((1:V.o2_voidim)');
tmp = unique(V.map_vA2_A2_inversemapping);
V.o2_multiplicity = histcounts(V.map_vA2_A2_inversemapping(:), [tmp; max(tmp)+1]);
V.map_A2_vA2 = @Voigt_map_A2_vA2;
V.ddot_vA2_vB2  = @Voigt_ddot_vA2_vB2;

%% map_vA4m_A4
switch dim
  case 2
    V.o4_vdim  =  6;
    V.o4_mapsequence2 = [...
      1 2 3 2 1 1;... % I = voigt{ij}
      1 2 3 3 3 2;... % K = voigt{kl}
      ];
    V.o4_mapsequence4 = zeros(4,6);
    V.o4_multiplicity_ddot   = [1 1 1]'*[1 1 2];
    V.o4_multiplicity_ddddot = [1 1 2]'*[1 1 2];
    switch version
      case 1,      o4_versionpermute = [1 2 3 4 5 6];
      otherwise,   fprintf('voigt_init(): dim==2, version can only be 1\n');   return;
    end
    V.o4_mapsequence2 = V.o4_mapsequence2(:,o4_versionpermute);
    for k = 1:V.o4_vdim
      I = V.o4_mapsequence2(1,k);
      J = V.o4_mapsequence2(2,k);
      V.o4_mapsequence4(:,k) = [V.o2_mapsequence(:,I); V.o2_mapsequence(:,J)];
    end
  case 3
    V.o4_vdim     = 21;
    V.o4_mapsequence2 = [...
      %1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1
      1 2 3 4 5 6 5 4 3 2 1 1 1 1 1 2 3 4 3 2 2;... % I = voigt{ij}
      1 2 3 4 5 6 6 6 6 6 6 5 4 3 2 3 4 5 5 5 4;... % K = voigt{kl}
      ];
    V.o4_mapsequence4 = zeros(4,21);
    %V.map_vA4m_A4_mapping = [...
    %  1 2 3 2 3 1 3 2 3 2 1 1 1 1 1 2 3 2 3 2 2;... % i                [1111 1122 1133 1123 1131 1112]
    %  1 2 3 3 1 2 1 3 3 2 1 1 1 1 1 2 3 3 3 2 2;... % j                [     2222 2233 2223 2231 2212]
    %  1 2 3 2 3 1 1 1 1 1 1 3 2 3 2 3 2 3 3 3 2;... % k                [          3333 3323 3331 3312]
    %  1 2 3 3 1 2 2 2 2 2 2 1 3 3 2 3 3 1 1 1 3;... % l                [               2323 2331 2312]
    %  1 2 3 4 5 6 5 4 3 2 1 1 1 1 1 2 3 4 3 2 2;... % I = voigt{ij}    [                    3131 3112]
    %  1 2 3 4 5 6 6 6 6 6 6 5 4 3 2 3 4 5 5 5 4;... % K = voigt{kl}    [ sym                     1212]
    %  ];
    V.o4_multiplicity_ddot   = [1 1 1 1 1 1]'*[1 1 1 2 2 2];
    V.o4_multiplicity_ddddot = [1 1 1 2 2 2]'*[1 1 1 2 2 2];
    switch version
      case 1,      o4_versionpermute = [1 2 3 4 5 6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21];
      case 2,      o4_versionpermute = [1 2 3 4 5 6 15 16 17 18  7 14 21 19  8 13 20  9 12 10 11];
      case 3,      o4_versionpermute = [1 2 3 4 5 6 15 14 13 12 11 16 21 20 10 17 19  9 18  8  7];
      otherwise,   fprintf('voigt_init(): dim==3, version can only be 1 or 2 or 3\n');  return;
    end
    V.o4_mapsequence2 = V.o4_mapsequence2(:,o4_versionpermute);
    for k = 1:V.o4_vdim
      I = V.o4_mapsequence2(1,k);
      J = V.o4_mapsequence2(2,k);
      V.o4_mapsequence4(:,k) = [V.o2_mapsequence(:,I); V.o2_mapsequence(:,J)];
    end
  otherwise
    fprintf('voigt_init(): dim can only be 2 or 3\n');
    return;
end

V.map_vA4m_A4_inversemapping.I = zeros(V.dim,V.dim,V.dim,V.dim);
V.map_vA4m_A4_inversemapping.J = zeros(V.dim,V.dim,V.dim,V.dim);
% V.map_vA4m_A4_vecmultiplicity = zeros(1,V.o4_vdim);
% V.map_vA4m_A4_matmultiplicity = zeros(V.o4_mdim,V.o4_mdim);
% for m = 1:V.o4_vdim
%   i_ = V.map_vA4m_A4_mapping(1,m);
%   j_ = V.map_vA4m_A4_mapping(2,m);
%   k_ = V.map_vA4m_A4_mapping(3,m);
%   l_ = V.map_vA4m_A4_mapping(4,m);
%   I_ = V.map_vA4m_A4_mapping(5,m);
%   J_ = V.map_vA4m_A4_mapping(6,m);
%   p_ = [...
%     i_ j_ k_ l_;...
%     i_ j_ l_ k_;...
%     j_ i_ k_ l_;...
%     j_ i_ l_ k_;...
%     k_ l_ i_ j_;...
%     k_ l_ j_ i_;...
%     l_ k_ i_ j_;...
%     l_ k_ j_ i_]; % all symmetric permutations
%   p_ = unique(p_, 'rows');
%   np = size(p_,1);
%   for n = 1:np
%     V.map_vA4m_A4_inversemapping.I(p_(n,1),p_(n,2),p_(n,3),p_(n,4)) = I_;
%     V.map_vA4m_A4_inversemapping.J(p_(n,1),p_(n,2),p_(n,3),p_(n,4)) = J_;
%   end
%   V.map_vA4m_A4_vecmultiplicity(1,m) = np;
%   if I_ == J_
%     V.map_vA4m_A4_matmultiplicity(I_,J_) = np;
%   else
%     V.map_vA4m_A4_matmultiplicity(I_,J_) = np/2;
%     V.map_vA4m_A4_matmultiplicity(J_,I_) = np/2;
%   end
% end
V.map_vA4m_A4       = @Voigt_map_vA4m_A4;
V.map_A4_vA4m       = @Voigt_map_A4_vA4m;
V.map_A4_vA4        = @Voigt_map_A4_vA4;
V.map_vA4_A4        = @Voigt_map_vA4_A4;
V.ddot_vA4m_vB2     = @Voigt_ddot_vA4m_vB2;
V.ddot_vA2_vB4m     = @Voigt_ddot_vA2_vB4m;
V.ddot_vA2_vB4m_vC2 = @Voigt_ddot_vA2_vB4m_vC2;
V.symotimes_A_B     = @Voigt_symotimes_A_B;
V.symotimes_A2_B2   = @Voigt_symotimes_A2_B2;
V.symotimes_vA2_vB2 = @Voigt_symotimes_vA2_vB2;
V.symodot_A2        = @Voigt_symodot_A2;
V.symodot_vA2       = @Voigt_symodot_vA2;
V.fourotimes_N      = @Voigt_fourotimes_N;
V.symfourotimes_N   = @Voigt_symfourotimes_N;
end


