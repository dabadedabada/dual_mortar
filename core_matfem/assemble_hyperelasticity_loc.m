function [M_loc, Kc_loc, Ks_loc, f_loc, res_loc, byprod] = assemble_hyperelasticity_loc(Xloc,uloc,FE,consts,problem)
% Inputs
%   Xloc   ... n_{nodes}    x dim           matrix ... point coordinates
%   FE     ...                                         Finite Element (shape functions + Gauss points & weights)
%   consts ...                                         physics
% Outputs
%   M  ...      mass    matrix  \int_{\Omega}           u   .      v   d\Omega
%   K  ... stiffness    matrix  \int_{\Omega}      grad(u)  . grad(v)  d\Omega
%   f  ... volume force vector  \int_{\Omega}      consts.? .      v   d\Omega
  function y = det3cell(x)
    y = x{1,1}.*x{2,2}.*x{3,3} + x{1,3}.*x{2,1}.*x{3,2} + x{1,2}.*x{2,3}.*x{3,1} - x{1,3}.*x{2,2}.*x{3,1} - x{1,1}.*x{2,3}.*x{3,2} - x{1,2}.*x{2,1}.*x{3,3};
  end
  function y = det33k(x)
    y = x(1,1,:).*x(2,2,:).*x(3,3,:) + x(1,3,:).*x(2,1,:).*x(3,2,:) + x(1,2,:).*x(2,3,:).*x(3,1,:) - x(1,3,:).*x(2,2,:).*x(3,1,:) - x(1,1,:).*x(2,3,:).*x(3,2,:) - x(1,2,:).*x(2,1,:).*x(3,3,:);
  end
  function y = inv3cell(x,detx)
    y = {...
      (x{2,2}.*x{3,3}-x{2,3}.*x{3,2})./detx, (x{1,3}.*x{3,2}-x{1,2}.*x{3,3})./detx, (x{1,2}.*x{2,3}-x{1,3}.*x{2,2})./detx;...
      (x{2,3}.*x{3,1}-x{2,1}.*x{3,3})./detx, (x{1,1}.*x{3,3}-x{1,3}.*x{3,1})./detx, (x{1,3}.*x{2,1}-x{1,1}.*x{2,3})./detx;...
      (x{2,1}.*x{3,2}-x{2,2}.*x{3,1})./detx, (x{1,2}.*x{3,1}-x{1,1}.*x{3,2})./detx, (x{1,1}.*x{2,2}-x{1,2}.*x{2,1})./detx};
  end
  function y = inv33k(x,detx)
    y = zeros(3,3,size(x,3));
    detx = reshape(detx,1,1,size(detx,2));
    y(1,1,:) = (x(2,2,:).*x(3,3,:)-x(2,3,:).*x(3,2,:))./detx;   y(1,2,:) = (x(1,3,:).*x(3,2,:)-x(1,2,:).*x(3,3,:))./detx;   y(1,3,:) = (x(1,2,:).*x(2,3,:)-x(1,3,:).*x(2,2,:))./detx;
    y(2,1,:) = (x(2,3,:).*x(3,1,:)-x(2,1,:).*x(3,3,:))./detx;   y(2,2,:) = (x(1,1,:).*x(3,3,:)-x(1,3,:).*x(3,1,:))./detx;   y(2,3,:) = (x(1,3,:).*x(2,1,:)-x(1,1,:).*x(2,3,:))./detx;
    y(3,1,:) = (x(2,1,:).*x(3,2,:)-x(2,2,:).*x(3,1,:))./detx;   y(3,2,:) = (x(1,2,:).*x(3,1,:)-x(1,1,:).*x(3,2,:))./detx;   y(3,3,:) = (x(1,1,:).*x(2,2,:)-x(1,2,:).*x(2,1,:))./detx;
  end

dim = size(FE.QP,1);   nqp = size(FE.QP,2);   nN  = size(FE.xi,1);
% dofs
N__loc    = FE.N(FE.QP);                   %  nN   x nqp        matrix
dNdxi_loc = FE.dNdxi(FE.QP);               %  nN   x nqp*dim    matrix [[dNdxi1(qp1) dNdxi1(qp2) ...] [dNdxi2(qp1) dNdxi2(qp2) ...] {dNdxi3(qp1) dNdxi3(qp2) ...}]

% i__loc_mat = 1:nN;
% I__loc_mat = (1:nN)'*ones(1,nN);             %  nN   x  nN        matrix (for local  I J combinations)
% J__loc_mat = I__loc_mat';                    %  nN   x  nN        matrix (for local  I J combinations)
% I__loc_vec = reshape( I__loc_mat,[],1);      %  nN^2 x   1        vector (for local  I J combinations)
% J__loc_vec = reshape( J__loc_mat,[],1);      %  nN^2 x   1        vector (for local  I J combinations)

% NI_loc_mat = N__loc_mat(I__loc_vec,:);       %  nN^2 x nqp        matrix
% NJ_loc_mat = N__loc_mat(J__loc_vec,:);       %  nN^2 x nqp        matrix
% W_loc      = FE.QW';                         %  nN   x 1          vector



%%% E  = consts.Young_modulus;
%%% nu = consts.Poissons_ratio;
%C__ = E/(1 - nu^2)*[1, nu, 0; nu, 1, 0; 0,   0, (1 - nu)/2];
switch dim
  case 2
    %     dNdxi1_loc = dNdxi_loc(:,1:end/2);   dNdxi2_loc = dNdxi_loc(:,end/2+1:  end  );
    %     J11_loc  = Xloc(1,:)*dNdxi1_loc;   J12_loc = Xloc(1,:)*dNdxi2_loc;
    %     J21_loc  = Xloc(2,:)*dNdxi1_loc;   J22_loc = Xloc(2,:)*dNdxi2_loc;
    %     detJ_loc = J11_loc.*J22_loc - J12_loc.*J21_loc;
    %     absdetJ_loc_mat = abs(detJ_loc);
    %     Jinv11_loc =  J22_loc./detJ_loc;   Jinv12_loc = -J12_loc./detJ_loc;
    %     Jinv21_loc = -J21_loc./detJ_loc;   Jinv22_loc =  J11_loc./detJ_loc;
    %     dNdX1_loc  = dNdxi1_loc .* (ones(nN,1)*Jinv11_loc) + dNdxi2_loc .* (ones(nN,1)*Jinv21_loc);
    %     dNdX2_loc  = dNdxi1_loc .* (ones(nN,1)*Jinv12_loc) + dNdxi2_loc .* (ones(nN,1)*Jinv22_loc);
    %     dNIdX1_loc_mat = dNdX1_loc(I__loc_vec,:);    dNIdX2_loc_mat = dNdX2_loc(I__loc_vec,:);
    %     dNJdX1_loc_mat = dNdX1_loc(J__loc_vec,:);    dNJdX2_loc_mat = dNdX2_loc(J__loc_vec,:);
    %     NI_NJ_loc_mat         = (    NI_loc_mat.*    NJ_loc_mat)*(absdetJ_loc_mat'.*W_loc);
    %     Ni_loc_mat            = N__loc_mat*(absdetJ_loc_mat'.*W_loc);
    %     dNIdX1_dNJdX1_loc_mat = (dNIdX1_loc_mat.*dNJdX1_loc_mat)*(absdetJ_loc_mat'.*W_loc);
    %     dNIdX1_dNJdX2_loc_mat = (dNIdX1_loc_mat.*dNJdX2_loc_mat)*(absdetJ_loc_mat'.*W_loc);
    %     dNIdX2_dNJdX1_loc_mat = (dNIdX2_loc_mat.*dNJdX1_loc_mat)*(absdetJ_loc_mat'.*W_loc);
    %     dNIdX2_dNJdX2_loc_mat = (dNIdX2_loc_mat.*dNJdX2_loc_mat)*(absdetJ_loc_mat'.*W_loc);
    %     MM_IJ    = consts.Density*NI_NJ_loc_mat;
    %     KK_IJ_11 = (E/(1 - nu^2))*(dNIdX1_dNJdX1_loc_mat            + dNIdX2_dNJdX2_loc_mat*((1-nu)/2));
    %     KK_IJ_22 = (E/(1 - nu^2))*(dNIdX1_dNJdX1_loc_mat*((1-nu)/2) + dNIdX2_dNJdX2_loc_mat           );
    %     KK_IJ_12 = (E/(1 - nu^2))*(dNIdX1_dNJdX2_loc_mat*nu         + dNIdX2_dNJdX1_loc_mat*((1-nu)/2));    %KK_IJ_21 = (E/(1 - nu^2))*(dNIdX2_dNJdX1_loc_mat*nu + ((1-nu)/2)*dNIdX1_dNJdX2_loc_mat);
    %     ff_i_1   = problem.vol(1)*Ni_loc_mat;
    %     ff_i_2   = problem.vol(2)*Ni_loc_mat;
    %     M  = sparse(...
    %       [2*(I__loc_vec-1)+1; 2*(I__loc_vec-1)+2],...
    %       [2*(J__loc_vec-1)+1; 2*(J__loc_vec-1)+2],...
    %       [             MM_IJ;              MM_IJ]);
    %     K = sparse(...
    %       [2*(I__loc_vec-1)+1; 2*(I__loc_vec-1)+1;2*(J__loc_vec-1)+2; 2*(I__loc_vec-1)+2],...
    %       [2*(J__loc_vec-1)+1; 2*(J__loc_vec-1)+2;2*(I__loc_vec-1)+1; 2*(J__loc_vec-1)+2],...
    %       [          KK_IJ_11;           KK_IJ_12;          KK_IJ_12;           KK_IJ_22]);
    %     f = sparse([2*(i__loc_mat-1)+1;2*(i__loc_mat-1)+2],[0*i__loc_mat+1;0*i__loc_mat+1],[ff_i_1;ff_i_2]);
  case 3
    % uloc, Xloc ... 3 x N
    % dNdxi_loc  ... N x 3 x q
    dNdxi_loc = permute(reshape(dNdxi_loc,nN,nqp,3),[1 3 2]);      % N . 3 . q
    dXdxi_loc = reshape(Xloc*reshape(dNdxi_loc,nN,3*nqp),3,3,nqp); % 3 . 3 . q
    detdXdxi_loc = reshape(det33k(dXdxi_loc),1,nqp);               % q
    absJW = abs(detdXdxi_loc)'.*FE.QW';                            % q
    dxidX_loc = inv33k(dXdxi_loc,detdXdxi_loc);                    % 3 . 3 . q
    dNdX_loc = pagemtimes(dNdxi_loc,dxidX_loc);                    % N . 3 . q
    Z_loc  = reshape(uloc*reshape(dNdX_loc,nN,3*nqp),3,3,nqp);     % 3 . 3 . q
    F_loc  = Z_loc;                                                % 3 . 3 . q
    F_loc(1,1,:) = F_loc(1,1,:)+1;   F_loc(2,2,:) = F_loc(2,2,:)+1;   F_loc(3,3,:) = F_loc(3,3,:)+1;
    C_loc  = pagemtimes(F_loc,'transpose',F_loc,'none') ;          % 3 . 3 . q
    C_tmp  = reshape([C_loc(1,1,:); C_loc(2,2,:); C_loc(3,3,:); C_loc(1,2,:); C_loc(2,3,:); C_loc(1,3,:)],6,nqp);
    V = voigt_init(3,2);
    switch problem.material_model                                  % 6 . 6 . q (CC),   6 . q (S)
      case {'arruda_boyce','arruda-boyce'}      ,    [CC_tmp,S_tmp] =    arruda_boyce(C_tmp,consts.lambda_L,consts.d,consts.mu);
      case {'blatz_ko','blatz-ko'}              ,    [CC_tmp,S_tmp] =        blatz_ko(C_tmp,consts.mu);
      case  'gent'                              ,    [CC_tmp,S_tmp] =            gent(C_tmp,consts.mu,consts.Jm,consts.dmu);
      case  'kirchhoff'                         ,    [CC_tmp,S_tmp] =       kirchhoff(C_tmp,consts.lambda,consts.mu);
      %case  'kirchhoff'                         ,    [vC4m_loc,vS_loc] =     n_kirchhoff(V,V.map_A2_vA2(C_loc),consts.lambda,consts.mu); S_tmp = vS_loc; CC_tmp = V.map_A4_vA4(V.map_vA4m_A4(vC4m_loc));
      case {'mooney_rivlin_2','mooney-rivlin-2'},    [CC_tmp,S_tmp] = mooney_rivlin_2(C_tmp,consts.c,consts.d);
      case {'mooney_rivlin_3','mooney-rivlin-3'},    [CC_tmp,S_tmp] = mooney_rivlin_3(C_tmp,consts.c,consts.d);
      case {'mooney_rivlin_5','mooney-rivlin-5'},    [CC_tmp,S_tmp] = mooney_rivlin_5(C_tmp,consts.c,consts.d);
      case {'mooney_rivlin_9','mooney-rivlin-9'},    [CC_tmp,S_tmp] = mooney_rivlin_9(C_tmp,consts.c,consts.d);
      %case {'neo_hookean_cmp','neo-hookean-cmp'},    [CC_tmp,S_tmp] = neo_hookean_cmp(C_tmp,consts.mu,consts.d);
      case {'neo_hookean_cmp','neo-hookean-cmp'},    [vC4m_loc,vS_loc] =     n_neo_hookean_cmp(V,V.map_A2_vA2(C_loc),consts.lambda,consts.mu); S_tmp = vS_loc; CC_tmp = V.map_A4_vA4(V.map_vA4m_A4(vC4m_loc));
      case {'neo_hookean_inc','neo-hookean-inc'},    [CC_tmp,S_tmp] = neo_hookean_inc(C_tmp,consts.lambda,consts.mu);
      case {'ogden_1','ogden-1'}                ,    [CC_tmp,S_tmp] =         ogden_1(C_tmp,consts.mu,consts.alpha,consts.d);
      case {'ogden_2','ogden-2'}                ,    [CC_tmp,S_tmp] =         ogden_2(C_tmp,consts.mu,consts.alpha,consts.d);
      case {'ogden_3','ogden-3'}                ,    [CC_tmp,S_tmp] =         ogden_3(C_tmp,consts.mu,consts.alpha,consts.d);
      case {'ogden_cmp_1','ogden-cmp-1'}        ,    [CC_tmp,S_tmp] =     ogden_cmp_1(C_tmp,consts.mu,consts.alpha,consts.beta);
      case {'ogden_cmp_2','ogden-cmp-2'}        ,    [CC_tmp,S_tmp] =     ogden_cmp_2(C_tmp,consts.mu,consts.alpha,consts.beta);
      case {'ogden_cmp_3','ogden-cmp-3'}        ,    [CC_tmp,S_tmp] =     ogden_cmp_3(C_tmp,consts.mu,consts.alpha,consts.beta);
      case  'yeoh'                              ,    [CC_tmp,S_tmp] =            yeoh(C_tmp,consts.mu,consts.alpha,consts.beta);
    end
    S_loc  = zeros(3,3,nqp);    % no Voigt
    CC_loc = zeros(6,6,nqp);    %    Voigt
    S_loc(1,1,:) = S_tmp(1,:);   S_loc(1,2,:) = S_tmp(4,:);   S_loc(1,3,:) = S_tmp(6,:);
    S_loc(2,1,:) = S_tmp(4,:);   S_loc(2,2,:) = S_tmp(2,:);   S_loc(2,3,:) = S_tmp(5,:);
    S_loc(3,1,:) = S_tmp(6,:);   S_loc(3,2,:) = S_tmp(5,:);   S_loc(3,3,:) = S_tmp(3,:);
    CC_loc(1,1,:) = CC_tmp( 1,:);   CC_loc(1,2,:) = CC_tmp( 7,:);   CC_loc(1,3,:) = CC_tmp(12,:);   CC_loc(1,4,:) = CC_tmp(16,:);   CC_loc(1,5,:) = CC_tmp(19,:);   CC_loc(1,6,:) = CC_tmp(21,:);
    CC_loc(2,1,:) = CC_tmp( 7,:);   CC_loc(2,2,:) = CC_tmp( 2,:);   CC_loc(2,3,:) = CC_tmp( 8,:);   CC_loc(2,4,:) = CC_tmp(13,:);   CC_loc(2,5,:) = CC_tmp(17,:);   CC_loc(2,6,:) = CC_tmp(20,:);
    CC_loc(3,1,:) = CC_tmp(12,:);   CC_loc(3,2,:) = CC_tmp( 8,:);   CC_loc(3,3,:) = CC_tmp( 3,:);   CC_loc(3,4,:) = CC_tmp( 9,:);   CC_loc(3,5,:) = CC_tmp(14,:);   CC_loc(3,6,:) = CC_tmp(18,:);
    CC_loc(4,1,:) = CC_tmp(16,:);   CC_loc(4,2,:) = CC_tmp(13,:);   CC_loc(4,3,:) = CC_tmp( 9,:);   CC_loc(4,4,:) = CC_tmp( 4,:);   CC_loc(4,5,:) = CC_tmp(10,:);   CC_loc(4,6,:) = CC_tmp(15,:);
    CC_loc(5,1,:) = CC_tmp(19,:);   CC_loc(5,2,:) = CC_tmp(17,:);   CC_loc(5,3,:) = CC_tmp(14,:);   CC_loc(5,4,:) = CC_tmp(10,:);   CC_loc(5,5,:) = CC_tmp( 5,:);   CC_loc(5,6,:) = CC_tmp(11,:);
    CC_loc(6,1,:) = CC_tmp(21,:);   CC_loc(6,2,:) = CC_tmp(20,:);   CC_loc(6,3,:) = CC_tmp(18,:);   CC_loc(6,4,:) = CC_tmp(15,:);   CC_loc(6,5,:) = CC_tmp(11,:);   CC_loc(6,6,:) = CC_tmp( 6,:);
    B_loc  = zeros(6,nN,3,nqp);                                    % 6 . N . q . 3
    B_loc(1,:,1,:) = reshape( (reshape(F_loc(1,1,:),1,nqp)) .* reshape(dNdX_loc(:,1,:),nN,nqp)                                                                   ,1,nN,nqp);
    B_loc(1,:,2,:) = reshape( (reshape(F_loc(2,1,:),1,nqp)) .* reshape(dNdX_loc(:,1,:),nN,nqp)                                                                   ,1,nN,nqp);
    B_loc(1,:,3,:) = reshape( (reshape(F_loc(3,1,:),1,nqp)) .* reshape(dNdX_loc(:,1,:),nN,nqp)                                                                   ,1,nN,nqp);
    B_loc(2,:,1,:) = reshape( (reshape(F_loc(1,2,:),1,nqp)) .* reshape(dNdX_loc(:,2,:),nN,nqp)                                                                   ,1,nN,nqp);
    B_loc(2,:,2,:) = reshape( (reshape(F_loc(2,2,:),1,nqp)) .* reshape(dNdX_loc(:,2,:),nN,nqp)                                                                   ,1,nN,nqp);
    B_loc(2,:,3,:) = reshape( (reshape(F_loc(3,2,:),1,nqp)) .* reshape(dNdX_loc(:,2,:),nN,nqp)                                                                   ,1,nN,nqp);
    B_loc(3,:,1,:) = reshape( (reshape(F_loc(1,3,:),1,nqp)) .* reshape(dNdX_loc(:,3,:),nN,nqp)                                                                   ,1,nN,nqp);
    B_loc(3,:,2,:) = reshape( (reshape(F_loc(2,3,:),1,nqp)) .* reshape(dNdX_loc(:,3,:),nN,nqp)                                                                   ,1,nN,nqp);
    B_loc(3,:,3,:) = reshape( (reshape(F_loc(3,3,:),1,nqp)) .* reshape(dNdX_loc(:,3,:),nN,nqp)                                                                   ,1,nN,nqp);
    B_loc(4,:,1,:) = reshape( (reshape(F_loc(1,1,:),1,nqp)) .* reshape(dNdX_loc(:,2,:),nN,nqp) + (reshape(F_loc(1,2,:),1,nqp)) .* reshape(dNdX_loc(:,1,:),nN,nqp),1,nN,nqp);
    B_loc(4,:,2,:) = reshape( (reshape(F_loc(2,1,:),1,nqp)) .* reshape(dNdX_loc(:,2,:),nN,nqp) + (reshape(F_loc(2,2,:),1,nqp)) .* reshape(dNdX_loc(:,1,:),nN,nqp),1,nN,nqp);
    B_loc(4,:,3,:) = reshape( (reshape(F_loc(3,1,:),1,nqp)) .* reshape(dNdX_loc(:,2,:),nN,nqp) + (reshape(F_loc(3,2,:),1,nqp)) .* reshape(dNdX_loc(:,1,:),nN,nqp),1,nN,nqp);
    B_loc(5,:,1,:) = reshape( (reshape(F_loc(1,2,:),1,nqp)) .* reshape(dNdX_loc(:,3,:),nN,nqp) + (reshape(F_loc(1,3,:),1,nqp)) .* reshape(dNdX_loc(:,2,:),nN,nqp),1,nN,nqp);
    B_loc(5,:,2,:) = reshape( (reshape(F_loc(2,2,:),1,nqp)) .* reshape(dNdX_loc(:,3,:),nN,nqp) + (reshape(F_loc(2,3,:),1,nqp)) .* reshape(dNdX_loc(:,2,:),nN,nqp),1,nN,nqp);
    B_loc(5,:,3,:) = reshape( (reshape(F_loc(3,2,:),1,nqp)) .* reshape(dNdX_loc(:,3,:),nN,nqp) + (reshape(F_loc(3,3,:),1,nqp)) .* reshape(dNdX_loc(:,2,:),nN,nqp),1,nN,nqp);
    B_loc(6,:,1,:) = reshape( (reshape(F_loc(1,1,:),1,nqp)) .* reshape(dNdX_loc(:,3,:),nN,nqp) + (reshape(F_loc(1,3,:),1,nqp)) .* reshape(dNdX_loc(:,1,:),nN,nqp),1,nN,nqp);
    B_loc(6,:,2,:) = reshape( (reshape(F_loc(2,1,:),1,nqp)) .* reshape(dNdX_loc(:,3,:),nN,nqp) + (reshape(F_loc(2,3,:),1,nqp)) .* reshape(dNdX_loc(:,1,:),nN,nqp),1,nN,nqp);
    B_loc(6,:,3,:) = reshape( (reshape(F_loc(3,1,:),1,nqp)) .* reshape(dNdX_loc(:,3,:),nN,nqp) + (reshape(F_loc(3,3,:),1,nqp)) .* reshape(dNdX_loc(:,1,:),nN,nqp),1,nN,nqp);
    Kc_loc  =      reshape(reshape(pagemtimes(reshape(B_loc,6,nN*3,nqp),'transpose', pagemtimes(CC_loc,'none',reshape(B_loc,6,nN*3,nqp),     'none'),'none'),(3*nN)^2,nqp)*absJW,3*nN,3*nN);
    Ks_loc  = kron(reshape(reshape(pagemtimes(                dNdX_loc ,     'none', pagemtimes( S_loc,'none',                 dNdX_loc,'transpose'),'none'),(  nN)^2,nqp)*absJW,  nN,  nN),eye(3));
    M_loc   = consts.density*kron((N__loc.*absJW')*N__loc',eye(3));
    floc    = kron(ones(nN,1),problem.volume_force');
    f_loc   = M_loc*floc;
    %e_loc   = (1/2)*(C_loc-eye(3));
    %e_locv  = [e_loc(1,1,:); e_loc(2,2,:); e_loc(3,3,:); 2*e_loc(1,2,:); 2*e_loc(2,3,:); 2*e_loc(1,3,:)];
    %S_locv  = pagemtimes(CC_loc,'none',e_locv,'none');
    S_locv  =  [S_loc(1,1,:); S_loc(2,2,:); S_loc(3,3,:); S_loc(1,2,:); S_loc(2,3,:); S_loc(1,3,:)];
    res_loc = reshape(pagemtimes(reshape(B_loc,6,3*nN,nqp),'transpose', S_locv,'none'),3*nN,nqp)*absJW;
    % permute
    ind = reshape([1:nN; nN+1:2*nN; 2*nN+1:3*nN],[],1);
    %Ks_loc = Ks_loc(ind,ind);
    Kc_loc  = Kc_loc(ind,ind);
    res_loc = res_loc(ind,:);
    byprod.S = [sum(S_loc(1,1,:)); sum(S_loc(2,2,:)); sum(S_loc(3,3,:)); sum(S_loc(1,2,:)); sum(S_loc(1,3,:)); sum(S_loc(2,3,:))]/nqp;
    % save('olda.mat','Kc_loc','Ks_loc','B_loc','dNdX_loc','S_loc','CC_loc','absJW');
    % a = 1;
    %M_loc  =  M_loc(ind,ind);
    %f_loc  =  f_loc(ind,  1);
%     ff_i_1   = problem.vol(1)*Ni_loc_mat;
%     ff_i_2   = problem.vol(2)*Ni_loc_mat;
%     ff_i_3   = problem.vol(3)*Ni_loc_mat;
%     M  = sparse(...
%       [3*(I__loc_vec-1)+1; 3*(I__loc_vec-1)+2; 3*(I__loc_vec-1)+3],...
%       [3*(J__loc_vec-1)+1; 3*(J__loc_vec-1)+2; 3*(J__loc_vec-1)+3],...
%       [             MM_IJ;              MM_IJ;              MM_IJ]);
%     K = sparse(...
%       [3*(I__loc_vec-1)+1; 3*(I__loc_vec-1)+2; 3*(I__loc_vec-1)+3; 3*(I__loc_vec-1)+1; 3*(J__loc_vec-1)+2; 3*(I__loc_vec-1)+1; 3*(J__loc_vec-1)+3; 3*(I__loc_vec-1)+2; 3*(J__loc_vec-1)+3],...
%       [3*(J__loc_vec-1)+1; 3*(J__loc_vec-1)+2; 3*(J__loc_vec-1)+3; 3*(J__loc_vec-1)+2; 3*(I__loc_vec-1)+1; 3*(J__loc_vec-1)+3; 3*(I__loc_vec-1)+1; 3*(J__loc_vec-1)+3; 3*(I__loc_vec-1)+2],...
%       [          KK_IJ_11;           KK_IJ_22;           KK_IJ_33;           KK_IJ_12;           KK_IJ_12;           KK_IJ_13;           KK_IJ_13;           KK_IJ_23;           KK_IJ_23]);
%     f = sparse([3*(i__loc_mat-1)+1;3*(i__loc_mat-1)+2;3*(i__loc_mat-1)+3],[0*i__loc_mat+1;0*i__loc_mat+1;0*i__loc_mat+1],[ff_i_1;ff_i_2;ff_i_3]);
end
end
