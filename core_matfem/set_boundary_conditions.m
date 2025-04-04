function [Bc] = set_boundary_conditions(mesh,problem,time,indtime,inu)

Bc.Bd  = [];
Bc.U0  = zeros(mesh.info.noddof*mesh.info.nodcount,1);

% Assemble dirichlet
dirichletnodes = false(mesh.info.nodcount,1);
for b_ = 1:size(problem.boundary.dirichlet,1)
  b     = problem.boundary.dirichlet{b_,1};
  form  = problem.boundary.dirichlet{b_,2};
  grade = problem.boundary.dirichlet{b_,3};
  nod   = reshape(unique(mesh.bou{b}.nod),[],1);
  coo   = mesh.nod.coo(nod,:);
  switch form
    case {[],'','0','[0,0,0]'}
      val = 0*coo;
    otherwise
      switch mesh.info.dim
        case 3
          replace_what   = {       'X',       'Y',       'Z',   'T','#'  };
          replace_to     = {'coo(:,1)','coo(:,2)','coo(:,3)','time','NaN'};
      end
      form           = replace(form,replace_what,replace_to);
      val            = eval(form);
  end
  iind = ones(mesh.info.noddof*mesh.info.nodcount,1);
  jind = 0*iind;
  jindcnt = 1;
  %% begin: not vectorized yet
  for i = 1:size(nod,1)
    switch mesh.info.noddof
      case 3
        if ~isnan(val(i,1)), Bc.U0(3*nod(i)-2,:) = val(i,1); else, iind(jindcnt+0) = 0; end
        if ~isnan(val(i,2)), Bc.U0(3*nod(i)-1,:) = val(i,2); else, iind(jindcnt+1) = 0; end
        if ~isnan(val(i,3)), Bc.U0(3*nod(i)  ,:) = val(i,3); else, iind(jindcnt+2) = 0; end
        jind(jindcnt+[0 1 2]) = 3*nod(i) - [2 1 0]; jindcnt = jindcnt + 3;
    end
  end
  %% end  : not vectorized yet
  % Bd
  ind = iind>0;
  Bc.Bd{b_} = sparse(1:sum(iind),jind(ind),1,sum(iind),mesh.info.noddof*nLf);
  add_dirichlet_dofs_logind = (sum(abs(Bd{b_}),1)>0)';
  add_dirichlet_nods_logind = add_dirichlet_dofs_logind(1:3:end) | add_dirichlet_dofs_logind(2:3:end) | add_dirichlet_dofs_logind(3:3:end);
  add_dirichlet_nods        = find(add_dirichlet_nods_logind);   dirichletnodes(add_dirichlet_nods,1) = true;
  tmp_map_logind = [add_dirichlet_dofs_logind(3*add_dirichlet_nods-2)    , add_dirichlet_dofs_logind(3*add_dirichlet_nods-1)    , add_dirichlet_dofs_logind(3*add_dirichlet_nods)    ];
  tmp_map_Bdrows = double(tmp_map_logind);             tmp_map_Bdrows(   tmp_map_logind(:,1),1) = find(sum(Bd{b_ }(:,3*add_dirichlet_nods-2),2));    tmp_map_Bdrows(   tmp_map_logind(:,2),2) = find(sum(Bd{b_ }(:,3*add_dirichlet_nods-1),2));     tmp_map_Bdrows(   tmp_map_logind(:,3),3) = find(sum(Bd{b_ }(:,3*add_dirichlet_nods),2));
  for b__ = 1:b_-1
    b__add_dirichlet_dofs_logind = (sum(abs(Bd{b__}),1)>0)';
    if nnz(b__add_dirichlet_dofs_logind & add_dirichlet_dofs_logind) > 0
      b__tmp_map_logind   = [b__add_dirichlet_dofs_logind(3*add_dirichlet_nods-2)    , b__add_dirichlet_dofs_logind(3*add_dirichlet_nods-1)    , b__add_dirichlet_dofs_logind(3*add_dirichlet_nods)    ];
      b__tmp_map_Bdrows   = double(b__tmp_map_logind);  b__tmp_map_Bdrows(b__tmp_map_logind(:,1),1) = find(sum(Bd{b__}(:,3*add_dirichlet_nods-2),2)); b__tmp_map_Bdrows(b__tmp_map_logind(:,2),2) = find(sum(Bd{b__}(:,3*add_dirichlet_nods-1),2));  b__tmp_map_Bdrows(b__tmp_map_logind(:,3),3) = find(sum(Bd{b__}(:,3*add_dirichlet_nods),2));
      tmp_delind          = b__tmp_map_logind(:,1) | b__tmp_map_logind(:,2) | b__tmp_map_logind(:,3);
      tmp_rows_to_enforce = unique(b__tmp_map_Bdrows(tmp_delind,:));   tmp_rows_to_enforce = tmp_rows_to_enforce(tmp_rows_to_enforce>0);
      tmp_rows_to_relent  = unique(   tmp_map_Bdrows(tmp_delind,:));   tmp_rows_to_relent  = tmp_rows_to_relent( tmp_rows_to_relent >0);
      tmp_rows_to_add     = tmp_map_Bdrows(tmp_delind,:);   tmp_rows_to_add = unique(tmp_rows_to_add(tmp_map_logind(tmp_delind,:) & b__tmp_map_logind(tmp_delind,:)));    tmp_rows_to_add = tmp_rows_to_add(tmp_rows_to_add>0);
      switch grade
        case 'enforce', Bd{b__}(tmp_rows_to_enforce,:) = [];
        case 'relent' , Bd{b_ }(tmp_rows_to_relent ,:) = [];
        case 'add',     Bd{b_ }(tmp_rows_to_add    ,:) = [];
        otherwise,      fprintf('tearing.m : unknown grade[%s] options are only{''enforce'',''relent'',''add''}\n',grade);
      end
    end
  end
end
Bc.U0 = U0;
Bc.Bd = cell2mat(Bd);
Bc.Bd = mat2cell(Bc.Bd,size(Bc.Bd,1),mesh.info.noddof*deco.sub2nodcount);

%% Contact
if isfield(problem.boundary,'contact')
  if isstruct(inu)
    U    = reshape(inu.u, 3,[])';
  else
    U    = reshape(inu, 3,[])';
  end
  N    = cell(size(problem.boundary.contact,1), deco.info.subcount);   T1 = N;   T2 = N;
  Dt   = cell(size(problem.boundary.contact,1), deco.info.subcount);
  SF   = cell(size(problem.boundary.contact,1), deco.info.subcount);
  gap  = cell(size(problem.boundary.contact,1), 1);
  IT   = cell(size(problem.boundary.contact,1), 1);
  ogo  = cell(size(problem.boundary.contact,1), 1);
  indA = cell(size(problem.boundary.contact,1), 1);
  globindsum = 0;
  x    = mesh.nod.coo + U;
  for ico = 1:size(problem.boundary.contact,1)
    ib     = problem.boundary.contact{ico}.bou;
    bounod = unique(mesh.bou{ib}.nod);
    logind = false(mesh.info.nodcount,1);% all global nodes
    logind(bounod) = true;
    for id_ = 1:size(problem.boundary.dirichlet)
      id = problem.boundary.dirichlet{id_,1};
      logind(unique(mesh.bou{id}.nod)) = false;
    end
    contnod           = find(logind);
    contnod_subcount  = deco.nod2subcount(logind,:);
    contnod_subcounts = unique(contnod_subcount);
    contsub           = 0*contnod;
    for icontnod_subcounts = 1:size(contnod_subcounts,1)
      tmpind = contnod_subcount==contnod_subcounts(icontnod_subcounts);
      sub = cell2mat(deco.nod2sub(contnod(tmpind,1)));
      contsub(tmpind) = sub(:,1);
    end
    switch problem.boundary.contact{ico}.type
      case 'sphere'
        sphere_centr  = eval(replace(problem.boundary.contact{ico}.param.centr,{'T','It'},{num2str(time),num2str(indtime+1)}));
        sphere_radius = str2double(problem.boundary.contact{ico}.param.rad);
        xs_x          = sphere_centr - x(logind,:);
        norm_xs_x     = sqrt(sum((xs_x).^2,2));
        normal        = xs_x./norm_xs_x;
        gap_          = norm_xs_x - sphere_radius*ones(size(normal,1),1);
      case 'cone'
        [gap_, ~, normal] = getDistanceFromCone( x(logind,:),problem.boundary.contact{ico}.param.cone );
      case 'plane'
        [gap_, ~, normal] = getDistanceFromPlane(x(logind,:),problem.boundary.contact{ico}.param.plane);
      otherwise
        disp('not programmed yet');
    end
    iind   = (1:size(normal,1))';
    for is = 1:1 %% deco.info.subcount
      gamma       = 1e-3; %7e10;
      logind      = (contsub == is);
      locind      = (1:nnz(logind))';
      globind     = deco.map_nodglo2loc{is}(contnod(logind)) + sum(deco.sub2nodcount(1:is-1,1));
      [tangent1,tangent2] = normal2tangent(normal);
      if isstruct(inu)
        inu.laCON = reshape(inu.laCON,3,[])';
        laCON     = inu.laCON(globindsum+locind,:);
      else
        laCON     = 0*normal;
      end
      globindsum   = globindsum + size(globind,1);
      laCONn       = sum(laCON.*normal,2);
      indA_        = laCONn - gamma*gap_ > 0;
      ogo_         = 0*normal;
      ogo_(indA_,1)= 1*gap_(indA_);
      s_           = -normal; % -xs_x./norm_xs_x
      f1_          = 0*s_;
      f2_          = 0*s_;
      valIT        = zeros(size(normal,1),9);
      valSF        = zeros(size(normal,1),9);
      if nnz(  indA_)>0
        valIT( indA_,4:6) = tangent1( indA_,:);        % T_A1
        valIT( indA_,7:9) = tangent2( indA_,:);        % T_A2
        valSF( indA_,1:3) = s_(       indA_,:);        % S
        valSF( indA_,4:6) = f1_(      indA_,:);        % F1
        valSF( indA_,7:9) = f2_(      indA_,:);        % F2
      end
      if nnz( ~indA_)>0
        valIT(~indA_,[1 5 9]) = 1;              % I_A
      end
      N{   ico,is} = sparse( (iind(logind))*[1 1 1]        , 3*deco.map_nodglo2loc{is}(contnod(logind))*[1 1 1]-[2 1 0], normal(  logind,:),   size(normal,1), 3*deco.sub2nodcount(is));
      T1{  ico,is} = sparse( (iind(logind))*[1 1 1]        , 3*deco.map_nodglo2loc{is}(contnod(logind))*[1 1 1]-[2 1 0], tangent1(logind,:),   size(normal,1), 3*deco.sub2nodcount(is));
      T2{  ico,is} = sparse( (iind(logind))*[1 1 1]        , 3*deco.map_nodglo2loc{is}(contnod(logind))*[1 1 1]-[2 1 0], tangent2(logind,:),   size(normal,1), 3*deco.sub2nodcount(is));
      gap{ ico}    = sparse(  iind(logind)                 , 1                                                         , gap_              ,   size(normal,1), 1                      );
      Dt{  ico,is} = sparse(3*iind(logind)* [1 1 1]-[2 1 0], 3*deco.map_nodglo2loc{is}(contnod(logind))*[1 1 1]-[2 1 0], ones(size(normal)), 3*size(normal,1), 3*deco.sub2nodcount(is));
      IT{  ico}    = sparse(3*locind-[2 2 2 1 1 1 0 0 0]   , 3*locind -[2 1 0 2 1 0 2 1 0]                             , valIT             , 3*size(normal,1), 3*size(normal,1)       );
      ogo{ ico}    = ogo_;
      SF{  ico,is} = sparse(3*locind-[2 2 2 1 1 1 0 0 0]   , 3*globind-[2 1 0 2 1 0 2 1 0]                             , valSF             , 3*size(normal,1), 3*deco.sub2nodcount(is));
      indA{ico}    = indA_;
    end
  end
  d           = mesh.info.noddof;
  tmp         = size(cell2mat(N),1);
  Tearing.N   = cell2mat(N );   Tearing.N  = mat2cell(Tearing.N ,  tmp,d*deco.sub2nodcount);
  Tearing.Dt  = cell2mat(Dt);   Tearing.Dt = mat2cell(Tearing.Dt,d*tmp,d*deco.sub2nodcount);
  Tearing.T1  = cell2mat(T1);   Tearing.T1 = mat2cell(Tearing.T1,  tmp,d*deco.sub2nodcount);
  Tearing.T2  = cell2mat(T2);   Tearing.T2 = mat2cell(Tearing.T2,  tmp,d*deco.sub2nodcount);
  Tearing.SF  = cell2mat(SF);   Tearing.SF = mat2cell(Tearing.SF,d*tmp,d*deco.sub2nodcount);
  Tearing.IT  = blkdiag(IT{:}); Tearing.IT = mat2cell(Tearing.IT,d*tmp,d*tmp              );
  Tearing.iA  = indA;
  Tearing.gap = gap;
  Tearing.ogo = ogo;
  Tearing.sphere_disp  = sphere_centr-eval(replace(problem.boundary.contact{1}.param.centr,{'T','It'},{'0','1'}));

else
  Tearing.N   = cell(1,deco.info.subcount); for i = 1:deco.info.subcount, Tearing.N{1,i} = zeros(0,mesh.info.noddof*deco.sub2nodcount(i)); end
  Tearing.gap = zeros(0,1);
end

end
