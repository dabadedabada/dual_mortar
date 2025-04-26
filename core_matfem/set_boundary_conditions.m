function [bc] = set_boundary_conditions(mesh,problem,time,indtime,inu)

bc.Bd  = cell(size(problem.boundary.dirichlet,1),1);
bc.U0  = zeros(mesh.info.noddof*mesh.info.nodcount,1);

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
        if ~isnan(val(i,1)), bc.U0(3*nod(i)-2,:) = val(i,1); else, iind(jindcnt+0) = 0; end
        if ~isnan(val(i,2)), bc.U0(3*nod(i)-1,:) = val(i,2); else, iind(jindcnt+1) = 0; end
        if ~isnan(val(i,3)), bc.U0(3*nod(i)  ,:) = val(i,3); else, iind(jindcnt+2) = 0; end
        jind(jindcnt+[0 1 2]) = 3*nod(i) - [2 1 0]; jindcnt = jindcnt + 3;
    end
  end
  %% end  : not vectorized yet
  % Bd
  ind = jind>0;
  bc.Bd{b_} = sparse(1:sum(ind),jind(ind),1,sum(ind),mesh.info.noddof*mesh.info.nodcount);
  add_dirichlet_dofs_logind = (sum(abs(bc.Bd{b_}),1)>0)';
  add_dirichlet_nods_logind = add_dirichlet_dofs_logind(1:3:end) | add_dirichlet_dofs_logind(2:3:end) | add_dirichlet_dofs_logind(3:3:end);
  add_dirichlet_nods        = find(add_dirichlet_nods_logind);   dirichletnodes(add_dirichlet_nods,1) = true;
  tmp_map_logind = [add_dirichlet_dofs_logind(3*add_dirichlet_nods-2)    , add_dirichlet_dofs_logind(3*add_dirichlet_nods-1)    , add_dirichlet_dofs_logind(3*add_dirichlet_nods)    ];
  tmp_map_Bdrows = double(tmp_map_logind);             tmp_map_Bdrows(   tmp_map_logind(:,1),1) = find(sum(bc.Bd{b_ }(:,3*add_dirichlet_nods-2),2));    tmp_map_Bdrows(   tmp_map_logind(:,2),2) = find(sum(bc.Bd{b_ }(:,3*add_dirichlet_nods-1),2));     tmp_map_Bdrows(   tmp_map_logind(:,3),3) = find(sum(bc.Bd{b_ }(:,3*add_dirichlet_nods),2));
  for b__ = 1:b_-1
    b__add_dirichlet_dofs_logind = (sum(abs(bc.Bd{b__}),1)>0)';
    if nnz(b__add_dirichlet_dofs_logind & add_dirichlet_dofs_logind) > 0
      b__tmp_map_logind   = [b__add_dirichlet_dofs_logind(3*add_dirichlet_nods-2)    , b__add_dirichlet_dofs_logind(3*add_dirichlet_nods-1)    , b__add_dirichlet_dofs_logind(3*add_dirichlet_nods)    ];
      b__tmp_map_Bdrows   = double(b__tmp_map_logind);  b__tmp_map_Bdrows(b__tmp_map_logind(:,1),1) = find(sum(bc.Bd{b__}(:,3*add_dirichlet_nods-2),2)); b__tmp_map_Bdrows(b__tmp_map_logind(:,2),2) = find(sum(bc.Bd{b__}(:,3*add_dirichlet_nods-1),2));  b__tmp_map_Bdrows(b__tmp_map_logind(:,3),3) = find(sum(bc.Bd{b__}(:,3*add_dirichlet_nods),2));
      tmp_delind          = b__tmp_map_logind(:,1) | b__tmp_map_logind(:,2) | b__tmp_map_logind(:,3);
      tmp_rows_to_enforce = unique(b__tmp_map_Bdrows(tmp_delind,:));   tmp_rows_to_enforce = tmp_rows_to_enforce(tmp_rows_to_enforce>0);
      tmp_rows_to_relent  = unique(   tmp_map_Bdrows(tmp_delind,:));   tmp_rows_to_relent  = tmp_rows_to_relent( tmp_rows_to_relent >0);
      tmp_rows_to_add     = tmp_map_Bdrows(tmp_delind,:);   tmp_rows_to_add = unique(tmp_rows_to_add(tmp_map_logind(tmp_delind,:) & b__tmp_map_logind(tmp_delind,:)));    tmp_rows_to_add = tmp_rows_to_add(tmp_rows_to_add>0);
      switch grade
        case 'enforce', bc.Bd{b__}(tmp_rows_to_enforce,:) = [];
        case 'relent' , bc.Bd{b_ }(tmp_rows_to_relent ,:) = [];
        case 'add'    , % todo    
        otherwise,      fprintf('tearing.m : unknown grade[%s] options are only{''enforce'',''relent'',''add''}\n',grade);
      end
    end
  end
end
bc.Bd = cell2mat(bc.Bd);
bc.Bd = mat2cell(bc.Bd,size(bc.Bd,1),mesh.info.noddof*mesh.info.nodcount);

end
