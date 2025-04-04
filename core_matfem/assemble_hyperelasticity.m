function [Discretization,Byproducts] = assemble_hyperelasticity(FE,mesh,consts,problem,u,time)

  function y = normcross32k(x)
    % |a x b|^2 = |a|^2 * |b|^2 -(a.b)^2
    y = sqrt( (x(1,1,:).^2+x(2,1,:).^2+x(3,1,:).^2).*(x(1,2,:).^2+x(2,2,:).^2+x(3,2,:).^2) - (x(1,1,:).*x(1,2,:) + x(2,1,:).*x(2,2,:) + x(3,1,:).*x(3,2,:)).^2 );
  end

Discretization = struct('M' ,[],'Ks',[],'Kc',[],'fext',[], 'fint', []);

Byproducts = struct(...
  'nod',struct('u',u),...
  'ele',struct('sigma',zeros(mesh.info.elecount,6)));


%%% DISCRETIZATION %%%
%fprintf('subdomain %03d',s);
% volume
Discretization.M      = sparse(mesh.info.noddof*mesh.info.nodcount);
Discretization.Kc     = sparse(mesh.info.noddof*mesh.info.nodcount);
Discretization.Ks     = sparse(mesh.info.noddof*mesh.info.nodcount);
Discretization.fint   = zeros( mesh.info.noddof*mesh.info.nodcount,1);
Discretization.fext   = zeros( mesh.info.noddof*mesh.info.nodcount,1);
n_e       = size(mesh.ele.nod,1);
n_loc     = size(mesh.ele.nod,2);
veclocdofcount = mesh.info.noddof*n_loc;
matlocdofcount = veclocdofcount^2;
iii = zeros(n_e*(mesh.info.noddof*n_loc)^2,1);
jjj = iii;   vvKc = iii;   vvKs = iii;   vvM = iii;

for e = 1:n_e
  elnod   = mesh.ele.nod(e,:);
  Xloc    = mesh.nod.coo(elnod,:)';
  uloc    = [u(3*elnod-2,:)'; u(3*elnod-1,:)'; u(3*elnod,:)'];
  ivecloc = reshape(mesh.info.noddof*(ones(mesh.info.noddof,1)*elnod-1)+(1:mesh.info.noddof)',[],1);
  imatloc = ivecloc*ones(1,mesh.info.noddof*size(elnod,2));
  [Mloc,Kcloc,Ksloc,fextloc,fintloc,bypr] = assemble_hyperelasticity_loc(Xloc,uloc,FE,consts,problem);
  Byproducts.ele.sigma(e,:) = bypr.S;
  iii( (e-1)*matlocdofcount+1:e*matlocdofcount)   = reshape(imatloc ,matlocdofcount,1);
  jjj( (e-1)*matlocdofcount+1:e*matlocdofcount)   = reshape(imatloc',matlocdofcount,1);
  vvKc((e-1)*matlocdofcount+1:e*matlocdofcount)   = reshape(   Kcloc,matlocdofcount,1);
  vvKs((e-1)*matlocdofcount+1:e*matlocdofcount)   = reshape(   Ksloc,matlocdofcount,1);
  vvM( (e-1)*matlocdofcount+1:e*matlocdofcount)   = reshape(    Mloc,matlocdofcount,1);
  Discretization.fext(ivecloc) = Discretization.fext(ivecloc) + fextloc;
  Discretization.fint(ivecloc) = Discretization.fint(ivecloc) + fintloc;
end
%fprintf(' ... sparse ... ');
Discretization.M  = sparse(iii,jjj,vvM );
Discretization.Kc = sparse(iii,jjj,vvKc);
Discretization.Ks = sparse(iii,jjj,vvKs);
%fprintf(' ... done\n');

%{
%--------------------------------------------------------------
d   = mesh.info.dim;
dm1 = mesh.info.dim-1;
for b_ = 1:size(problem.boundary.neumann,1)
  b     = problem.boundary.neumann{b_,1};
  form  = problem.boundary.neumann{b_,2};
  nod   = reshape(unique(mesh.bou{b}.nod),[],1);
  nod2sub      = deco.nod2sub(     nod,:);
  nod2subcount = deco.nod2subcount(nod,:);
  coo   = mesh.nod.coo(nod,:);
  switch form
    case {[],'','0','[0,0,0]'}
      val = 0*coo;
    otherwise
      vars  = problem.variables;
      vars_tmpnames = cell(size(vars,1),1);
      switch mesh.info.dim
        case 2
          for i = 1:size(vars,1)
            replace_what = {       'X',       'Y',   'T','#'  ,problem.constants{:,1},problem.variables{1:i-1,1}};
            replace_to   = {'coo(:,1)','coo(:,2)','time','NaN',problem.constants{:,2},vars_tmpnames{    1:i-1,1}};
            vars{i,2}    = eval(replace(vars{i,2},replace_what,replace_to));
            vars_tmpnames{i} = ['vars{',num2str(i),',2}'];
          end
          replace_what   = {       'X',       'Y',   'T','#'  ,problem.constants{:,1},problem.variables{:,1}};
          replace_to     = {'coo(:,1)','coo(:,2)','time','NaN',problem.constants{:,2},vars_tmpnames{    :,1}};
        case 3
          for i = 1:size(vars,1)
            replace_what = {       'X',       'Y',       'Z',   'T','#'  ,problem.constants{:,1},problem.variables{1:i-1,1}};
            replace_to   = {'coo(:,1)','coo(:,2)','coo(:,3)','time','NaN',problem.constants{:,2},vars_tmpnames{    1:i-1,1}};
            vars{i,2}    = eval(replace(vars{i,2},replace_what,replace_to));
            vars_tmpnames{i} = ['vars{',num2str(i),',2}'];
          end
          replace_what   = {       'X',       'Y',       'Z',   'T','#'  ,problem.constants{:,1},problem.variables{:,1}};
          replace_to     = {'coo(:,1)','coo(:,2)','coo(:,3)','time','NaN',problem.constants{:,2},vars_tmpnames{    :,1}};
      end
      form           = replace(['0*X+',form],replace_what,replace_to);
      val            = eval(form);
  end
  switch mesh.info.dim
    case 2
      pressure_map = sparse([nod;nod    ],[0*nod+1;0*nod+2        ],[val(:,1);val(:,2)         ]);
    case 3
      pressure_map = sparse([nod;nod;nod],[0*nod+1;0*nod+2;0*nod+3],[val(:,1);val(:,2);val(:,3)]);
  end
  FEbou = fe_init(mesh.bou{b}.type,6);
  N__loc    = FEbou.N(FEbou.QP);
  dNdxi_loc = FEbou.dNdxi(FEbou.QP);
  %dimbou = size(FEbou.QP,1);
  nqpbou = size(FEbou.QP,2);
  nNbou  = size(FEbou.xi,1);
  dNdxi_loc = permute(reshape(dNdxi_loc,nNbou,nqpbou,dm1),[1 3 2]);               % N . 2 . q
  for e_ = 1:size(mesh.bou{b}.ele,1)
    elnod = mesh.bou{b}.nod(e_,:);
    e = mesh.bou{b}.ele(e_,:);
    s = deco.ele2sub(e);
    Xloc = mesh.nod.coo(elnod,:)';
    ivecloc = reshape(mesh.info.noddof*(ones(mesh.info.noddof,1)*elnod-1)+(1:mesh.info.noddof)',[],1);
    dXdxi_loc = reshape(Xloc*reshape(dNdxi_loc,nNbou,dm1*nqpbou),d,dm1,nqpbou);  % 3 . 2 . q
    switch mesh.info.dim
      case 3
        detdXdxi_loc = reshape(normcross32k(dXdxi_loc),1,nqpbou);                % q
      case 2
        detdXdxi_loc = reshape(norm21k(dXdxi_loc),1,nqpbou);                     % q
    end
    absJW = abs(detdXdxi_loc)'.*FEbou.QW';                                       % q
    floc = reshape(((N__loc*absJW).*pressure_map(elnod,:))',[],1);
    if nnz(deco.map_dofglo2loc{s}(ivecloc)) > 0
      Discretization.f{s}(deco.map_dofglo2loc{s}(ivecloc)) = Discretization.f{s}(deco.map_dofglo2loc{s}(ivecloc)) + floc;
    end
  end
end
%}
end
