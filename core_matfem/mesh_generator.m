function mesh = mesh_generator(L,ne,etype)
% input:
%   L         ... lenghts                x  y [z]
%   ne        ...   elecountments/subdomain  nx ny [nz]
%   etype     ... elecountment type {tria3,tria6,quad4,quad8,tetr4,tetr10,hexa8,hexa20,hexa27}
% output:
%   mesh.inf  ... info                struct {noddof,nodcount,nelecount}
%       .nod  ... nodes               struct {coo}
%       .elecount  ... elecountments            struct {nod,type}
%       .bou  ... boundary elecountments   cell of struct {nod,type}
%
%  element nodes numbering see ../doc/fe.lyx

switch etype
  %                dim      #local_nodes   boundary_type      #face_nodes #boundaries 
  case 'tria3' ,   d = 2;   nnl =  3;      btype = 'line2';   nbn = 2;     nb = 4;
  case 'tria6' ,   d = 2;   nnl =  6;      btype = 'line3';   nbn = 3;     nb = 4;
  case 'quad4' ,   d = 2;   nnl =  4;      btype = 'line2';   nbn = 2;     nb = 4;
  case 'quad8' ,   d = 2;   nnl =  8;      btype = 'line3';   nbn = 3;     nb = 4;
  case 'quad9' ,   d = 2;   nnl =  9;      btype = 'line3';   nbn = 3;     nb = 4;
  case 'tetr4' ,   d = 3;   nnl =  4;      btype = 'tria3';   nbn = 3;     nb = 6;
  case 'tetr10',   d = 3;   nnl = 10;      btype = 'tria6';   nbn = 3;     nb = 6;
  case 'hexa8' ,   d = 3;   nnl =  8;      btype = 'quad4';   nbn = 4;     nb = 6;
  case 'hexa20',   d = 3;   nnl = 20;      btype = 'quad8';   nbn = 8;     nb = 6;
  case 'hexa27',   d = 3;   nnl = 27;      btype = 'quad9';   nbn = 9;     nb = 6;
  otherwise, error(['Unknown elecountment type: ',etype]);
end
if size( L,2) == 1, L  =  L*ones(1,d); end;   if d == 2,  L = [ L 1]; end
if size(ne,2) == 1, ne = ne*ones(1,d); end;   if d == 2, ne = [ne 1]; end
nn  = prod(ne+1);                                                          % number of nodes
nbe = [ne(2)*ne(3)*[1;1];  ne(1)*ne(3)*[1;1];  (d-2)*ne(1)*ne(2)*[1;1]];   % number of boundary elecountments

mesh = struct(...
  'info',struct('noddof',[],'nodcount',[],'elecount',[],'boucounts',nb,'boucount',[],'extent',zeros(d,2)),...
  'nod' ,struct('coo',[]),...
  'ele' ,struct('nod',[],'type',etype),...
  'bou' ,[]);   mesh.bou = cell(nb,1);
for i = 1:nb, mesh.bou{i} = struct('nod',[],'ele',[],'type',btype); end
switch etype
  case {'quad4','tria3'}
    [xx,yy] = meshgrid( 0:L(1)/ne(1):L(1), 0:L(2)/ne(2):L(2));
    [ii,jj] = meshgrid( 0:ne(1)          , 0:ne(2)          );
    mesh.nod.coo = [reshape(permute(xx,[2 1]),[],1) reshape(permute(yy,[2 1]),[],1)];
    nod = reshape(1:size(mesh.nod.coo,1),ne(1)+1,ne(2)+1,1);
    pne = prod(ne);
    ele = reshape(1:pne,ne(1),ne(2));
    elm = mod(ii(1:end-1,1:end-1)+jj(1:end-1,1:end-1),2)==0; % element mirrored
    % elements
    % n1 n4
    % n2 n3
    n1 = reshape(nod(1:end-1,1:end-1),[],1);
    n2 = reshape(nod(2:end  ,1:end-1),[],1);
    n3 = reshape(nod(2:end  ,2:end  ),[],1);
    n4 = reshape(nod(1:end-1,2:end  ),[],1);
    switch etype
      case 'quad4'
        mesh.ele.nod = [n1 n2 n3 n4];
      case 'tria3'
        mesh.ele.nod = zeros(size(n1,1),6);
        ind = reshape( elm,[],1);
        mesh.ele.nod(ind,:) = [[n1(ind,:) n2(ind,:) n3(ind,:)] [n1(ind,:) n3(ind,:) n4(ind,:)]];
        ind = reshape(~elm,[],1);
        mesh.ele.nod(ind,:) = [[n1(ind,:) n2(ind,:) n4(ind,:)] [n2(ind,:) n3(ind,:) n4(ind,:)]];
        mesh.ele.nod = [mesh.ele.nod(:,1:3); mesh.ele.nod(:,4:6)];
    end
    % boundary elements
    switch etype
      case 'quad4'
        mesh.bou{1}.nod = [reshape(nod(      1,2:end  ),[],1) reshape(nod(      1,1:end-1),[],1)];   mesh.bou{1}.ele = reshape(ele(    1,1:end),[],1); %  1 (x == 0)
        mesh.bou{2}.nod = [reshape(nod(    end,1:end-1),[],1) reshape(nod(    end,2:end  ),[],1)];   mesh.bou{2}.ele = reshape(ele(  end,1:end),[],1); %  2 (x == L)
        mesh.bou{3}.nod = [reshape(nod(1:end-1,      1),[],1) reshape(nod(2:end  ,      1),[],1)];   mesh.bou{3}.ele = reshape(ele(1:end,    1),[],1); %  3 (y == 0)
        mesh.bou{4}.nod = [reshape(nod(2:end  ,    end),[],1) reshape(nod(1:end-1,    end),[],1)];   mesh.bou{4}.ele = reshape(ele(1:end,  end),[],1); %  4 (y == L)
      case 'tria3'
        mesh.bou{1}.nod = [reshape(nod(      1,2:end  ),[],1) reshape(nod(      1,1:end-1),[],1)];   mesh.bou{1}.ele = reshape(ele(    1,1:end),[],1)+reshape(elm(    1,1:end),[],1)*pne; %  1 (x == 0)
        mesh.bou{2}.nod = [reshape(nod(    end,1:end-1),[],1) reshape(nod(    end,2:end  ),[],1)];   mesh.bou{2}.ele = reshape(ele(  end,1:end),[],1)+reshape(elm(  end,1:end),[],1)*pne; %  2 (x == L)
        mesh.bou{3}.nod = [reshape(nod(1:end-1,      1),[],1) reshape(nod(2:end  ,      1),[],1)];   mesh.bou{3}.ele = reshape(ele(1:end,    1),[],1)+reshape(elm(1:end,    1),[],1)*pne; %  3 (y == 0)
        mesh.bou{4}.nod = [reshape(nod(2:end  ,    end),[],1) reshape(nod(1:end-1,    end),[],1)];   mesh.bou{4}.ele = reshape(ele(1:end,  end),[],1)+reshape(elm(1:end,  end),[],1)*pne; %  4 (y == L)        
    end
  case {'tria6','quad8','quad9'}
    [xx,yy] = meshgrid( 0:L(1)/(2*ne(1)):L(1), 0:L(2)/(2*ne(2)):L(2));
    [ii,jj] = meshgrid( 0:ne(1)              , 0:ne(2)              );
    mesh.nod.coo = [reshape(permute(xx,[2 1]),[],1) reshape(permute(yy,[2 1]),[],1)];
    delcoo = false(size(mesh.nod.coo,1),1);
    nod = reshape(1:size(mesh.nod.coo,1),2*ne(1)+1,2*ne(2)+1,1);
    pne = prod(ne);
    ele = reshape(1:pne,ne(1),ne(2));
    elm = mod(ii(1:end-1,1:end-1)+jj(1:end-1,1:end-1),2)==0; % element mirrored
    % elements
    %  n1 n8 n4
    %  n5 n9 n7
    %  n2 n6 n3
    n1 = reshape(nod(1:2:end-2,1:2:end-2),[],1);
    n2 = reshape(nod(1:2:end-2,3:2:end  ),[],1);
    n3 = reshape(nod(3:2:end  ,3:2:end  ),[],1);
    n4 = reshape(nod(3:2:end  ,1:2:end-2),[],1);
    n5 = reshape(nod(1:2:end-2,2:2:end-1),[],1);
    n6 = reshape(nod(2:2:end-1,3:2:end  ),[],1);
    n7 = reshape(nod(3:2:end  ,2:2:end-1),[],1);
    n8 = reshape(nod(2:2:end-1,1:2:end-2),[],1);
    n9 = reshape(nod(2:2:end-1,2:2:end-1),[],1);
    switch etype
      case 'quad8'
        mesh.ele.nod = [n1 n2 n3 n4 n5 n6 n7 n8];
        delcoo(n9) = true;
      case 'quad9'
        mesh.ele.nod = [n1 n2 n3 n4 n5 n6 n7 n8 n9];
      case 'tria6'
        mesh.ele.nod = zeros(size(n1,1),12);
        ind = reshape(mod(ii(1:end-1,1:end-1)+jj(1:end-1,1:end-1),2),[],1)==0;
        mesh.ele.nod(ind,:) = [[n1(ind,:) n2(ind,:) n3(ind,:) n5(ind,:) n6(ind,:) n9(ind,:)] [n1(ind,:) n3(ind,:) n4(ind,:) n9(ind,:) n7(ind,:) n8(ind,:)]];
        ind = reshape(mod(ii(1:end-1,1:end-1)+jj(1:end-1,1:end-1),2),[],1)==1;
        mesh.ele.nod(ind,:) = [[n1(ind,:) n2(ind,:) n4(ind,:) n5(ind,:) n9(ind,:) n8(ind,:)] [n2(ind,:) n3(ind,:) n4(ind,:) n6(ind,:) n7(ind,:) n9(ind,:)]];
        mesh.ele.nod = [mesh.ele.nod(:,1:6); mesh.ele.nod(:,7:12)];
    end
    % boundary elements
    switch etype
      case {'quad8','quad9'}
        mesh.bou{1}.nod = [reshape(nod(        1,3:2:end  ),[],1) reshape(nod(        1,1:2:end-2),[],1) reshape(nod(        1,2:2:end-1),[],1)];   mesh.bou{1}.ele = reshape(ele(    1,1:end),[],1); %  1 (x == 0)
        mesh.bou{2}.nod = [reshape(nod(      end,1:2:end-2),[],1) reshape(nod(      end,3:2:end  ),[],1) reshape(nod(      end,2:2:end-1),[],1)];   mesh.bou{2}.ele = reshape(ele(  end,1:end),[],1); %  2 (x == L)
        mesh.bou{3}.nod = [reshape(nod(1:2:end-2,        1),[],1) reshape(nod(3:2:end  ,        1),[],1) reshape(nod(2:2:end-1,        1),[],1)];   mesh.bou{3}.ele = reshape(ele(1:end,    1),[],1); %  3 (y == 0)
        mesh.bou{4}.nod = [reshape(nod(3:2:end  ,      end),[],1) reshape(nod(1:2:end-2,      end),[],1) reshape(nod(2:2:end-1,      end),[],1)];   mesh.bou{4}.ele = reshape(ele(1:end,  end),[],1); %  4 (y == L)
      case 'tria6'
        mesh.bou{1}.nod = [reshape(nod(        1,3:2:end  ),[],1) reshape(nod(        1,1:2:end-2),[],1) reshape(nod(        1,2:2:end-1),[],1)];   mesh.bou{1}.ele = reshape(ele(    1,1:end),[],1)+reshape(elm(    1,1:end),[],1)*pne; %  1 (x == 0)
        mesh.bou{2}.nod = [reshape(nod(      end,1:2:end-2),[],1) reshape(nod(      end,3:2:end  ),[],1) reshape(nod(      end,2:2:end-1),[],1)];   mesh.bou{2}.ele = reshape(ele(  end,1:end),[],1)+reshape(elm(  end,1:end),[],1)*pne; %  2 (x == L)
        mesh.bou{3}.nod = [reshape(nod(1:2:end-2,        1),[],1) reshape(nod(3:2:end  ,        1),[],1) reshape(nod(2:2:end-1,        1),[],1)];   mesh.bou{3}.ele = reshape(ele(1:end,    1),[],1)+reshape(elm(1:end,    1),[],1)*pne; %  3 (y == 0)
        mesh.bou{4}.nod = [reshape(nod(3:2:end  ,      end),[],1) reshape(nod(1:2:end-2,      end),[],1) reshape(nod(2:2:end-1,      end),[],1)];   mesh.bou{4}.ele = reshape(ele(1:end,  end),[],1)+reshape(elm(1:end,  end),[],1)*pne; %  4 (y == L)        
    end
    % delete redundant coordinates
    switch etype
      case 'quad8'
        map = zeros(size(delcoo,1),1);   map(~delcoo,1) = 1:sum(~delcoo);
        mesh.nod.coo(delcoo,:) = [];
        mesh.ele.nod    = map(   mesh.ele.nod);
        mesh.bou{1}.nod = map(mesh.bou{1}.nod);
        mesh.bou{2}.nod = map(mesh.bou{2}.nod);
        mesh.bou{3}.nod = map(mesh.bou{3}.nod);
        mesh.bou{4}.nod = map(mesh.bou{4}.nod);
    end
  case {'hexa8','tetr4'}
    [xx,yy,zz] = meshgrid( 0:L(1)/ne(1):L(1), 0:L(2)/ne(2):L(2), 0:L(3)/ne(3):L(3));   xx = permute(xx,[2 1 3]);   yy = permute(yy,[2 1 3]);   zz = permute(zz,[2 1 3]);
    [ii,jj,kk] = meshgrid( 0:ne(1)          , 0:ne(2)          , 0:ne(3)          );   ii = permute(ii,[2 1 3]);   jj = permute(jj,[2 1 3]);   kk = permute(kk,[2 1 3]);
    mesh.nod.coo = [reshape(xx,[],1) reshape(yy,[],1) reshape(zz,[],1)];
    delcoo = false(size(mesh.nod.coo,1),1);
    nod = reshape(1:size(mesh.nod.coo,1),ne(1)+1,ne(2)+1,ne(3)+1);
    % elements
    %  z=0     z=1
    % n1 n4 | n5 n8
    % n2 n3 | n6 n7
    pne = prod(ne);
    ele = reshape(1:pne,ne(1),ne(2),ne(3));
    elm = mod(ii(1:end-1,1:end-1,1:end-1)+jj(1:end-1,1:end-1,1:end-1)+kk(1:end-1,1:end-1,1:end-1),2)==0; % element mirrored
    n1 = reshape(nod(1:end-1,1:end-1,1:end-1),[],1);
    n2 = reshape(nod(2:end  ,1:end-1,1:end-1),[],1);
    n3 = reshape(nod(2:end  ,2:end  ,1:end-1),[],1);
    n4 = reshape(nod(1:end-1,2:end  ,1:end-1),[],1);
    n5 = reshape(nod(1:end-1,1:end-1,2:end  ),[],1);
    n6 = reshape(nod(2:end  ,1:end-1,2:end  ),[],1);
    n7 = reshape(nod(2:end  ,2:end  ,2:end  ),[],1);
    n8 = reshape(nod(1:end-1,2:end  ,2:end  ),[],1);
    switch etype
      case 'hexa8'
        mesh.ele.nod = [n1 n2 n3 n4 n5 n6 n7 n8];
      case 'tetr4'
        tmp = zeros(size(n1,1),20);
        ind = reshape(mod(ii(1:end-1,1:end-1,1:end-1)+jj(1:end-1,1:end-1,1:end-1)+kk(1:end-1,1:end-1,1:end-1),2),[],1)==0;
        tmp(ind,:) = [[n1(ind,:) n2(ind,:) n3(ind,:) n6(ind,:)] [n1(ind,:) n3(ind,:) n4(ind,:) n8(ind,:)] [n1(ind,:) n5(ind,:) n6(ind,:) n8(ind,:)] [n3(ind,:) n6(ind,:) n7(ind,:) n8(ind,:)] [n1(ind,:) n6(ind,:) n3(ind,:) n8(ind,:)]];
        ind = ~ind; %ind = reshape(mod(ii(1:end-1,1:end-1,1:end-1)+jj(1:end-1,1:end-1,1:end-1)+kk(1:end-1,1:end-1,1:end-1),2),[],1)==1;
        tmp(ind,:) = [[n4(ind,:) n1(ind,:) n2(ind,:) n5(ind,:)] [n2(ind,:) n3(ind,:) n4(ind,:) n7(ind,:)] [n5(ind,:) n6(ind,:) n2(ind,:) n7(ind,:)] [n5(ind,:) n8(ind,:) n7(ind,:) n4(ind,:)] [n2(ind,:) n5(ind,:) n7(ind,:) n4(ind,:)]];
        mesh.ele.nod = zeros(5*size(n1,1),4);
        mesh.ele.nod(1:5:end,:) = tmp(:, 1: 4);
        mesh.ele.nod(2:5:end,:) = tmp(:, 5: 8);
        mesh.ele.nod(3:5:end,:) = tmp(:, 9:12);
        mesh.ele.nod(4:5:end,:) = tmp(:,13:16);
        mesh.ele.nod(5:5:end,:) = tmp(:,17:20);
    end
    % boundary elements
    switch etype
      case 'hexa8'
        n1 = reshape(nod(1      ,1:end-1,1:end-1),[],1); n4 = reshape(nod(1      ,2:end  ,1:end-1),[],1); n8 = reshape(nod(1      ,2:end  ,2:end  ),[],1); n5 = reshape(nod(1      ,1:end-1,2:end  ),[],1);   mesh.bou{1}.nod = [n1 n4 n8 n5];   mesh.bou{1}.ele = reshape(ele(1    ,1:end,1:end),[],1);   %  1 (x == 0)
        n2 = reshape(nod(  end  ,1:end-1,1:end-1),[],1); n3 = reshape(nod(  end  ,2:end  ,1:end-1),[],1); n7 = reshape(nod(  end  ,2:end  ,2:end  ),[],1); n6 = reshape(nod(  end  ,1:end-1,2:end  ),[],1);   mesh.bou{2}.nod = [n6 n7 n3 n2];   mesh.bou{2}.ele = reshape(ele(  end,1:end,1:end),[],1);   %  2 (x == L)
        n1 = reshape(nod(1:end-1,1      ,1:end-1),[],1); n2 = reshape(nod(2:end  ,1      ,1:end-1),[],1); n6 = reshape(nod(2:end  ,1      ,2:end  ),[],1); n5 = reshape(nod(1:end-1,1      ,2:end  ),[],1);   mesh.bou{3}.nod = [n1 n2 n6 n5];   mesh.bou{3}.ele = reshape(ele(1:end,1    ,1:end),[],1);   %  3 (y == 0)
        n4 = reshape(nod(1:end-1,  end  ,1:end-1),[],1); n3 = reshape(nod(2:end  ,  end  ,1:end-1),[],1); n7 = reshape(nod(2:end  ,  end  ,2:end  ),[],1); n8 = reshape(nod(1:end-1,  end  ,2:end  ),[],1);   mesh.bou{4}.nod = [n8 n7 n3 n4];   mesh.bou{4}.ele = reshape(ele(1:end,  end,1:end),[],1);   %  4 (y == L)
        n1 = reshape(nod(1:end-1,1:end-1,1      ),[],1); n2 = reshape(nod(2:end  ,1:end-1,1      ),[],1); n3 = reshape(nod(2:end  ,2:end  ,1      ),[],1); n4 = reshape(nod(1:end-1,2:end  ,1      ),[],1);   mesh.bou{5}.nod = [n1 n2 n3 n4];   mesh.bou{5}.ele = reshape(ele(1:end,1:end,1    ),[],1);   %  5 (z == 0)
        n5 = reshape(nod(1:end-1,1:end-1,  end  ),[],1); n6 = reshape(nod(2:end  ,1:end-1,  end  ),[],1); n7 = reshape(nod(2:end  ,2:end  ,  end  ),[],1); n8 = reshape(nod(1:end-1,2:end  ,  end  ),[],1);   mesh.bou{6}.nod = [n8 n7 n6 n5];   mesh.bou{6}.ele = reshape(ele(1:end,1:end,  end),[],1);   %  6 (z == L)
      case 'tetr4'
        eo = mod(ii(1:end-1,1:end-1,1:end-1)+jj(1:end-1,1:end-1,1:end-1)+kk(1:end-1,1:end-1,1:end-1),2);
        n1 = reshape(nod(1      ,1:end-1,1:end-1),[],1); n4 = reshape(nod(1      ,2:end  ,1:end-1),[],1); n8 = reshape(nod(1      ,2:end  ,2:end  ),[],1); n5 = reshape(nod(1      ,1:end-1,2:end  ),[],1);   ind = reshape(elm(  1,  :,  :),[],1);   tmp = zeros(size(n1,1),6);   tmp( ind,:) = [[n4( ind,:) n8( ind,:) n1( ind,:)] [n5( ind,:) n1( ind,:) n8( ind,:)]];   tmp(~ind,:) = [[n1(~ind,:) n4(~ind,:) n5(~ind,:)] [n8(~ind,:) n5(~ind,:) n4(~ind,:)]];   mesh.bou{1}.nod = [tmp(:,1:3); tmp(:,4:6)];   mesh.bou{1}.ele = [reshape((ele(  1,  :,  :)-1)*5+2*eo(  1,  :,  :)+1*~eo(  1,  :,  :),[],1); reshape((ele(  1,  :,  :)-1)*5+3*eo(  1,  :,  :)+4*~eo(  1,  :,  :),[],1)];   %  1 (x == 0)
        n2 = reshape(nod(  end  ,1:end-1,1:end-1),[],1); n3 = reshape(nod(  end  ,2:end  ,1:end-1),[],1); n7 = reshape(nod(  end  ,2:end  ,2:end  ),[],1); n6 = reshape(nod(  end  ,1:end-1,2:end  ),[],1);   ind = reshape(elm(end,  :,  :),[],1);   tmp = zeros(size(n2,1),6);   tmp( ind,:) = [[n2( ind,:) n6( ind,:) n3( ind,:)] [n7( ind,:) n3( ind,:) n6( ind,:)]];   tmp(~ind,:) = [[n3(~ind,:) n2(~ind,:) n7(~ind,:)] [n6(~ind,:) n7(~ind,:) n2(~ind,:)]];   mesh.bou{2}.nod = [tmp(:,1:3); tmp(:,4:6)];   mesh.bou{2}.ele = [reshape((ele(end,  :,  :)-1)*5+1*eo(end,  :,  :)+2*~eo(end,  :,  :),[],1); reshape((ele(end,  :,  :)-1)*5+4*eo(end,  :,  :)+3*~eo(end,  :,  :),[],1)];   %  2 (x == L)
        n1 = reshape(nod(1:end-1,1      ,1:end-1),[],1); n2 = reshape(nod(2:end  ,1      ,1:end-1),[],1); n6 = reshape(nod(2:end  ,1      ,2:end  ),[],1); n5 = reshape(nod(1:end-1,1      ,2:end  ),[],1);   ind = reshape(elm(  :,  1,  :),[],1);   tmp = zeros(size(n1,1),6);   tmp( ind,:) = [[n2( ind,:) n1( ind,:) n6( ind,:)] [n5( ind,:) n6( ind,:) n1( ind,:)]];   tmp(~ind,:) = [[n1(~ind,:) n5(~ind,:) n2(~ind,:)] [n6(~ind,:) n2(~ind,:) n5(~ind,:)]];   mesh.bou{3}.nod = [tmp(:,1:3); tmp(:,4:6)];   mesh.bou{3}.ele = [reshape((ele(  :,  1,  :)-1)*5+1*eo(  :,  1,  :)+1*~eo(  :,  1,  :),[],1); reshape((ele(  :,  1,  :)-1)*5+3*eo(  :,  1,  :)+3*~eo(  :,  1,  :),[],1)];   %  3 (y == 0)
        n4 = reshape(nod(1:end-1,  end  ,1:end-1),[],1); n3 = reshape(nod(2:end  ,  end  ,1:end-1),[],1); n7 = reshape(nod(2:end  ,  end  ,2:end  ),[],1); n8 = reshape(nod(1:end-1,  end  ,2:end  ),[],1);   ind = reshape(elm(  :,end,  :),[],1);   tmp = zeros(size(n3,1),6);   tmp( ind,:) = [[n4( ind,:) n3( ind,:) n8( ind,:)] [n7( ind,:) n8( ind,:) n3( ind,:)]];   tmp(~ind,:) = [[n3(~ind,:) n7(~ind,:) n4(~ind,:)] [n8(~ind,:) n4(~ind,:) n7(~ind,:)]];   mesh.bou{4}.nod = [tmp(:,1:3); tmp(:,4:6)];   mesh.bou{4}.ele = [reshape((ele(  :,end,  :)-1)*5+2*eo(  :,end,  :)+2*~eo(  :,end,  :),[],1); reshape((ele(  :,end,  :)-1)*5+4*eo(  :,end,  :)+4*~eo(  :,end,  :),[],1)];   %  4 (y == L)
        n1 = reshape(nod(1:end-1,1:end-1,1      ),[],1); n2 = reshape(nod(2:end  ,1:end-1,1      ),[],1); n3 = reshape(nod(2:end  ,2:end  ,1      ),[],1); n4 = reshape(nod(1:end-1,2:end  ,1      ),[],1);   ind = reshape(elm(  :,  :,  1),[],1);   tmp = zeros(size(n1,1),6);   tmp( ind,:) = [[n2( ind,:) n3( ind,:) n1( ind,:)] [n4( ind,:) n1( ind,:) n3( ind,:)]];   tmp(~ind,:) = [[n1(~ind,:) n2(~ind,:) n4(~ind,:)] [n3(~ind,:) n4(~ind,:) n2(~ind,:)]];   mesh.bou{5}.nod = [tmp(:,1:3); tmp(:,4:6)];   mesh.bou{5}.ele = [reshape((ele(  :,  :,  1)-1)*5+1*eo(  :,  :,  1)+1*~eo(  :,  :,  1),[],1); reshape((ele(  :,  :,  1)-1)*5+2*eo(  :,  :,  1)+2*~eo(  :,  :,  1),[],1)];   %  5 (z == 0)
        n5 = reshape(nod(1:end-1,1:end-1,  end  ),[],1); n6 = reshape(nod(2:end  ,1:end-1,  end  ),[],1); n7 = reshape(nod(2:end  ,2:end  ,  end  ),[],1); n8 = reshape(nod(1:end-1,2:end  ,  end  ),[],1);   ind = reshape(elm(  :,  :,end),[],1);   tmp = zeros(size(n5,1),6);   tmp( ind,:) = [[n5( ind,:) n8( ind,:) n6( ind,:)] [n7( ind,:) n6( ind,:) n8( ind,:)]];   tmp(~ind,:) = [[n6(~ind,:) n5(~ind,:) n7(~ind,:)] [n8(~ind,:) n7(~ind,:) n5(~ind,:)]];   mesh.bou{6}.nod = [tmp(:,1:3); tmp(:,4:6)];   mesh.bou{6}.ele = [reshape((ele(  :,  :,end)-1)*5+3*eo(  :,  :,end)+3*~eo(  :,  :,end),[],1); reshape((ele(  :,  :,end)-1)*5+4*eo(  :,  :,end)+4*~eo(  :,  :,end),[],1)];   %  6 (z == L)
    end
  case {'hexa20','hexa27','tetr10'}
    [xx,yy,zz] = meshgrid( 0:L(1)/(2*ne(1)):L(1), 0:L(2)/(2*ne(2)):L(2), 0:L(3)/(2*ne(3)):L(3));   xx = permute(xx,[2 1 3]);   yy = permute(yy,[2 1 3]);   zz = permute(zz,[2 1 3]);
    [ii,jj,kk] = meshgrid( 0:ne(1)              , 0:ne(2)              , 0:ne(3)              );   ii = permute(ii,[2 1 3]);   jj = permute(jj,[2 1 3]);   kk = permute(kk,[2 1 3]);
    mesh.nod.coo = [reshape(xx,[],1) reshape(yy,[],1) reshape(zz,[],1)];
    delcoo = false(size(mesh.nod.coo,1),1);
    nod = reshape(1:size(mesh.nod.coo,1),2*ne(1)+1,2*ne(2)+1,2*ne(3)+1);
    % elements
    %      z=0           z=0.5         z=1
    %  n01 n09 n02   n17 n23 n18   n05 n13 n06
    %  n12 n21 n10   n26 n27 n24   n16 n22 n14   ---> z
    %  n04 n11 n03   n20 n25 n19   n08 n15 n07
    pne = prod(ne);
    ele = reshape(1:pne,ne(1),ne(2),ne(3));
    elm = mod(ii(1:end-1,1:end-1,1:end-1)+jj(1:end-1,1:end-1,1:end-1)+kk(1:end-1,1:end-1,1:end-1),2)==0; % element mirrored
    n01 = reshape(nod(1:2:end-2,1:2:end-2,1:2:end-2),[],1);
    n02 = reshape(nod(3:2:end  ,1:2:end-2,1:2:end-2),[],1);
    n03 = reshape(nod(3:2:end  ,3:2:end  ,1:2:end-2),[],1);
    n04 = reshape(nod(1:2:end-2,3:2:end  ,1:2:end-2),[],1);
    n05 = reshape(nod(1:2:end-2,1:2:end-2,3:2:end  ),[],1);
    n06 = reshape(nod(3:2:end  ,1:2:end-2,3:2:end  ),[],1);
    n07 = reshape(nod(3:2:end  ,3:2:end  ,3:2:end  ),[],1);
    n08 = reshape(nod(1:2:end-2,3:2:end  ,3:2:end  ),[],1);
    n09 = reshape(nod(2:2:end-1,1:2:end-2,1:2:end-2),[],1);
    n10 = reshape(nod(3:2:end  ,2:2:end-1,1:2:end-2),[],1);
    n11 = reshape(nod(2:2:end-1,3:2:end  ,1:2:end-2),[],1);
    n12 = reshape(nod(1:2:end-2,2:2:end-1,1:2:end-2),[],1);
    n13 = reshape(nod(2:2:end-1,1:2:end-2,3:2:end  ),[],1);
    n14 = reshape(nod(3:2:end  ,2:2:end-1,3:2:end  ),[],1);
    n15 = reshape(nod(2:2:end-1,3:2:end  ,3:2:end  ),[],1);
    n16 = reshape(nod(1:2:end-2,2:2:end-1,3:2:end  ),[],1);
    n17 = reshape(nod(1:2:end-2,1:2:end-2,2:2:end-1),[],1);
    n18 = reshape(nod(3:2:end  ,1:2:end-2,2:2:end-1),[],1);
    n19 = reshape(nod(3:2:end  ,3:2:end  ,2:2:end-1),[],1);
    n20 = reshape(nod(1:2:end-2,3:2:end  ,2:2:end-1),[],1);
    n21 = reshape(nod(2:2:end-1,2:2:end-1,1:2:end-2),[],1);
    n22 = reshape(nod(2:2:end-1,2:2:end-1,3:2:end  ),[],1);
    n23 = reshape(nod(2:2:end-1,1:2:end-2,2:2:end-1),[],1);
    n24 = reshape(nod(3:2:end  ,2:2:end-1,2:2:end-1),[],1);
    n25 = reshape(nod(2:2:end-1,3:2:end  ,2:2:end-1),[],1);
    n26 = reshape(nod(1:2:end-2,2:2:end-1,2:2:end-1),[],1);
    n27 = reshape(nod(2:2:end-1,2:2:end-1,2:2:end-1),[],1);
    switch etype
      case 'hexa20'
        mesh.ele.nod = [ n01 n02 n03 n04 n05 n06 n07 n08 n09 n10 n11 n12 n13 n14 n15 n16 n17 n18 n19 n20];
        delcoo([n21;n22;n23;n24;n25;n26;n27]) = true;
      case 'hexa27'
        mesh.ele.nod = [ n01 n02 n03 n04 n05 n06 n07 n08 n09 n10 n11 n12 n13 n14 n15 n16 n17 n18 n19 n20 n21 n22 n23 n24 n25 n26 n27];
      case 'tetr10'
        tmp = zeros(size(n01,1),50);
        ind = reshape(mod(ii(1:end-1,1:end-1,1:end-1)+jj(1:end-1,1:end-1,1:end-1)+kk(1:end-1,1:end-1,1:end-1),2),[],1)==0;
        tmp( ind,:) = [[n01(ind,:) n02(ind,:) n03(ind,:) n06(ind,:) n09(ind,:) n10(ind,:) n21(ind,:) n23(ind,:) n18(ind,:) n24(ind,:)] [n01(ind,:) n03(ind,:) n04(ind,:) n08(ind,:) n21(ind,:) n11(ind,:) n12(ind,:) n26(ind,:) n25(ind,:) n20(ind,:)] [n01(ind,:) n05(ind,:) n06(ind,:) n08(ind,:) n17(ind,:) n13(ind,:) n23(ind,:) n26(ind,:) n16(ind,:) n22(ind,:)] [n03(ind,:) n06(ind,:) n07(ind,:) n08(ind,:) n24(ind,:) n14(ind,:) n19(ind,:) n25(ind,:) n22(ind,:) n15(ind,:)] [n01(ind,:) n06(ind,:) n03(ind,:) n08(ind,:) n23(ind,:) n24(ind,:) n21(ind,:) n26(ind,:) n22(ind,:) n25(ind,:)]];
        tmp(~ind,:) = [[n04(ind,:) n01(ind,:) n02(ind,:) n05(ind,:) n12(ind,:) n09(ind,:) n21(ind,:) n26(ind,:) n17(ind,:) n23(ind,:)] [n02(ind,:) n03(ind,:) n04(ind,:) n07(ind,:) n10(ind,:) n11(ind,:) n21(ind,:) n24(ind,:) n19(ind,:) n25(ind,:)] [n05(ind,:) n06(ind,:) n02(ind,:) n07(ind,:) n13(ind,:) n18(ind,:) n23(ind,:) n22(ind,:) n14(ind,:) n24(ind,:)] [n05(ind,:) n08(ind,:) n07(ind,:) n04(ind,:) n16(ind,:) n15(ind,:) n22(ind,:) n26(ind,:) n20(ind,:) n25(ind,:)] [n02(ind,:) n05(ind,:) n07(ind,:) n04(ind,:) n23(ind,:) n22(ind,:) n24(ind,:) n21(ind,:) n26(ind,:) n25(ind,:)]];
        mesh.ele.nod = zeros(5*size(n01,1),10);
        mesh.ele.nod(1:5:end,:) = tmp(:, 1:10);
        mesh.ele.nod(2:5:end,:) = tmp(:,11:20);
        mesh.ele.nod(3:5:end,:) = tmp(:,21:30);
        mesh.ele.nod(4:5:end,:) = tmp(:,31:40);
        mesh.ele.nod(5:5:end,:) = tmp(:,41:50);
        delcoo(n27) = true;
    end
    % boundary elements
    switch etype
      case 'hexa20'
        n01 = reshape(nod(1        ,1:2:end-2,1:2:end-2),[],1); n04 = reshape(nod(1        ,3:2:end  ,1:2:end-2),[],1); n08 = reshape(nod(1        ,3:2:end  ,3:2:end  ),[],1); n05 = reshape(nod(1        ,1:2:end-2,3:2:end  ),[],1); n12 = reshape(nod(1        ,2:2:end-1,1:2:end-2),[],1); n20 = reshape(nod(1        ,3:2:end  ,2:2:end-1),[],1); n16 = reshape(nod(1        ,2:2:end-1,3:2:end  ),[],1); n17 = reshape(nod(1        ,1:2:end-2,2:2:end-1),[],1);                                                          mesh.bou{1}.nod = [n01 n04 n08 n05 n12 n20 n16 n17    ];                   mesh.bou{1}.ele = reshape(ele(1    ,1:end,1:end),[],1);   %  1 (x == 0)
        n02 = reshape(nod(    end  ,1:2:end-2,1:2:end-2),[],1); n03 = reshape(nod(    end  ,3:2:end  ,1:2:end-2),[],1); n07 = reshape(nod(    end  ,3:2:end  ,3:2:end  ),[],1); n06 = reshape(nod(    end  ,1:2:end-2,3:2:end  ),[],1); n10 = reshape(nod(    end  ,2:2:end-1,1:2:end-2),[],1); n19 = reshape(nod(    end  ,3:2:end  ,2:2:end-1),[],1); n14 = reshape(nod(    end  ,2:2:end-1,3:2:end  ),[],1); n18 = reshape(nod(    end  ,1:2:end-2,2:2:end-1),[],1);                                                          mesh.bou{2}.nod = [n06 n07 n03 n02 n14 n19 n10 n18    ];                   mesh.bou{2}.ele = reshape(ele(  end,1:end,1:end),[],1);   %  2 (x == L)
        n01 = reshape(nod(1:2:end-1,1        ,1:2:end-2),[],1); n02 = reshape(nod(3:2:end  ,1        ,1:2:end-2),[],1); n06 = reshape(nod(3:2:end  ,1        ,3:2:end  ),[],1); n05 = reshape(nod(1:2:end-2,1        ,3:2:end  ),[],1); n09 = reshape(nod(2:2:end-1,1        ,1:2:end-2),[],1); n18 = reshape(nod(3:2:end  ,1        ,2:2:end-1),[],1); n13 = reshape(nod(2:2:end-1,1        ,3:2:end  ),[],1); n17 = reshape(nod(1:2:end-2,1        ,2:2:end-1),[],1);                                                          mesh.bou{3}.nod = [n01 n02 n06 n05 n09 n18 n13 n17    ];                   mesh.bou{3}.ele = reshape(ele(1:end,1    ,1:end),[],1);   %  3 (y == 0)
        n04 = reshape(nod(1:2:end-1,    end  ,1:2:end-2),[],1); n03 = reshape(nod(3:2:end  ,    end  ,1:2:end-2),[],1); n07 = reshape(nod(3:2:end  ,    end  ,3:2:end  ),[],1); n08 = reshape(nod(1:2:end-2,    end  ,3:2:end  ),[],1); n11 = reshape(nod(2:2:end-1,    end  ,1:2:end-2),[],1); n19 = reshape(nod(3:2:end  ,    end  ,2:2:end-1),[],1); n15 = reshape(nod(2:2:end-1,    end  ,3:2:end  ),[],1); n20 = reshape(nod(1:2:end-2,    end  ,2:2:end-1),[],1);                                                          mesh.bou{4}.nod = [n08 n07 n03 n04 n15 n19 n11 n20    ];                   mesh.bou{4}.ele = reshape(ele(1:end,  end,1:end),[],1);   %  4 (y == L)
        n01 = reshape(nod(1:2:end-1,1:2:end-1,1        ),[],1); n02 = reshape(nod(3:2:end  ,1:2:end-2,1        ),[],1); n03 = reshape(nod(3:2:end  ,3:2:end  ,1        ),[],1); n04 = reshape(nod(1:2:end-2,3:2:end  ,1        ),[],1); n09 = reshape(nod(2:2:end-1,1:2:end-2,1        ),[],1); n10 = reshape(nod(3:2:end  ,2:2:end-1,1        ),[],1); n11 = reshape(nod(2:2:end-1,3:2:end  ,1        ),[],1); n12 = reshape(nod(1:2:end-2,2:2:end-1,1        ),[],1);                                                          mesh.bou{5}.nod = [n01 n02 n03 n04 n09 n10 n11 n12    ];                   mesh.bou{5}.ele = reshape(ele(1:end,1:end,1    ),[],1);   %  5 (z == 0)
        n05 = reshape(nod(1:2:end-1,1:2:end-1,    end  ),[],1); n06 = reshape(nod(3:2:end  ,1:2:end-2,    end  ),[],1); n07 = reshape(nod(3:2:end  ,3:2:end  ,    end  ),[],1); n08 = reshape(nod(1:2:end-2,3:2:end  ,    end  ),[],1); n13 = reshape(nod(2:2:end-1,1:2:end-2,    end  ),[],1); n14 = reshape(nod(3:2:end  ,2:2:end-1,    end  ),[],1); n15 = reshape(nod(2:2:end-1,3:2:end  ,    end  ),[],1); n16 = reshape(nod(1:2:end-2,2:2:end-1,    end  ),[],1);                                                          mesh.bou{6}.nod = [n08 n07 n06 n05 n15 n14 n13 n16    ];                   mesh.bou{6}.ele = reshape(ele(1:end,1:end,  end),[],1);   %  6 (z == L)
      case 'hexa27'
        n01 = reshape(nod(1        ,1:2:end-2,1:2:end-2),[],1); n04 = reshape(nod(1        ,3:2:end  ,1:2:end-2),[],1); n08 = reshape(nod(1        ,3:2:end  ,3:2:end  ),[],1); n05 = reshape(nod(1        ,1:2:end-2,3:2:end  ),[],1); n12 = reshape(nod(1        ,2:2:end-1,1:2:end-2),[],1); n20 = reshape(nod(1        ,3:2:end  ,2:2:end-1),[],1); n16 = reshape(nod(1        ,2:2:end-1,3:2:end  ),[],1); n17 = reshape(nod(1        ,1:2:end-2,2:2:end-1),[],1); n26 = reshape(nod(1        ,2:2:end-1,2:2:end-1),[],1);  mesh.bou{1}.nod = [n01 n04 n08 n05 n12 n20 n16 n17 n26];                   mesh.bou{1}.ele = reshape(ele(1    ,1:end,1:end),[],1);   %  1 (x == 0)
        n02 = reshape(nod(    end  ,1:2:end-2,1:2:end-2),[],1); n03 = reshape(nod(    end  ,3:2:end  ,1:2:end-2),[],1); n07 = reshape(nod(    end  ,3:2:end  ,3:2:end  ),[],1); n06 = reshape(nod(    end  ,1:2:end-2,3:2:end  ),[],1); n10 = reshape(nod(    end  ,2:2:end-1,1:2:end-2),[],1); n19 = reshape(nod(    end  ,3:2:end  ,2:2:end-1),[],1); n14 = reshape(nod(    end  ,2:2:end-1,3:2:end  ),[],1); n18 = reshape(nod(    end  ,1:2:end-2,2:2:end-1),[],1); n24 = reshape(nod(    end  ,2:2:end-1,2:2:end-1),[],1);  mesh.bou{2}.nod = [n06 n07 n03 n02 n14 n19 n10 n18 n24];                   mesh.bou{2}.ele = reshape(ele(  end,1:end,1:end),[],1);   %  2 (x == L)
        n01 = reshape(nod(1:2:end-1,1        ,1:2:end-2),[],1); n02 = reshape(nod(3:2:end  ,1        ,1:2:end-2),[],1); n06 = reshape(nod(3:2:end  ,1        ,3:2:end  ),[],1); n05 = reshape(nod(1:2:end-2,1        ,3:2:end  ),[],1); n09 = reshape(nod(2:2:end-1,1        ,1:2:end-2),[],1); n18 = reshape(nod(3:2:end  ,1        ,2:2:end-1),[],1); n13 = reshape(nod(2:2:end-1,1        ,3:2:end  ),[],1); n17 = reshape(nod(1:2:end-2,1        ,2:2:end-1),[],1); n23 = reshape(nod(2:2:end-1,1        ,2:2:end-1),[],1);  mesh.bou{3}.nod = [n01 n02 n06 n05 n09 n18 n13 n17 n23];                   mesh.bou{3}.ele = reshape(ele(1:end,1    ,1:end),[],1);   %  3 (y == 0)
        n04 = reshape(nod(1:2:end-1,    end  ,1:2:end-2),[],1); n03 = reshape(nod(3:2:end  ,    end  ,1:2:end-2),[],1); n07 = reshape(nod(3:2:end  ,    end  ,3:2:end  ),[],1); n08 = reshape(nod(1:2:end-2,    end  ,3:2:end  ),[],1); n11 = reshape(nod(2:2:end-1,    end  ,1:2:end-2),[],1); n19 = reshape(nod(3:2:end  ,    end  ,2:2:end-1),[],1); n15 = reshape(nod(2:2:end-1,    end  ,3:2:end  ),[],1); n20 = reshape(nod(1:2:end-2,    end  ,2:2:end-1),[],1); n25 = reshape(nod(2:2:end-1,    end  ,2:2:end-1),[],1);  mesh.bou{4}.nod = [n08 n07 n03 n04 n15 n19 n11 n20 n25];                   mesh.bou{4}.ele = reshape(ele(1:end,  end,1:end),[],1);   %  4 (y == L)
        n01 = reshape(nod(1:2:end-1,1:2:end-1,1        ),[],1); n02 = reshape(nod(3:2:end  ,1:2:end-2,1        ),[],1); n03 = reshape(nod(3:2:end  ,3:2:end  ,1        ),[],1); n04 = reshape(nod(1:2:end-2,3:2:end  ,1        ),[],1); n09 = reshape(nod(2:2:end-1,1:2:end-2,1        ),[],1); n10 = reshape(nod(3:2:end  ,2:2:end-1,1        ),[],1); n11 = reshape(nod(2:2:end-1,3:2:end  ,1        ),[],1); n12 = reshape(nod(1:2:end-2,2:2:end-1,1        ),[],1); n21 = reshape(nod(2:2:end-1,2:2:end-1,1        ),[],1);  mesh.bou{5}.nod = [n01 n02 n03 n04 n09 n10 n11 n12 n21];                   mesh.bou{5}.ele = reshape(ele(1:end,1:end,1    ),[],1);   %  5 (z == 0)
        n05 = reshape(nod(1:2:end-1,1:2:end-1,    end  ),[],1); n06 = reshape(nod(3:2:end  ,1:2:end-2,    end  ),[],1); n07 = reshape(nod(3:2:end  ,3:2:end  ,    end  ),[],1); n08 = reshape(nod(1:2:end-2,3:2:end  ,    end  ),[],1); n13 = reshape(nod(2:2:end-1,1:2:end-2,    end  ),[],1); n14 = reshape(nod(3:2:end  ,2:2:end-1,    end  ),[],1); n15 = reshape(nod(2:2:end-1,3:2:end  ,    end  ),[],1); n16 = reshape(nod(1:2:end-2,2:2:end-1,    end  ),[],1); n22 = reshape(nod(2:2:end-1,2:2:end-1,    end  ),[],1);  mesh.bou{6}.nod = [n08 n07 n06 n05 n15 n14 n13 n16 n22];                   mesh.bou{6}.ele = reshape(ele(1:end,1:end,  end),[],1);   %  6 (z == L)
      case 'tetr10'
        n01 = reshape(nod(1        ,1:2:end-2,1:2:end-2),[],1); n04 = reshape(nod(1        ,3:2:end  ,1:2:end-2),[],1); n08 = reshape(nod(1        ,3:2:end  ,3:2:end  ),[],1); n05 = reshape(nod(1        ,1:2:end-2,3:2:end  ),[],1); n12 = reshape(nod(1        ,2:2:end-1,1:2:end-2),[],1); n20 = reshape(nod(1        ,3:2:end  ,2:2:end-1),[],1); n16 = reshape(nod(1        ,2:2:end-1,3:2:end  ),[],1); n17 = reshape(nod(1        ,1:2:end-2,2:2:end-1),[],1); n26 = reshape(nod(1        ,2:2:end-1,2:2:end-1),[],1);   ind = reshape(elm(  1,  :,  :),[],1);   mesh.bou{1}.nod = zeros(size(n01,1),12);   mesh.bou{1}.nod(ind,:) = [[n04(ind,:) n08(ind,:) n01(ind,:) n20(ind,:) n26(ind,:) n12(ind,:)] [n05(ind,:) n01(ind,:) n08(ind,:) n17(ind,:) n26(ind,:) n16(ind,:)]];   mesh.bou{1}.nod(~ind,:) = [[n01(~ind,:) n04(~ind,:) n05(~ind,:) n12(~ind,:) n26(~ind,:) n17(~ind,:)] [n08(~ind,:) n05(~ind,:) n04(~ind,:) n16(~ind,:) n26(~ind,:) n20(~ind,:)]];   mesh.bou{1}.nod = [mesh.bou{1}.nod(:,1:6);mesh.bou{1}.nod(:,7:12)];   mesh.bou{1}.ele = [reshape(ele(1    ,1:end,1:end),[],1)+1*pne*ind+0*pne*~ind; reshape(ele(1    ,1:end,1:end),[],1)+2*pne*ind+3*pne*~ind];   %  1 (x == 0)
        n02 = reshape(nod(    end  ,1:2:end-2,1:2:end-2),[],1); n03 = reshape(nod(    end  ,3:2:end  ,1:2:end-2),[],1); n07 = reshape(nod(    end  ,3:2:end  ,3:2:end  ),[],1); n06 = reshape(nod(    end  ,1:2:end-2,3:2:end  ),[],1); n10 = reshape(nod(    end  ,2:2:end-1,1:2:end-2),[],1); n19 = reshape(nod(    end  ,3:2:end  ,2:2:end-1),[],1); n14 = reshape(nod(    end  ,2:2:end-1,3:2:end  ),[],1); n18 = reshape(nod(    end  ,1:2:end-2,2:2:end-1),[],1); n24 = reshape(nod(    end  ,2:2:end-1,2:2:end-1),[],1);   ind = reshape(elm(end,  :,  :),[],1);   mesh.bou{2}.nod = zeros(size(n02,1),12);   mesh.bou{2}.nod(ind,:) = [[n02(ind,:) n06(ind,:) n03(ind,:) n18(ind,:) n24(ind,:) n10(ind,:)] [n07(ind,:) n03(ind,:) n06(ind,:) n19(ind,:) n24(ind,:) n14(ind,:)]];   mesh.bou{2}.nod(~ind,:) = [[n03(~ind,:) n02(~ind,:) n07(~ind,:) n10(~ind,:) n24(~ind,:) n19(~ind,:)] [n06(~ind,:) n07(~ind,:) n02(~ind,:) n14(~ind,:) n24(~ind,:) n18(~ind,:)]];   mesh.bou{2}.nod = [mesh.bou{2}.nod(:,1:6);mesh.bou{2}.nod(:,7:12)];   mesh.bou{2}.ele = [reshape(ele(  end,1:end,1:end),[],1)+0*pne*ind+1*pne*~ind; reshape(ele(  end,1:end,1:end),[],1)+3*pne*ind+2*pne*~ind];   %  2 (x == L)
        n01 = reshape(nod(1:2:end-1,1        ,1:2:end-2),[],1); n02 = reshape(nod(3:2:end  ,1        ,1:2:end-2),[],1); n06 = reshape(nod(3:2:end  ,1        ,3:2:end  ),[],1); n05 = reshape(nod(1:2:end-2,1        ,3:2:end  ),[],1); n09 = reshape(nod(2:2:end-1,1        ,1:2:end-2),[],1); n18 = reshape(nod(3:2:end  ,1        ,2:2:end-1),[],1); n13 = reshape(nod(2:2:end-1,1        ,3:2:end  ),[],1); n17 = reshape(nod(1:2:end-2,1        ,2:2:end-1),[],1); n23 = reshape(nod(2:2:end-1,1        ,2:2:end-1),[],1);   ind = reshape(elm(  :,  1,  :),[],1);   mesh.bou{3}.nod = zeros(size(n01,1),12);   mesh.bou{3}.nod(ind,:) = [[n02(ind,:) n01(ind,:) n06(ind,:) n09(ind,:) n23(ind,:) n18(ind,:)] [n05(ind,:) n06(ind,:) n01(ind,:) n13(ind,:) n23(ind,:) n17(ind,:)]];   mesh.bou{3}.nod(~ind,:) = [[n01(~ind,:) n05(~ind,:) n02(~ind,:) n17(~ind,:) n23(~ind,:) n09(~ind,:)] [n06(~ind,:) n02(~ind,:) n05(~ind,:) n18(~ind,:) n23(~ind,:) n13(~ind,:)]];   mesh.bou{3}.nod = [mesh.bou{3}.nod(:,1:6);mesh.bou{3}.nod(:,7:12)];   mesh.bou{3}.ele = [reshape(ele(1:end,1    ,1:end),[],1)+0*pne*ind+0*pne*~ind; reshape(ele(1:end,1    ,1:end),[],1)+2*pne*ind+2*pne*~ind];   %  3 (y == 0)
        n04 = reshape(nod(1:2:end-1,    end  ,1:2:end-2),[],1); n03 = reshape(nod(3:2:end  ,    end  ,1:2:end-2),[],1); n07 = reshape(nod(3:2:end  ,    end  ,3:2:end  ),[],1); n08 = reshape(nod(1:2:end-2,    end  ,3:2:end  ),[],1); n11 = reshape(nod(2:2:end-1,    end  ,1:2:end-2),[],1); n19 = reshape(nod(3:2:end  ,    end  ,2:2:end-1),[],1); n15 = reshape(nod(2:2:end-1,    end  ,3:2:end  ),[],1); n20 = reshape(nod(1:2:end-2,    end  ,2:2:end-1),[],1); n25 = reshape(nod(2:2:end-1,    end  ,2:2:end-1),[],1);   ind = reshape(elm(  :,end,  :),[],1);   mesh.bou{4}.nod = zeros(size(n03,1),12);   mesh.bou{4}.nod(ind,:) = [[n04(ind,:) n03(ind,:) n08(ind,:) n11(ind,:) n25(ind,:) n20(ind,:)] [n07(ind,:) n08(ind,:) n03(ind,:) n15(ind,:) n25(ind,:) n19(ind,:)]];   mesh.bou{4}.nod(~ind,:) = [[n03(~ind,:) n07(~ind,:) n04(~ind,:) n19(~ind,:) n25(~ind,:) n11(~ind,:)] [n08(~ind,:) n04(~ind,:) n07(~ind,:) n20(~ind,:) n25(~ind,:) n15(~ind,:)]];   mesh.bou{4}.nod = [mesh.bou{4}.nod(:,1:6);mesh.bou{4}.nod(:,7:12)];   mesh.bou{4}.ele = [reshape(ele(1:end,  end,1:end),[],1)+1*pne*ind+1*pne*~ind; reshape(ele(1:end,  end,1:end),[],1)+3*pne*ind+3*pne*~ind];   %  4 (y == L)
        n01 = reshape(nod(1:2:end-1,1:2:end-1,1        ),[],1); n02 = reshape(nod(3:2:end  ,1:2:end-2,1        ),[],1); n03 = reshape(nod(3:2:end  ,3:2:end  ,1        ),[],1); n04 = reshape(nod(1:2:end-2,3:2:end  ,1        ),[],1); n09 = reshape(nod(2:2:end-1,1:2:end-2,1        ),[],1); n10 = reshape(nod(3:2:end  ,2:2:end-1,1        ),[],1); n11 = reshape(nod(2:2:end-1,3:2:end  ,1        ),[],1); n12 = reshape(nod(1:2:end-2,2:2:end-1,1        ),[],1); n21 = reshape(nod(2:2:end-1,2:2:end-1,1        ),[],1);   ind = reshape(elm(  :,  :,  1),[],1);   mesh.bou{5}.nod = zeros(size(n01,1),12);   mesh.bou{5}.nod(ind,:) = [[n02(ind,:) n03(ind,:) n01(ind,:) n10(ind,:) n21(ind,:) n09(ind,:)] [n04(ind,:) n01(ind,:) n03(ind,:) n12(ind,:) n21(ind,:) n11(ind,:)]];   mesh.bou{5}.nod(~ind,:) = [[n01(~ind,:) n02(~ind,:) n04(~ind,:) n09(~ind,:) n21(~ind,:) n12(~ind,:)] [n03(~ind,:) n04(~ind,:) n02(~ind,:) n11(~ind,:) n21(~ind,:) n10(~ind,:)]];   mesh.bou{5}.nod = [mesh.bou{5}.nod(:,1:6);mesh.bou{5}.nod(:,7:12)];   mesh.bou{5}.ele = [reshape(ele(1:end,1:end,1    ),[],1)+0*pne*ind+0*pne*~ind; reshape(ele(1:end,1:end,1    ),[],1)+1*pne*ind+1*pne*~ind];   %  5 (z == 0)
        n05 = reshape(nod(1:2:end-1,1:2:end-1,    end  ),[],1); n06 = reshape(nod(3:2:end  ,1:2:end-2,    end  ),[],1); n07 = reshape(nod(3:2:end  ,3:2:end  ,    end  ),[],1); n08 = reshape(nod(1:2:end-2,3:2:end  ,    end  ),[],1); n13 = reshape(nod(2:2:end-1,1:2:end-2,    end  ),[],1); n14 = reshape(nod(3:2:end  ,2:2:end-1,    end  ),[],1); n15 = reshape(nod(2:2:end-1,3:2:end  ,    end  ),[],1); n16 = reshape(nod(1:2:end-2,2:2:end-1,    end  ),[],1); n22 = reshape(nod(2:2:end-1,2:2:end-1,    end  ),[],1);   ind = reshape(elm(  :,  :,end),[],1);   mesh.bou{6}.nod = zeros(size(n05,1),12);   mesh.bou{6}.nod(ind,:) = [[n05(ind,:) n08(ind,:) n06(ind,:) n16(ind,:) n22(ind,:) n13(ind,:)] [n07(ind,:) n06(ind,:) n08(ind,:) n14(ind,:) n22(ind,:) n15(ind,:)]];   mesh.bou{6}.nod(~ind,:) = [[n06(~ind,:) n05(~ind,:) n07(~ind,:) n13(~ind,:) n22(~ind,:) n14(~ind,:)] [n08(~ind,:) n07(~ind,:) n05(~ind,:) n15(~ind,:) n22(~ind,:) n16(~ind,:)]];   mesh.bou{6}.nod = [mesh.bou{6}.nod(:,1:6);mesh.bou{6}.nod(:,7:12)];   mesh.bou{6}.ele = [reshape(ele(1:end,1:end,  end),[],1)+2*pne*ind+2*pne*~ind; reshape(ele(1:end,1:end,  end),[],1)+3*pne*ind+3*pne*~ind];   %  6 (z == L)
    end
    % delete redundant coordinates
    switch etype
      case {'hexa20','tetr10'}
        map = zeros(size(delcoo,1),1);   map(~delcoo,1) = 1:sum(~delcoo);
        mesh.nod.coo(delcoo,:) = [];                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
        mesh.ele.nod    = map(   mesh.ele.nod);   if size(mesh.ele.nod   ,2)==1, mesh.ele.nod    = reshape(mesh.ele.nod   ,1,[]); end
        mesh.bou{1}.nod = map(mesh.bou{1}.nod);   if size(mesh.bou{1}.nod,2)==1, mesh.bou{1}.nod = reshape(mesh.bou{1}.nod,1,[]); end
        mesh.bou{2}.nod = map(mesh.bou{2}.nod);   if size(mesh.bou{2}.nod,2)==1, mesh.bou{2}.nod = reshape(mesh.bou{2}.nod,1,[]); end
        mesh.bou{3}.nod = map(mesh.bou{3}.nod);   if size(mesh.bou{3}.nod,2)==1, mesh.bou{3}.nod = reshape(mesh.bou{3}.nod,1,[]); end
        mesh.bou{4}.nod = map(mesh.bou{4}.nod);   if size(mesh.bou{4}.nod,2)==1, mesh.bou{4}.nod = reshape(mesh.bou{4}.nod,1,[]); end
        mesh.bou{5}.nod = map(mesh.bou{5}.nod);   if size(mesh.bou{5}.nod,2)==1, mesh.bou{5}.nod = reshape(mesh.bou{5}.nod,1,[]); end
        mesh.bou{6}.nod = map(mesh.bou{6}.nod);   if size(mesh.bou{6}.nod,2)==1, mesh.bou{6}.nod = reshape(mesh.bou{6}.nod,1,[]); end
    end
end
mesh.info.noddof      = d;
mesh.info.dim         = d;
mesh.info.nodcount    = size(mesh.nod.coo,1);
mesh.info.elecount    = size(mesh.ele.nod,1);
mesh.info.boucount    = zeros(mesh.info.boucounts,1);  for i = 1:mesh.info.boucount,   mesh.info.boucounts(i) = size(mesh.bou{i}.nod,1);  end
mesh.info.extent(1,:) = [0 L(1)]; mesh.info.extent(2,:) = [0 L(2)];    if d>2,   mesh.info.extent(3,:) = [0 L(3)]; end

