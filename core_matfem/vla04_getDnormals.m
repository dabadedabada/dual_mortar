function Dn = vla04_getDnormals(mesh,u)
nn = size(mesh.coo,1);
Dn = zeros(nn,3,nn,3);

switch mesh.type
  case 'quad4', centrxi = [0   0  ]';  nodxi = [-1 -1; 1 -1; 1 1; -1 1];
  case 'tria3', centrxi = [0.5 0.5]';  nodxi = [0 0; 1 0; 0 1];
end
FE    = fe_init(mesh.type);
dNdxi_cen = FE.dNdxi(centrxi);
dNdxi_nod = zeros(mesh.info.nN,2,mesh.info.nN);
switch mesh.type
  case 'quad4', dNdxi_nod(:,:,1) = FE.dNdxi(nodxi(1,:)');  dNdxi_nod(:,:,2) = FE.dNdxi(nodxi(2,:)');  dNdxi_nod(:,:,3) = FE.dNdxi(nodxi(3,:)');  dNdxi_nod(:,:,4) = FE.dNdxi(nodxi(4,:)');
  case 'tria3', dNdxi_nod(:,:,1) = FE.dNdxi(nodxi(1,:)');  dNdxi_nod(:,:,2) = FE.dNdxi(nodxi(2,:)');  dNdxi_nod(:,:,3) = FE.dNdxi(nodxi(3,:)');
end

coo   = mesh.coo+reshape(u,3,[])';
cntrornod = 'nod';
Dn_hat = zeros(mesh.info.nele,mesh.info.nN_ele,3,nn,3);
for i = 1:mesh.info.nele
  for j = 1:mesh.info.nN
    switch mesh.type
      case 'quad4'
        % (A3) 
        dNd1x = [dNdxi_nod(:,1,1)'; dNdxi_nod(:,1,2)'; dNdxi_nod(:,1,3)'; dNdxi_nod(:,1,4)']*coo; % sum(k,1,4,dNd1*x)
        dNd2x = [dNdxi_nod(:,2,1)'; dNdxi_nod(:,2,2)'; dNdxi_nod(:,2,3)'; dNdxi_nod(:,2,4)']*coo; % sum(k,1,4,dNd2*x)
        dNd1  = (dNdxi_nod(:,1,1) + dNdxi_nod(:,1,2) + dNdxi_nod(:,1,3) + dNdxi_nod(:,1,4))';     % sum(k,1,4,dNd1)
        dNd2  = (dNdxi_nod(:,2,1) + dNdxi_nod(:,2,2) + dNdxi_nod(:,2,3) + dNdxi_nod(:,2,4))';     % sum(k,1,4,dNd2)
        
      case 'tria3'

    end
    %Dn_hat(i,:,j,m) = 1;
  end
end


normals_ele = zeros(nn,3);
hist = zeros(nn,1);
for e = 1:mesh.info.nele
  nods = mesh.nod(e,:);
  coos = coo(nods,:);
  switch cntrornod
    case 'nod'
      switch mesh.type
        case 'quad4'
          t01 = [dNdxi_nod(:,1,1)'*coos; dNdxi_nod(:,1,2)'*coos; dNdxi_nod(:,1,3)'*coos; dNdxi_nod(:,1,4)'*coos];
          t02 = [dNdxi_nod(:,2,1)'*coos; dNdxi_nod(:,2,2)'*coos; dNdxi_nod(:,2,3)'*coos; dNdxi_nod(:,2,4)'*coos];
          n0  = -[cross(t01(1,:),t02(1,:)); cross(t01(2,:),t02(2,:)); cross(t01(3,:),t02(3,:)); cross(t01(4,:),t02(4,:))];
          n0  = n0./sqrt(sum(n0.^2,2));
          normals_ele(nods(1),:) = normals_ele(nods(1),:) + n0(1,:);
          normals_ele(nods(2),:) = normals_ele(nods(2),:) + n0(2,:);
          normals_ele(nods(3),:) = normals_ele(nods(3),:) + n0(3,:);
          normals_ele(nods(4),:) = normals_ele(nods(4),:) + n0(4,:);
        case 'tria3'
          t01 = [dNdxi_nod(:,1,1)'*coos; dNdxi_nod(:,1,2)'*coos; dNdxi_nod(:,1,3)'*coos];
          t02 = [dNdxi_nod(:,2,1)'*coos; dNdxi_nod(:,2,2)'*coos; dNdxi_nod(:,2,3)'*coos];
          n0  = -[cross(t01(1,:),t02(1,:)); cross(t01(2,:),t02(2,:)); cross(t01(3,:),t02(3,:))];
          n0  = n0./sqrt(sum(n0.^2,2));
          normals_ele(nods(1),:) = normals_ele(nods(1),:) + n0(1,:);
          normals_ele(nods(2),:) = normals_ele(nods(2),:) + n0(2,:);
          normals_ele(nods(3),:) = normals_ele(nods(3),:) + n0(3,:);
      end
    otherwise
      t01 = dNdxi_cen(:,1)'*coos;
      t02 = dNdxi_cen(:,2)'*coos;
      n0  = -cross(t01,t02);   n0 = n0/norm(n0);
      normals_ele(nods,:) = normals_ele(nods,:) + n0;
  end
  hist(nods) = hist(nods) + 1;
end
normals_ele = normals_ele./hist;
normals_ele = normals_ele./sqrt(sum(normals_ele.^2,2));
n = normals_ele;