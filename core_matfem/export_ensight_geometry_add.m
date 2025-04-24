function export_ensight_geometry_add(name,mesh,ico)

% % ALTER .case
% fid = fopen(['results/',name,'.case'],'w');
% fprintf(fid,'#Date\n');
% fprintf(fid,'# EnSight Gold Model\n');
% fprintf(fid,'\n');
% fprintf(fid,'FORMAT\n');
% fprintf(fid,'type:                   ensight gold\n');
% fprintf(fid,'\n');
% fprintf(fid,'GEOMETRY\n');
% fprintf(fid,'model:                  %s.geo\n',name);
% fprintf(fid,'\n');
% fprintf(fid,'VARIABLES\n');
% if nargin == 5
%   for i = 1: size(exports,1)
%     switch exports{i,2}
%       case 'nod.vec'   , fprintf(fid,'vector      per node:       1 1 %s %s_%s.Nvec***\n',exports{i,1},name,exports{i,1});
%       case 'ele.tensym', fprintf(fid,'tensor symm per element:    1 1 %s %s_%s.Eten***\n',exports{i,1},name,exports{i,1});
%       case 'ele.sca'   , fprintf(fid,'scalar      per element:    1 1 %s %s_%s.Esca***\n',exports{i,1},name,exports{i,1});
%       otherwise        , fprintf('mesh_export_geometry: unknown export type:%s\n',exports{i,2});
%     end
%   end
% else
%   fprintf(fid,'vector      per node:    1 1 displacement %s_Displacement.Nvec***\n',name);
%   fprintf(fid,'tensor symm per element: 1 1 stress       %s_Stress.Eten***\n'      ,name);
% end
% fprintf(fid,'TIME\n');
% fprintf(fid,'time set:              1\n');
% fprintf(fid,'number of steps:       %d\n',timestepping.steps+1);
% fprintf(fid,'filename start number: 0\n');
% fprintf(fid,'filename increment:    1\n');
% fprintf(fid,'time values:');
% for k = 1:timestepping.steps+1
%   fprintf(fid,' %d',timestepping.times(k));
%   if (k > 0) && (mod(k,10) == 0)
%     fprintf(fid,'\n');
%   end
% end
% fprintf(fid,'\n');
% fclose(fid);

% ALTER .geo
fid = fopen(['results/',name,'.geo'],'a');
fprintf(fid,'part\n');
fprintf(fid,'         %d\n',ico+1);
fprintf(fid,'body ->> 0  0  0  %d\n',ico+1);
fprintf(fid,'coordinates\n');
fprintf(fid,'%10d\n',mesh.info.nodcount);
for n = 1:(mesh.info.nodcount),   fprintf(fid,'%10d\n',n);   end
for n = 1:mesh.info.nodcount,   fprintf(fid,' %1.5e\n',mesh.nod.coo(n,1)); end
for n = 1:mesh.info.nodcount,   fprintf(fid,' %1.5e\n',mesh.nod.coo(n,2)); end
for n = 1:mesh.info.nodcount,   fprintf(fid,' %1.5e\n',mesh.nod.coo(n,3)); end

switch mesh.ele.type
  case 'quad4',    fprintf(fid,'quad4\n%10d\n' ,mesh.info.elecount);    for e = 1:mesh.info.elecount,   fprintf(fid,'%10d\n',e);   end;    for e = 1:mesh.info.elecount,   tmp = mesh.ele.nod(e,:); fprintf(fid,'%10d%10d%10d%10d\n'                                                                ,tmp(1),tmp(2),tmp(3),tmp(4)                                                                                                                           );   end
  case 'quad8',    fprintf(fid,'quad8\n%10d\n' ,mesh.info.elecount);    for e = 1:mesh.info.elecount,   fprintf(fid,'%10d\n',e);   end;    for e = 1:mesh.info.elecount,   tmp = mesh.ele.nod(e,:); fprintf(fid,'%10d%10d%10d%10d%10d%10d%10d%10d\n'                                                ,tmp(1),tmp(2),tmp(3),tmp(4),tmp(5),tmp(6),tmp(7),tmp(8)                                                                                               );   end
  case 'tetr4',    fprintf(fid,'tetra4\n%10d\n',mesh.info.elecount);    for e = 1:mesh.info.elecount,   fprintf(fid,'%10d\n',e);   end;    for e = 1:mesh.info.elecount,   tmp = mesh.ele.nod(e,:); fprintf(fid,'%10d%10d%10d%10d\n'                                                                ,tmp(1),tmp(2),tmp(3),tmp(4)                                                                                                                           );   end
  case 'hexa8',    fprintf(fid,'hexa8\n%10d\n' ,mesh.info.elecount);    for e = 1:mesh.info.elecount,   fprintf(fid,'%10d\n',e);   end;    for e = 1:mesh.info.elecount,   tmp = mesh.ele.nod(e,:); fprintf(fid,'%10d%10d%10d%10d%10d%10d%10d%10d\n'                                                ,tmp(1),tmp(2),tmp(3),tmp(4),tmp(5),tmp(6),tmp(7),tmp(8)                                                                                               );   end
  case 'hexa20',   fprintf(fid,'hexa20\n%10d\n',mesh.info.elecount);    for e = 1:mesh.info.elecount,   fprintf(fid,'%10d\n',e);   end;    for e = 1:mesh.info.elecount,   tmp = mesh.ele.nod(e,:); fprintf(fid,'%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d\n',tmp(1),tmp(2),tmp(3),tmp(4),tmp(5),tmp(6),tmp(7),tmp(8),tmp(9),tmp(10),tmp(11),tmp(12),tmp(13),tmp(14),tmp(15),tmp(16),tmp(17),tmp(18),tmp(19),tmp(20));   end
  otherwise,    error(['Element type (',mesh.ele.type,') not programmed yet']);
end
fclose(fid);

 