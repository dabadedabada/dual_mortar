function export_ensight_result(name,result,description,time_string)

if isempty(time_string),  time_string = '';   end
switch description.eletype
  case {'quad4','quad8','tria3','tria6'}
    dim = 2;
  otherwise
    dim = 3;
end
switch description.eletype
  case {'tetr4','tetr10'}
    description_ensight = ['tetra',description.eletype(5:end)];
  otherwise
    description_ensight = description.eletype;
end
switch description.type
  case 'nod.sca'
    fid = fopen(sprintf('results/%s_%s.Nsca%s',name,description.name,time_string),'w');
    fprintf(fid,'This is the description of the EnSight Gold geometry, MODEL:%s\npart\n         1\ncoordinates\n',name);
    for i = 1:length(result),           fprintf(fid,' %+12.5e\n',full(result(i)));   end
  case 'ele.sca'
    fid = fopen(sprintf('results/%s_%s.Esca%s',name,description.name,time_string),'w');
    fprintf(fid,'This is the description of the EnSight Gold geometry, MODEL:%s\npart\n         1\n%s\n'         ,name,description_ensight);
    for i = 1:length(result),           fprintf(fid,' %+12.5e\n',result(i));   end
  case 'nod.vec'
    fid = fopen(sprintf('results/%s_%s.Nvec%s',name,description.name,time_string),'w');
    fprintf(fid,'This is the description of the EnSight Gold geometry, MODEL:%s\npart\n         1\ncoordinates\n',name);
    switch dim
      case 2
        for i = 1:2:length(result)-1,   fprintf(fid,' %+12.5e\n',result(i));   end
        for i = 2:2:length(result)  ,   fprintf(fid,' %+12.5e\n',result(i));   end
        for i = 2:2:length(result)  ,   fprintf(fid,' %+12.5e\n',0        );   end
      case 3
        for i = 1:3:length(result)-2,   fprintf(fid,' %+12.5e\n',result(i));   end
        for i = 2:3:length(result)-1,   fprintf(fid,' %+12.5e\n',result(i));   end
        for i = 3:3:length(result)  ,   fprintf(fid,' %+12.5e\n',result(i));   end
    end
  case 'ele.vec'
    fid = fopen(sprintf('results/%s_%s.Evec%s',name,description.name,time_string),'w');
    fprintf(fid,'This is the description of the EnSight Gold geometry, MODEL:%s\npart\n         1\n%s\n'         ,name,description_ensighte);
    switch dim
      case 2
        for i = 1:2:length(result)-1,   fprintf(fid,' %+12.5e\n',result(i));   end
        for i = 2:2:length(result)  ,   fprintf(fid,' %+12.5e\n',result(i));   end
        for i = 2:2:length(result)  ,   fprintf(fid,' %+12.5e\n',0        );   end
      case 3
        for i = 1:3:length(result)-2,   fprintf(fid,' %+12.5e\n',result(i));   end
        for i = 2:3:length(result)-1,   fprintf(fid,' %+12.5e\n',result(i));   end
        for i = 3:3:length(result)  ,   fprintf(fid,' %+12.5e\n',result(i));   end
    end
  case 'nod.tensym' % TODO
    fid = fopen(sprintf('results/%s_%s.Nten%s',name,description.name,time_string),'w');
    fprintf(fid,'This is the description of the EnSight Gold geometry, MODEL:%s\npart\n         1\ncoordinates\n',name);
    for i = 1:6:length(result)-5,   fprintf(fid,' %+12.5e\n',result(i));   end
    for i = 2:6:length(result)-4,   fprintf(fid,' %+12.5e\n',result(i));   end
    for i = 3:6:length(result)-3,   fprintf(fid,' %+12.5e\n',result(i));   end
    for i = 4:6:length(result)-2,   fprintf(fid,' %+12.5e\n',result(i));   end
    for i = 5:6:length(result)-1,   fprintf(fid,' %+12.5e\n',result(i));   end
    for i = 6:6:length(result)  ,   fprintf(fid,' %+12.5e\n',result(i));   end
  case 'ele.tensym'
    fid = fopen(sprintf('results/%s_%s.Eten%s',name,description.name,time_string),'w');
    fprintf(fid,'This is the description of the EnSight Gold geometry, MODEL:%s\npart\n         1\n\n%s\n'         ,name,description_ensight);
    switch dim
      case 2
        for i = 1:3:length(result)-2,   fprintf(fid,' %+12.5e\n',result(i));   end
        for i = 2:3:length(result)-1,   fprintf(fid,' %+12.5e\n',result(i));   end
        for i = 3:3:length(result)  ,   fprintf(fid,' %+12.5e\n',result(i));   end
        for i = 3:3:length(result)  ,   fprintf(fid,' %+12.5e\n',0        );   end
        for i = 3:3:length(result)  ,   fprintf(fid,' %+12.5e\n',0        );   end
        for i = 3:3:length(result)  ,   fprintf(fid,' %+12.5e\n',0        );   end
      case 3
        for i = 1:6:length(result)-5,   fprintf(fid,' %+12.5e\n',result(i));   end
        for i = 2:6:length(result)-4,   fprintf(fid,' %+12.5e\n',result(i));   end
        for i = 3:6:length(result)-3,   fprintf(fid,' %+12.5e\n',result(i));   end
        for i = 4:6:length(result)-2,   fprintf(fid,' %+12.5e\n',result(i));   end
        for i = 5:6:length(result)-1,   fprintf(fid,' %+12.5e\n',result(i));   end
        for i = 6:6:length(result)  ,   fprintf(fid,' %+12.5e\n',result(i));   end
    end
  otherwise
    fprintf('ERROR - mesh_export_result: undefined element (%s)\n',description_ensight);
end
fclose(fid);

