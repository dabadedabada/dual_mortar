function loc_nod = init_local_nods(map_nod, glob_nod)

n = length(map_nod);
loc_nod = zeros(size(glob_nod));

for i=1:n
  loc_nod(glob_nod==map_nod(i))=i;
end