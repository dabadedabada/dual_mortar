function delt_Jcell = linear_Jcell(cell, delt_cell)

term1 = cross(cell(2,:)-cell(1,:),cell(3,:)-cell(1,:), 2)/norm(cross(cell(2,:)-cell(1,:),cell(3,:)-cell(1,:), 2));
term2 = cross(delt_cell(2,:)-delt_cell(1,:),cell(3,:)-cell(1,:), 2);
term3 = cross(cell(2,:)-cell(1,:),delt_cell(3,:)-delt_cell(1,:), 2);
delt_Jcell = dot(term1, term2+term3);