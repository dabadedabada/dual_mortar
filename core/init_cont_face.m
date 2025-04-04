function cont_face = init_cont_face(mesh, face)
% Init contact face

% Init finite element
cont_face.type = mesh.bou{face}.type;
cont_face.fe = fe_init(mesh.bou{face}.type);

% Rows of maps are local (face) indexes, values are global indexes
cont_face.map.ele = mesh.bou{face}.ele;
cont_face.map.nod = unique(mesh.bou{face}.nod(:));

% Coordinates but for local node indexes
cont_face.coo = mesh.nod.coo(cont_face.map.nod,:);
% On each row (local element index) local nodes indexes in anti-clockwise
cont_face.nod = init_local_nods(cont_face.map.nod, mesh.bou{face}.nod);

% Sizes
cont_face.info.nN = length(cont_face.map.nod);  % Number of nodes on face
cont_face.info.nele = size(cont_face.map.ele,1); % Number of ele on face
cont_face.info.nN_ele = size(cont_face.nod,2);  % % Number of nodes of ele