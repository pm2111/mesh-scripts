%%inspect mesh quality 

M = read_CHASTE('D:\ARVC meshing automatic\patients\patient18\chaste_gmsh\HEART');
%M.celltype =5;
[edge_length,radius_ratio] = meshQuality(M,'edgelengths','radiusratio');
figure()
hist(edge_length)

M = read_CHASTE('D:\meshes\DTI24\HEART');
[edge_length24,radius_ratio24] = meshQuality(M,'edgelengths','radiusratio');

figure()
hist(edge_length24)

edges_long = find(edge_length >0.07);

edges_long24 = find(edge_length24 >0.07);