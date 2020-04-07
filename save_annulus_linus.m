addpath(genpath('/home/scratch/petnov/Dropbox'));

mesh = read_CHASTE('/data/sim3d/webinar_example/data/2Dmeshcm');

fibrosis = mesh.xyz(:,1) >0;
indices = find(fibrosis ==1);
indices_chaste = indices-1;

figure()
hold on
hplotMESH(mesh)
scatter(mesh.xyz(indices,1),mesh.xyz(indices,2))



fid = fopen('D:/meshes/fibrosis_annulus_indices.txt','w');
fprintf(fid,'%d\n',indices_chaste);
fclose(fid)

fid = fopen('D:/meshes/fibrosis_annulus_x_coords.txt','w');
fprintf(fid,'%5.7f \n',[mesh.xyz(indices,1)]);
fclose(fid)
fid = fopen('D:/meshes/fibrosis_annulus_y_coords.txt','w');
fprintf(fid,'%5.7f \n',[mesh.xyz(indices,2)]);
fclose(fid)

save indices.txt [mesh.xyz(indices,1), mesh.xyz(indices,2)]  -ascii

write_VTK(mesh_working,'/data/sim3d/webinar_example/annulus.vtk');

plotMESH(mesh1)
