addpath(genpath('C:/users/petnov/Dropbox/shared - copy/'))
enableVTK;
regions =load('D:/ARVC meshing automatic/patients/patient06/regions_full2.txt');

mesh = read_CHASTE('D:/ARVC meshing automatic/patients/patient06/chaste06/HEART');
mesh.xyzREGIONS=regions(:,1);
mesh.celltype =10;

write_VTK(mesh, 'D:/ARVC meshing automatic/patients/patient06/fibrosis.vtk','binary');