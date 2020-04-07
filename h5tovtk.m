%file to read h5 files from chaste sims and put them into vtk format

%LINUX filepaths
addpath(genpath('/users/petnov/Dropbox/'));


fileinfo = h5info('/data/sim3d/louie/results.h5');

activation_times = transpose(h5read(fileinfo.Filename,strcat('/',fileinfo.Datasets(3).Name)));

permutation = load('/data/sim3d/louie/permutation.txt');
activation_times = activation_times(permutation+1); %activation times are in the CHASTE node ordering

mesh = read_VTK('/data/sim3d/louie/HEART_with_UPS.vtk'); %the mesh is in VTK node ordering

mesh.xyzUPS2 = activation_times(1:size(mesh.xyz,1));

write_VTK_UNSTRUCTURED_GRID(mesh,'/data/sim3d/louie/vtk_test.vtk','binary');


%WINDOWS filepaths

addpath(genpath('D:\shared\'));

res_folder = 'D:\sim3d\patient02\';

mesh_path = 'D:\ARVC meshing automatic\patients\patient02\chaste_ernesto\HEART.vtk';

fileinfo = h5info(strcat(res_folder,'results.h5'));

activation_times = transpose(h5read(fileinfo.Filename,strcat('/',fileinfo.Datasets(3).Name)));

permutation = load(strcat(res_folder,'permutation.txt'));
activation_times = activation_times(permutation+1); %activation times are in the CHASTE node ordering

mesh = read_VTK(mesh_path); %the mesh is in VTK node ordering

mesh.xyzUPS = activation_times(1:size(mesh.xyz,1));

write_VTK_UNSTRUCTURED_GRID(mesh,strcat(res_folder,'HEARTwithUPs.vtk','binary');