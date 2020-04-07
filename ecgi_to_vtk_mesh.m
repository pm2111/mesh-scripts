addpath(genpath('C:/users/petnov/Dropbox/'));
ecgi = load('D:/ChrisReconstructions/ARVC9/ARVC9_Baseline.mat');

mesh = struct();
mesh.xyz = ecgi.arvc9_baseline.nodes;
mesh.tri = ecgi.arvc9_baseline.mesh;
mesh.xyzActs=ecgi.arvc9_baseline.nodeparam.timestamps.activationTimes-ecgi.arvc9_baseline.nodeparam.timestamps.activationStart;
mesh.Celltype = 10;

write_VTK(mesh,'D:/ARVC meshing automatic/patients/patient09/ecgi.vtk','binary');
