%eikonal plotter
addpath(genpath('/home/scratch/petnov/Dropbox/shared - copy/IO/'));
addpath(genpath('/home/scratch/petnov/Dropbox/shared - copy/'));
addpath(genpath('/home/scratch/petnov/Dropbox/mesh scripts/'));


ATs=csvread('/data/Dphil/Eikonal/metaData/ActivationTimesDTI003_coarse_EndoSpeed120_EpiSpeed1x.csv');
nodes=csvread('/data/Dphil/Eikonal/metaData/DTI003_coarse_xyz.csv');
tri=csvread('/data/Dphil/Eikonal/metaData/DTI003_coarse_tri.csv');
mesh=struct;

mesh.tri=tri;
mesh.xyz=nodes;

mesh.xyzATs=ATs;

figure()
hplotMESH(mesh)
figure()
h = patch('vertices', mesh.xyz, 'faces', mesh.tri, 'facecolor', 'interp', 'cdata', ATs, 'edgecolor', 'none') ;
colorbar()