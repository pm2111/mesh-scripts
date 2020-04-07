MPP_version = 'arvc:1.0';

try
isMPPdir = @(d) isdir(d) && isdir( fullfile( d , 'MeshPersonalizationPipeline' ) );
MPP_source_FOLDER = '';
if ~isMPPdir( MPP_source_FOLDER ), MPP_source_FOLDER = fileparts(fileparts(which('MeshPersonalizationPipeline.m')));  end
if ~isMPPdir( MPP_source_FOLDER ), MPP_source_FOLDER = fileparts(fileparts(mfilename('fullpath')));                   end
if ~isMPPdir( MPP_source_FOLDER ), MPP_source_FOLDER = fileparts(mfilename('fullpath'));                              end
if ~isMPPdir( MPP_source_FOLDER ), MPP_source_FOLDER = pwd;                                                           end
if ~isMPPdir( MPP_source_FOLDER ), MPP_source_FOLDER = fileparts(pwd);                                                end
if ~isMPPdir( MPP_source_FOLDER ), MPP_source_FOLDER = 'E:\Dropbox\shared\';                                          end
if ~isMPPdir( MPP_source_FOLDER ), MPP_source_FOLDER = 'C:\Dropbox\Vigente\shared\';                                  end
if ~isMPPdir( MPP_source_FOLDER ), MPP_source_FOLDER = '/users/petnov/Dropbox/shared';                      end
clearvars('isMPPdir');
MPP_source_FOLDER = fullfile( MPP_source_FOLDER , 'MeshPersonalizationPipeline' );
end

fprintf('MPP_source_FOLDER is :   "%s"\n' , MPP_source_FOLDER );

TORSO_MODEL_DIR = fullfile(MPP_source_FOLDER,'TORSO');
fprintf('TORSO_FOLDER is :        "%s"\n' , TORSO_MODEL_DIR );

try
try,OrbitPanZoom('never');end
ww__=warning;warning('off','all');rmpath(genpath('E:\Dropbox\'));warning(ww__);clearvars('ww__');
ww__=warning;warning('off','all');rmpath(genpath('C:\Dropbox\'));warning(ww__);clearvars('ww__');
ww__=warning;warning('off','all');rmpath(genpath(MPP_source_FOLDER));warning(ww__);clearvars('ww__');

addpath(           MPP_source_FOLDER );
addpath( fullfile( MPP_source_FOLDER , 'MPP_tools/' ) );

addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'thirdParty/' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'thirdParty/export_fig' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'IO/' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'LieAlgebra/' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'MESH/' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'MESHES/' ) ); 
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'MESHES/jigsaw/' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'MESHES/tetgen/' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'MESHES/carve/' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'Tools/' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'uiTools/' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'DICOM/' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'Image3D/' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'polygons/' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'POLYLINE/' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'CardiacAnalysisTools/' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'OPTIM/' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'MatrixCalculus/' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'uiTools/OrbitPanZoom/' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'euclidean_distance_A2B/' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'drawContours/' ) );
clearvars('MPP_source_FOLDER');
enableVTK;
end

% permutation = load('/data/sim3d/ORd/permutation.txt'); %the nodes are renumbered by chaste!!
% 
% P=permutation+1; %matlab starts ordering from 1

mesh = read_CHASTE('/data/temp/chaste/HEART');

fibers = read_CHASTE('/data/temp/chaste/HEART.ortho');

write_VTK(mesh,'/data/temp/fibers.vtk');
write_VTK_UNSTRUCTURED_GRID( mesh.triORTHO , '/data/temp/fibers.vtk' )

mesh_vtk = read_VTK('/data/sim3d/louie/HEART_with_UPS.vtk'); 

fileinfo = hdf5info('/data/sim3d/ORd/results.h5');

upstroke = h5read('/data/sim3d/ORd/results.h5',fileinfo.GroupHierarchy.Datasets(3).Name);

upstroke_heart = upstroke(1,1:dummy); %select only the heart activation times!! 

%upstroke_heart = upstroke(P);

upstroke_heart = upstroke(P);

mesh_vtk.xyzATs=  transpose(upstroke_heart);

write_VTK(mesh_vtk,'/data/sim3d/ORd/louieORdATshighstim.vtk','binary');
%write_VTK(mesh_vtk,'/data/temp/louieORdATs_ordered.vtk','binary');
patch('Faces',mesh.tri,'Vertices','FaceVertexCData',upstroke_heart,mesh.xyz,'Facecolor','interp')
