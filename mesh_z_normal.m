clear all
%reorient mesh so that normal to base plane is in the z direction

addpath('/users/petnov/Dropbox/shared/IO');
addpath('/users/petnov/Dropbox/mesh scripts/');
addpath('/users/petnov/Dropbox/shared/MESHES/');

%% Setting Up General settings!
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


%% import our ventricular surfaces

mesh = read_VTK('/home/scratch/meshes ready/patient02/HEART_0.04.vtk');

[EPI,LV,RV,LID,MAT] = HEARTparts(mesh);



mesh_normal = transform(mesh,MAT);

write_VTK_UNSTRUCTURED_GRID(mesh,'/home/scratch/meshes ready/patient02/HEART_0.04_normal.vtk')