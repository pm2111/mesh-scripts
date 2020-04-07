%save the torso into chaste format
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
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'MESHES/' ) ); enableVTK;
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
end

torso = read_CHASTE('D:\ARVC meshing automatic\patients\patient02\TORSO');


electrodes = load('D:\ARVC meshing automatic\patients\patient02\mpp\ECG_ELECTRODES.mat');

E = [ electrodes.ELECTRODES.LA ; electrodes.ELECTRODES.RA ; electrodes.ELECTRODES.LL ; electrodes.ELECTRODES.RL ; electrodes.ELECTRODES.V1 ; electrodes.ELECTRODES.V2 ; electrodes.ELECTRODES.V3 ; electrodes.ELECTRODES.V4 ; electrodes.ELECTRODES.V5 ; electrodes.ELECTRODES.V6 ];
E=E/10; %chaste is in cm
figure()
patch('vertices',torso.xyz,'faces',torso.face,'edgecolor','none'); camlight
alpha(0.3)
hplot3d( E , '*m' ,'LineWidth',4,'MarkerSize',10 );
axis(objbounds)
heart_surf = read_VTK('D:\ARVC meshing automatic\patients\patient02\HEART_0.20.vtk');%+ 1
hold all
patch('vertices',heart_surf.xyz/10,'faces',heart_surf.tri,'edgecolor','none','Facecolor','r'); 
axis off
axis equal
camlight

write_VTK(torso,'/data/meshes/patient02/TORSO_cm_pat2.vtk');

S = Mesh(rand(1000,3),'convexhull');

V = gmsh(S, ['v']);
