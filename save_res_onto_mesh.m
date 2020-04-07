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

permutation = load('permutation.txt'); %the nodes are renumbered by chaste!!

P=permutation+1; %matlab starts ordering from 1

mesh = read_CHASTE('chaste_tiny/HEART');

%mesh_vtk = read_VTK('/data/sim3d/louie/HEART_with_UPS.vtk'); 

fileinfo = hdf5info('results.h5');

upstroke = h5read('results.h5',fileinfo.GroupHierarchy.Datasets(3).Name);
N=31083;
t=7;
output = h5read('results.h5','/Data',[1,1,t] ,[1,N,1]);

upstroke_heart = upstroke(1,:,1); %select only the heart activation times!! 

%upstroke_heart = upstroke(P);

upstroke_heart = upstroke(P);

mesh.xyzATs=  upstroke_heart;

write_VTK(mesh,'activation_small_mesh.vtk','binary');
%write_VTK(mesh_vtk,'/data/temp/louieORdATs_ordered.vtk','binary');
 figure()
 
 for t=1:5:100
    output = h5read('results.h5','/Data',[1,1,t] ,[1,N,1]);
   
%patch(mesh.xyz(:,1),mesh.xyz(:,2),mesh.xyz(:,3),upstroke_heart)
patch('Faces',mesh.face,'Vertices',mesh.xyz,'FaceVertexCData',transpose(output(P)),'FaceColor','interp');
caxis([-85 100])

colorbar
pause
end
hist(output(P'))