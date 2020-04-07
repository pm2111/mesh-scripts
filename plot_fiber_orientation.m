MPP_version = 'arvc:1.0';


try
isMPPdir = @(d) isdir(d) && isdir( fullfile( d , 'MeshPersonalizationPipeline' ) );
MPP_source_FOLDER = '';
if ~isMPPdir( MPP_source_FOLDER ), MPP_source_FOLDER = fileparts(fileparts(which('MeshPersonalizationPipeline.m')));  end
if ~isMPPdir( MPP_source_FOLDER ), MPP_source_FOLDER = fileparts(fileparts(mfilename('fullpath')));                   end
if ~isMPPdir( MPP_source_FOLDER ), MPP_source_FOLDER = fileparts(mfilename('fullpath'));                              end
if ~isMPPdir( MPP_source_FOLDER ), MPP_source_FOLDER = pwd;                                                           end
if ~isMPPdir( MPP_source_FOLDER ), MPP_source_FOLDER = fileparts(pwd);                                                end
if ~isMPPdir( MPP_source_FOLDER ), MPP_source_FOLDER = 'D:\shared\';                                          end
clearvars('isMPPdir');
MPP_source_FOLDER = fullfile( MPP_source_FOLDER , 'MeshPersonalizationPipeline' );
end

fprintf('MPP_source_FOLDER is :   "%s"\n' , MPP_source_FOLDER );

TORSO_MODEL_DIR = fullfile(MPP_source_FOLDER,'TORSO');
fprintf('TORSO_FOLDER is :        "%s"\n' , TORSO_MODEL_DIR );

try
try,OrbitPanZoom('never');end
ww__=warning;warning('off','all');rmpath(genpath('E:\Dropbox\'));warning(ww__);clearvars('ww__');
ww__=warning;warning('off','all');rmpath(genpath(MPP_source_FOLDER));warning(ww__);clearvars('ww__');

addpath(           MPP_source_FOLDER );
addpath( fullfile( MPP_source_FOLDER , 'MPP_tools\' ) );

addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'thirdParty\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'thirdParty\export_fig' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'IO\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'LieAlgebra\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'MESH\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'MESHES\' ) ); enableVTK;
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'MESHES\jigsaw\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'MESHES\tetgen\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'MESHES\carve\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'Tools\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'uiTools\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'DICOM\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'Image3D\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'polygons\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'POLYLINE\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'CardiacAnalysisTools\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'OPTIM\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'MatrixCalculus\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'uiTools\OrbitPanZoom\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'euclidean_distance_A2B\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'drawContours\' ) );
addpath d:\shared\Tools\parseargs\; %important, this was missed last time!!!
clearvars('MPP_source_FOLDER');
end

cd D:\ARVC' meshing automatic'\patients\patient05\chaste_tiny\
mesh = read_CHASTE('HEART');
apexbase = dlmread('apexbase');

fibers = dlmread('HEART.ortho','\t',[1 0 size(mesh.tri,1) 8]);
permutation = load('permutation.txt'); %the nodes are renumbered by chaste!!

cd ..

P=permutation+1; %matlab starts ordering from 1



principal_fibers = fibers(:,1:3);

apexbase_vec = apexbase(1,:) - apexbase(2,:);

apexbase_vec = apexbase_vec /norm(apexbase_vec);
dot_prod = mtimes(principal_fibers , apexbase_vec');

angle = acos(dot_prod)*180/pi;
angle_normalised = angle/max(angle);
colors = zeros(size(angle_normalised,1),3);
colors(:,1) = angle_normalised;

fiber_nodes = ones(size(mesh.xyz,1),1);
fiber_nodes(mesh.tri(:,1)) = angle(mesh.tri(:,1));
fiber_nodes(mesh.tri(:,2)) = angle(mesh.tri(:,2));
fiber_nodes(mesh.tri(:,3)) = angle(mesh.tri(:,3));
fiber_nodes(mesh.tri(:,4)) = angle(mesh.tri(:,4));

%need 1 datapoint per vertex!! find a way to assign these efficiently!

%write_VTK(mesh_vtk,'/data/temp/louieORdATs_ordered.vtk','binary');
 figure()
 camlight
patch('Faces',mesh.face,'Vertices',mesh.xyz,'FaceVertexCData',fiber_nodes(P),'FaceColor','interp','EdgeColor','none');
%patch('Faces',mesh.face,'Vertices',mesh.xyz,'FaceColor',colors,'EdgeColor','none');

caxis([-85 100])

colorbar

hist(output(P'))