
%%%%Windows

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
if ~isMPPdir( MPP_source_FOLDER ), MPP_source_FOLDER = 'Z:\Dropbox\shared\';                      end
if ~isMPPdir( MPP_source_FOLDER ), MPP_source_FOLDER = 'E:\Dropbox\Vigente\shared\';                                  end
clearvars('isMPPdir');
MPP_source_FOLDER = fullfile( MPP_source_FOLDER , 'MeshPersonalizationPipeline' );
end
try
set(0,'DefaultFigureCreateFcn','factory');
fprintf('',cellfun(@(p)isempty(strfind(p,'Dropbox'))||find(rmpath(p),1),strsplit(path,';','CollapseDelimiters',true)));
fprintf('',cellfun(@(p)isempty(strfind(p,MPP_source_FOLDER))||find(rmpath(p),1),strsplit(path,';','CollapseDelimiters',true)));
addpath(                      MPP_source_FOLDER );
addpath( fullfile(            MPP_source_FOLDER ,   'MPP_tools\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'thirdParty\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'thirdParty\export_fig' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'IO\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'LieAlgebra\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'MESH\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'MESHES\' ) ); 
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'MESHES\jigsaw\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'MESHES\tetgen\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'MESHES\carve\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'MESHES\gmsh\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'Tools\' ) );
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'Tools\parseargs' ) );
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
enableVTK;
end

addpath(genpath('/home/scratch/petnov/Dropbox'));

mesh = read_CHASTE('D:\ARVC meshing automatic\patients\patient02\chaste_ernesto\HEART');


control_mesh = read_VTK('D:\meshes\control cohort\DTI003_V_wall_tickness.vtk');

length = size(control_mesh.xyz,1);

add = ones(length,3);
add_long(:,1) = add(:,1) *-20;
add_long(:,2) = add(:,2)*0;
add_long(:,3) = add(:,3)*10;

add_long(:,1) = add(:,1) *-40;
add_long(:,2) = add(:,2)*60;
add_long(:,3) = add(:,3)*60;


control_mesh.xyz= control_mesh.xyz + add_long;


write_VTK(control_mesh,'D:\meshes\control cohort modified\DTI003_V_wall_tickness.vtk','binary');






