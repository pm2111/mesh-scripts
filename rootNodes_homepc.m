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
if ~isMPPdir( MPP_source_FOLDER ), MPP_source_FOLDER = 'E:\Dropbox\Vigente\shared\'; 
if ~isMPPdir( MPP_source_FOLDER ), MPP_source_FOLDER = 'C:\Users\peter\Dropbox\shared\'; end
end
clear('isMPPdir');
MPP_source_FOLDER = fullfile( MPP_source_FOLDER , 'MeshPersonalizationPipeline' );
% end

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
addpath( fullfile( fileparts( MPP_source_FOLDER ) , 'MESHES\' ) ); enableVTK;
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
end
fprintf('MPP_source_FOLDER is :   "%s"\n' , MPP_source_FOLDER );
TORSO_MODEL_DIR = fullfile(MPP_source_FOLDER,'TORSO');
fprintf('TORSO_FOLDER is :        "%s"\n' , TORSO_MODEL_DIR );
clear('MPP_source_FOLDER');


TORSO_MODEL_DIR = 'C:\Dphil\patient02\chaste_ernesto\';
H0 = read_VTK( fullfile( TORSO_MODEL_DIR , 'HEART_with_UPS_louie.vtk' ) ); %this is louie's heart!

LV0_roots = [2012680,623088,983785,2359613];
RV0_roots = [1396217, 1226484, 2333921];



figure()
plotMESH( H0 , 'ne' );
hplot3d( H0.xyz( LV0_roots ,:) , 'o1kb10' );
hplot3d( H0.xyz( RV0_roots ,:) , 'o1kr10' );
headlight


%%

H1 = read_VTK(  'C:\Dphil\patient02\chaste_ernesto\HEART.vtk'); %our target heart, perhaps need full torso??

plotMESH( H0 ,'ne');
hplotMESH( H1 , 'ne' ,'FaceColor',[1 1 1]*0.6)
headlight; axis(objbounds);


%%

u = loadv( 'C:\Dphil\patient02\chaste_ernesto\Hmapping' , 'u' ); %transfer matrix u
%this is a function converting coordinates from H0 to H1

if 0
HD = H0;
HD.xyz = u( HD.xyz );

plotMESH( HD ,'ne');
hplotMESH( H1 , 'ne' ,'FaceColor',[1 1 1]*0.6)
headlight; axis(objbounds);
end


%%

LV1_roots = vtkClosestPoint( Mesh(H1) , u( H0.xyz( LV0_roots ,:) ) );
RV1_roots = vtkClosestPoint( Mesh(H1) , u( H0.xyz( RV0_roots ,:) ) );

figure()
plotMESH( H1 , 'ne' );
hplot3d( H1.xyz( LV1_roots ,:) , 'o1kb10' );
hplot3d( H1.xyz( RV1_roots ,:) , 'o1kr10' );
headlight;