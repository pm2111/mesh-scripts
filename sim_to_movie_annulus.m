%%%% ECGi to movie%%%%%%%%%%%

% need to load the arvc pat data(here pat 14 & baseline)
close all
clear all

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

cd  D:\ARVC' meshing automatic'\patients\patient05\TestHeart05Roots\
permutation = load('permutation.txt'); %the nodes are renumbered by chaste!!
P=permutation+1; %matlab starts ordering from 1
mesh = read_CHASTE('../chaste_gmsh/HEART');
N =size(permutation,1);
nr_frames = 200;
h = figure();
cd D:\ARVC' meshing automatic'\patients\patient05\results\annulus_indexed_fibrosis_two_node\
fileinfo = hdf5info('results.h5');
for i=1:nr_frames
    output = h5read('results.h5','/Data',[1,1,i] ,[1,N,1]);
    %set(h, 'Visible', 'off');
    patch('Faces',mesh.tri,'Vertices',mesh.xyz,'FaceVertexCData',transpose(output),'FaceColor','interp','EdgeColor','none');
    axis off
    htext =text(0,-0.6,sprintf('time=%f',i));

    F(i) = getframe(h);
    delete(htext);
end


fig = figure;
movie(fig,F,1)

v = VideoWriter('D:\sim3d\annulus_indexed_fibrosis_1node','Motion JPEG AVI');

open(v)
writeVideo(v,F);

close(v)

