
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

mesh = read_CHASTE('D:\ARVC meshing automatic\patients\patient02\chaste_ernesto\TORSO');

write_VTK(mesh,'D:\ARVC meshing automatic\patients\patient02\TORSO_vol.vtk','binary');

%find the rows which contain the RL electrode in the tetra indices
indices = find(mesh.tri(:,1) == 9463839);
indices2 =  find(mesh.tri(:,2) == 9463839);
indices3 = find(mesh.tri(:,3) == 9463839);
indices4 = find(mesh.tri(:,4) == 9463839);

all_indices = vertcat(indices, indices2,indices3,indices4);

all_mesh_indices =vertcat(mesh.tri(indices,:),mesh.tri(indices2,:),mesh.tri(indices3,:),mesh.tri(indices4,:));


all_mesh_indices_unique = unique(all_mesh_indices(:));

%find the coordinates of the unique nodes, which make up the tetrahedra
%touching the RL node
slab_xyz = mesh.xyz(all_mesh_indices_unique,:);

%reindex the nodes so that they can later be plotted, starting from 1
for i=1:size(all_mesh_indices_unique,1)
    all_mesh_indices(all_mesh_indices == all_mesh_indices_unique(i)) =i;
end

atts = mesh.triATTS(all_mesh_indices_unique);

%plot the slab containing all the tetra touching the RL node
figure()
hold on
patch('vertices', slab_xyz, 'faces',all_mesh_indices,'edgecolor','none')
camlight
alpha(0.3)
scatter3(slab_xyz(8,1),slab_xyz(8,2),slab_xyz(8,3)) %plot the point corresponding to the RL node (change this for each patient!!)
