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



cd D:\ARVC' meshing automatic'\patients\patient05\



mesh = read_CHASTE('chaste_gmsh/HEART');
cd  D:\ARVC' meshing automatic'\patients\patient05\TestHeart05Roots\
permutation = load('permutation.txt'); %the nodes are renumbered by chaste!!
P=permutation+1; %matlab starts ordering from 1

%mesh_vtk = read_VTK('/data/sim3d/louie/HEART_with_UPS.vtk'); 
cd  D:\ARVC' meshing automatic'\patients\patient05\results\TestHeart05Full\
fileinfo = hdf5info('results.h5');

upstroke = h5read('results.h5',fileinfo.GroupHierarchy.Datasets(3).Name);
N=size(permutation,1);
t=7;
output = h5read('results.h5','/Data',[1,1,t] ,[1,N,1]);

upstroke_heart = upstroke(1,:,1); %select only the heart activation times!! 

%upstroke_heart = upstroke(P);

upstroke_heart = upstroke(P);

mesh.xyzATs=  upstroke_heart;

write_VTK(mesh,'activation_small_mesh.vtk','binary');
%write_VTK(mesh_vtk,'/data/temp/louieORdATs_ordered.vtk','binary');
 figure()
 camlight
 hold all

 for t=1:20:200
    output = h5read('results.h5','/Data',[1,1,t] ,[1,N,1]);
       
%patch(mesh.xyz(:,1),mesh.xyz(:,2),mesh.xyz(:,3),upstroke_heart)
patch('Faces',mesh.face,'Vertices',mesh.xyz,'FaceVertexCData',transpose(output(P)),'FaceColor','interp','EdgeColor','none');
%patch('Faces',mesh.face,'Vertices',mesh.xyz,'FaceVertexCData',transpose(upstroke(P)),'FaceColor','interp','EdgeColor','none');
%scatter3(mesh.xyz(ids_unstable,1),mesh.xyz(ids_unstable,2),mesh.xyz(ids_unstable,3))
colorbar
%caxis([-95 -80])
axis off
colorbar
pause
end
hist(output(P'))


%investigate min voltage instability
ids_unstable = find(output < -90);
[unstable_tri_idx_row,unstable_tri_idx_col] = arrayfun(@(x)find(mesh.tri ==x,1),ids_unstable);

euclidian_dist = [];
for i=1:size(unstable_tri_idx,2)
    if unstable_tri_idx_col(i)>1
        euclidian_dist(i) = norm([ mesh.xyz(mesh.tri(unstable_tri_idx_row(i),unstable_tri_idx_col(i)),1), mesh.xyz(mesh.tri(unstable_tri_idx_row(i),unstable_tri_idx_col(i)),2), mesh.xyz(mesh.tri(unstable_tri_idx_row(i),unstable_tri_idx_col(i)),3)]- [ mesh.xyz(mesh.tri(unstable_tri_idx_row(i),unstable_tri_idx_col(i)-1) ,1), mesh.xyz(mesh.tri(unstable_tri_idx_row(i),unstable_tri_idx_col(i)-1),2), mesh.xyz(mesh.tri(unstable_tri_idx_row(i),unstable_tri_idx_col(i)-1),3)]);
    else 
        euclidian_dist(i) = norm([ mesh.xyz(mesh.tri(unstable_tri_idx_row(i),unstable_tri_idx_col(i)),1), mesh.xyz(mesh.tri(unstable_tri_idx_row(i),unstable_tri_idx_col(i)),2), mesh.xyz(mesh.tri(unstable_tri_idx_row(i),unstable_tri_idx_col(i)),3)]- [ mesh.xyz(mesh.tri(unstable_tri_idx_row(i),unstable_tri_idx_col(i)+1) ,1), mesh.xyz(mesh.tri(unstable_tri_idx_row(i),unstable_tri_idx_col(i)+1),2), mesh.xyz(mesh.tri(unstable_tri_idx_row(i),unstable_tri_idx_col(i)+1),3)]);
    end
    
end

