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



mesh = read_CHASTE('D:\meshes\data\2Dmeshcm');

%mesh_vtk = read_VTK('/data/sim3d/louie/HEART_with_UPS.vtk'); 
cd  D:\ARVC' meshing automatic'\patients\patient05\results\annulus_indexed_fibrosis_new_stim\
fileinfo = hdf5info('results.h5');
permutation = load('permutation.txt');
distances_to_min_index = load('distances.txt');
permutation = permutation+1;
upstroke = h5read('results.h5',fileinfo.GroupHierarchy.Datasets(7).Name);
upstroke_heart = upstroke(1,:,1); %select only the heart activation times!! 
upstroke_heart = upstroke_heart(permutation);
N=size(mesh.xyz,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fist_act_time,first_index_to_activate] = min(upstroke_heart);
for k=1:size(mesh.xyz,1)
%     [row,col]=find(mesh.tri ==k);
%     target_index = mesh.tri(row(1),col(1));
% 
%     indices_neighbors=[];
%     for i=1:size(row,1)
%         indices_neighbors = [indices_neighbors,mesh.tri(row(i),1),mesh.tri(row(i),2),mesh.tri(row(i),3)];
% 
%     end
% % 
%     unique_neighbors = unique(indices_neighbors);
%     eliminate = find(unique_neighbors ==target_index);
%     unique_neighbor_indices = unique_neighbors([1:eliminate-1 eliminate+1:end]);
% 
%     unique_neighbor_upstroke_times = upstroke_heart(unique_neighbor_indices);

    CV(k)= distances_to_min_index(k)*1000/(upstroke_heart(k)- upstroke_heart(first_index_to_activate));%cm/s
    
    %     larger_than_zero = find(CVs<0);
%     try
%         CV(k) =max(CVs(larger_than_zero));
%     catch
%         CV(k) =0;
%     end
end
    CV(find(CV==inf)) =0.0;
fileinfo = hdf5info('results.h5');

close all
figure()
%camlight
hold on
%htext =text(0,-0.6,sprintf('time=%f',t));

 for t=1:10:1000
output = h5read('results.h5','/Data',[1,1,t] ,[1,N,1]);
%delete(htext)    
%patch(mesh.xyz(:,1),mesh.xyz(:,2),mesh.xyz(:,3),upstroke_heart')
patch('Faces',mesh.tri,'Vertices',mesh.xyz,'FaceVertexCData',transpose(upstroke_heart),'FaceColor','interp','EdgeColor','none');
patch('Faces',mesh.tri,'Vertices',mesh.xyz,'FaceVertexCData',transpose(output(permutation')),'FaceColor','interp','EdgeColor','none');
patch('Faces',mesh.tri,'Vertices',mesh.xyz,'FaceVertexCData',distances_to_min_index,'FaceColor','interp','EdgeColor','none');
%patch('Faces',mesh.tri,'Vertices',mesh.xyz,'FaceVertexCData',CV','FaceColor','interp','EdgeColor','none');
caxis([0 100])
title('conduction velocity [ms]')
set(gca,'FontSize',21)
%htext =text(0,-0.6,sprintf('time=%f',t));
%patch('Faces',mesh.face,'Vertices',mesh.xyz,'FaceVertexCData',transpose(upstroke(P)),'FaceColor','interp','EdgeColor','none');
%scatter3(mesh.xyz(ids_unstable,1),mesh.xyz(ids_unstable,2),mesh.xyz(ids_unstable,3))
colorbar
%caxis([-95 -80])
axis off
pause
end
hisct(output(P'))


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

