%clear all

%%%%Linux
addpath('/users/petnov/Dropbox/shared/IO');
addpath('/users/petnov/Dropbox/mesh scripts/');

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



%% import our ventricular surfaces

%mesh = read_VTK('/home/scratch/meshes ready/patient02/HEART_0.04.vtk');
%%%%IN LINUX
mesh = read_CHASTE('/data/meshes/patient02/HEART');

filename = '/home/scratch/meshes ready/patient02/regions17/';

region_names = dir('/home/scratch/meshes ready/patient02/regions17/');

%%IN WINDOWS

mesh = read_VTK('D:\ARVC meshing automatic\patients\patient02\HEART.vtk');


region_names = dir('D:\ARVC meshing automatic\patients\patient02\regions17\');

filename = 'D:\ARVC meshing automatic\patients\patient02\regions17\';



regions = cell(1);
regions_complete = cell(1);
%rearrange all regions into one matrix
labelled_all = [];
regions_known_list_idx = [];
regions_known_list_number = [];
%have all labelled elements in a list, with column 1 indicating their index
%and column 2 indicating their region number
for i=1:size(region_names,1)-2
    labelled_all = vertcat(labelled_all,csvread(strcat(filename,region_names(i+2).name),1));
    regions{i} =  csvread(strcat(filename,region_names(i+2).name),1);
    regions_complete{i} = csvread(strcat(filename,region_names(i+2).name),1);
    region_sizes{i} = size(regions{i},1);
    regions_known_list_idx = vertcat(regions_known_list_idx,regions_complete{1,i}(:,1));
    length = size(regions_complete{1,i},1);
    regions_known_list_number = vertcat(regions_known_list_number,ones(length,1)*i);
end


region_start_idx{1} = 1;
for i=2:size(region_names,1)-2
    region_start_idx{i} = region_sizes{i} +region_sizes{i-1};
end
tic

total_dummy=0;
uncategorised_indices = [];
uncategorised_positions = [];
for i=1:size(mesh.xyz,1)
    
    %check if its in labelled data
    dummy=0;
    idx = find(labelled_all(:,1) == i);
    if isempty(idx)==0
        dummy =1;
    end
    
    if dummy ==0
        uncategorised_positions = vertcat(uncategorised_positions, mesh.xyz(i,1:3));
        uncategorised_indices = vertcat(uncategorised_indices,i);
    end
    total_dummy = total_dummy +dummy;
end



[idx,dist] = knnsearch(labelled_all(:,2:4),mesh.xyz(uncategorised_indices,1:3)); %find the index in labelled_all which is closesst to point in uncategorised indices 

indices_original_mesh = labelled_all(idx); %get indices of original mesh

%find the label of the labelled data which are closest to the unlabelled
%points
%region_full = zeros(size(indices_original_mesh));

for i =1:size(indices_original_mesh)
    
        belongs_to =  find(regions_known_list_idx ==indices_original_mesh(i));
        group = regions_known_list_number(belongs_to);
        regions_known_list_idx(end+1) = uncategorised_indices(i);
        regions_known_list_number(end+1) = group;
        
    

%     
%     regions_complete{1,region_id}(end+1,1) = uncategorised_indices(i);
%     regions_complete{1,region_id}(end,2) = uncategorised_positions(i,1);
%     regions_complete{1,region_id}(end,3) = uncategorised_positions(i,2);
%     regions_complete{1,region_id}(end,4) = uncategorised_positions(i,3);    
end

save(region_full,'D:\ARVC meshing automatic\patients\patient02\list of regions volumetric.mat');

save(strcat(filename(1:end-10),'patient aha17 regions LV and 12 regions RV volumetrised.mat'), 'regions_complete');

toc
% for j=1:size(regions_complete,2)
%     dummy = [];
%     %find triangles belonging to the jth region, having the ith element present
%     for i=1:size(regions_complete{1,j})
%         triangle_row = find(mesh.tri(:,1) == regions_complete{1,j}(i,1));
%         dummy = vertcat(dummy, mesh.tri(triangle_row,1:3));
%     end
%     triangles_region{1,j} = dummy; 
% end
%plot the resulting aha plot 
 colors = {'r','g','b','y','m','c','w','k','r','g','b','y','m','c','w','k','r','g','b','y','m','c','w','k','r','g','b','y','m','c','w','k','r','g','b','y','m','c','w','k'};
 figure()
 hold on
 axis off
     patch('vertices', mesh.xyz, 'faces', mesh.tri, 'facecolor', colors{3}, 'edgecolor', 'none') ;
     camlight 
     alpha(0.3)
     for i=1:size(regions_complete,2)
      scatter3(regions_complete{1,i}(:,2),regions_complete{1,i}(:,3),regions_complete{1,i}(:,4),colors{i});
     end
 j =1;
     %patch('vertices', mesh.xyz, 'faces', triangles_region{1,j}, 'facecolor', colors{j}, 'edgecolor', 'none') ;
j=2;
% 
% for i=1:size(regions,2)
% figure()
%  scatter3(regions_complete{1,i}(:,2),regions_complete{1,i}(:,3),regions_complete{1,i}(:,4),colors{i},'filled');
% end
%  %patch(regions_complete{1,i}(:,2),regions_complete{1,i}(:,3),regions_complete{1,i}(:,4),colors{i});
