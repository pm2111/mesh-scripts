%% Reset all
clear all;
close all;
clc;
addpath(genpath('/home/scratch/petnov/Dropox/shared - copy/IO/'));
addpath('/home/scratch/petnov/Dropbox/shared - copy/IO');

%% Move to working directory
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
% cd('new Tomasz');
%% Read the ECGs into 'info'
data_path = '/data/sim3d/Eikonal/chaste06_coarse/';
new_path = 'eikonal06';
if ~isfolder(new_path)
mkdir([data_path    new_path])
end
addpath(data_path)
%% Load coarse and fine meshes
fileName = 'HEART';
mesh=read_CHASTE(strcat(data_path,fileName));
fileName_final = 'eikonal06';


%% Calculate the edges of the fine mesh
elements = mesh.tri;
elements_1 = sort(elements(:, [1 2]), 2);
elements_2 = sort(elements(:, [1 3]), 2);
elements_3 = sort(elements(:, [1 4]), 2);
elements_4 = sort(elements(:, [2 3]), 2);
elements_5 = sort(elements(:, [2 4]), 2);
elements_6 = sort(elements(:, [3 4]), 2);
clear elements;
edges = unique([elements_1; elements_2; elements_3; elements_4; elements_5; elements_6], 'rows');
clear elements_1;
clear elements_2;
clear elements_3;
clear elements_4;
clear elements_5;
clear elements_6;
nodesA = mesh.xyz(edges(:, 1), :);
nodesB = mesh.xyz(edges(:, 2), :);
edgesCenter = [(nodesA(:, 1)+nodesB(:, 1))/2, (nodesA(:, 2)+nodesB(:, 2))/2, (nodesA(:, 3)+nodesB(:, 3))/2];
clear nodesA;
clear nodesB;
mesh.edges = edges;
mesh.edgeCenters = edgesCenter;
clear edges;
clear edgesCenter;
%% Calculate the centres of all the tetrahedrons in the fine mesh
fineTetrahedronCenters = (mesh.xyz(mesh.tri(:, 1), :)+ mesh.xyz(mesh.tri(:, 2), :) + mesh.xyz(mesh.tri(:, 3), :)+ mesh.xyz(mesh.tri(:, 4), :))./4;
mesh.centers = fineTetrahedronCenters;

%% Create fine mesh faces
rvNodes = dlmread([fileName '.rv'], ' ', 0, 0)+1;
lvNodes = dlmread([fileName '.lv'], ' ', 0, 0)+1;
epiNodes = dlmread([fileName '.epi'], ' ', 0, 0)+1;
fineFaces = dlmread([fileName '.face'], ' ', 1, 1)+1;

epiFaces = fineFaces(all(ismember(fineFaces(:, :),epiNodes), 2), :);
lvFaces = fineFaces(all(ismember(fineFaces(:, :),lvNodes), 2), :);
rvFaces = fineFaces(all(ismember(fineFaces(:, :),rvNodes), 2), :);

mesh.epiface = epiFaces;
mesh.lvface = lvFaces;
mesh.rvface = rvFaces;

%electrode positions
electrodes=load(strcat(data_path,'ECG_ELECTRODES.mat'));
el_names=fieldnames(electrodes.ELECTRODES);
pos=[];
for i=1:size(el_names,1)
    pos=vertcat(pos,electrodes.ELECTRODES.(el_names{i}));
end
pos=pos;
csvwrite(strcat(data_path,   new_path, '/' ,fileName_final, '_electrodePositions.csv'),pos/10);

%roots
%roots_pos=dlmread(strcat(data_path,'start_pos.txt'),',');
% 
% %check the plot below to see if correct roots are found
% figure()
% patch('vertices',mesh.xyz,'faces',mesh.face,'facecolor','b','edgecolor','none');
% camlight
% hold on;
% scatter3(mesh.xyz(roots,1),mesh.xyz(roots,2),mesh.xyz(roots,3),35,'filled','r')
% caxis([0 20])
% 

distance=@(x1,x2) ((x1(:,1)-x2(:,1)).^2+(x1(:,2)-x2(:,2)).^2+(x1(:,3)-x2(:,3)).^2).^0.5;

%roots_pos(5,:)=[6.463, -9.567 , -3.958];
% 
% roots=[];
% for i=1:7
%     [~,id]=min(distance(mesh.xyz,roots_pos(i,:)));
%     roots(i)=id;
% end
% 
% 
% 
% dlmwrite([data_path new_path '/' fileName_final '_fine_rootNodes5.csv'],roots([1:5,7]),'precision',10);
% C = {'k','b','r','g','y','m',[.8 .2 .6]}; % Cell array of colros.
% figure()
% patch('vertices',mesh.xyz,'faces',mesh.face,'facecolor','b','edgecolor','none');
% camlight
% hold on;
% for i=1:7
%     if i~=6
% scatter3(mesh.xyz(roots(i),1),mesh.xyz(roots(i),2),mesh.xyz(roots(i),3),205,[C{i}],'filled')
% camlight
%     end
% end



% 
% %lge regions
% fib=load(strcat(data_path,'regions_full_light.txt'));
% fib=fib(:,1);
% levels=unique(fib);
% scaling=[];
% for i=1:size(levels,1)
%     idx=find(fib==levels(end+1-i));
%     scaling(idx)=i-1;
% end
% 
% 
% dlmwrite(strcat(data_path,   new_path, '/', fileName_final, '_fine_nodeScaling.csv'),scaling,'precision',10);






%% To csv for the Python code
% Fine
dlmwrite(strcat(data_path,   new_path, '/', fileName_final, '_coarse_tetrahedronCenters.csv'),mesh.centers,'precision',10);
dlmwrite( strcat(data_path,   new_path, '/', fileName_final, '_coarse_xyz.csv'),mesh.xyz,'precision',10);
dlmwrite(strcat(data_path,   new_path, '/', fileName_final, '_coarse_tri.csv'),mesh.tri,'precision',10);
dlmwrite(strcat(data_path,   new_path, '/', fileName_final, '_coarse_edgesCenter.csv'),mesh.edgeCenters,'precision',10);
dlmwrite(strcat(data_path,   new_path, '/', fileName_final, '_coarse_epiface.csv'),mesh.epiface,'precision',10);
dlmwrite( strcat(data_path,   new_path, '/', fileName_final, '_coarse_lvface.csv'),mesh.lvface,'precision',10);
dlmwrite(strcat(data_path,   new_path, '/', fileName_final, '_coarse_rvface.csv'),mesh.rvface,'precision',10);
dlmwrite( strcat(data_path,   new_path, '/', fileName_final,'_coarse_edges.csv'),mesh.edges,'precision',10);

%fibres - need to copy the original fibers on the amazon cluster and run
%the fiber script on there
% EOF