close all;
clear all;
%Visualisation tool for the 3D Meshes from CMR imaging
%           -V.xyz are positions of nodes 
%           -V.face are elements of the body surface (triangles)
%           -V.heartface are the elements on the heart face (triangles)
%           -Can plot activation time on top of mesh

 set(gcf,'renderer','opengl')

activation = load('/home/scratch/Firefox_downloads/ResultsData.mat');

mesh = load('/home/scratch/Firefox_downloads/ARVC14Mesh.mat'); 

logical_index = logical(mesh.V.Heartface); 

%Below are ways to plot the surface nodes and edges
%patch('vertices', mesh.V.xyz, 'faces',mesh.V.Heartface +1,'facecolor', 'r','edgecolor','none')

% patch('vertices', mesh.V.xyz, 'faces',mesh.V.face+1,'facecolor', 'interp','cdata',activation.ATMap,'edgecolor','none')
% patch('vertices', mesh.V.xyz, 'faces',mesh.V.Heartface+1,'facecolor', 'interp','cdata',activation.ATMap,'edgecolor','none') %here we interpolate the results from the sims and plot them on top
% figure
% patch('vertices', mesh.V.xyz, 'faces',mesh.V.Heartface+1,'facecolor', 'interp','cdata',activation.ATMap,'edgecolor','none')
% camlight
% patch('vertices', mesh.V.xyz, 'faces',mesh.V.Heartface+1,'facecolor', 'interp','cdata',activation.APDMap,'edgecolor','none')%need +1 as indexing starting from zero in the mesh making software
% hold off
% patch('vertices', mesh.V.xyz, 'faces',mesh.V.Heartface+1,'facecolor', 'interp','cdata',activation.APDMap,'edgecolor','none')
 figure
 patch('vertices', mesh.V.xyz, 'faces',mesh.V.Heartface+1,'facecolor', 'interp','cdata',activation.APDMap,'edgecolor','none')
% caxis([250 350])
