addpath(genpath('C:\users\petnov\Dropbox\'));/data/ecgi

cd D:\ARVC' meshing automatic'\patients\patient18\




%mesh_vtk = read_VTK('/data/sim3d/louie/HEART_with_UPS.vtk'); 
mesh = read_CHASTE('chaste_gmsh\HEART');

fileinfo = hdf5info('results\TestHeart18LVRVLGE\results_maps_mono.h5');
permutation = load('results\TestHeart18LVRVLGE\permutation.txt');
permutation = permutation+1;
upstroke = h5read('results\TestHeart18LVRVLGE\results_maps_mono.h5',fileinfo.GroupHierarchy.Datasets(1).Name);
upstroke_bi = h5read('results\TestHeart18LVRVLGE\results_maps.h5',fileinfo.GroupHierarchy.Datasets(1).Name);

upstroke_control=h5read('results\TestHeart18ControlECG\results_maps.h5',fileinfo.GroupHierarchy.Datasets(1).Name);

upstroke_heart2 = upstroke(1,permutation,2)-800; %select only the heart activation times!! 
upstroke_heart2bi = upstroke_bi(1,permutation,2)-800; %select only the heart activation times!! 

upstroke_control_heart2 = upstroke_control(1,permutation,1); %select only the heart activation times!! 

N=size(mesh.xyz,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fileinfo = hdf5info('results.h5');
t=3;
close all
%figure()
%camlight
% hold on
% %htext =text(0,-0.6,sprintf('time=%f',t));
% 
%  for t=1:7
% output = h5read('results_maps.h5','/Data',[1,1,t] ,[1,N,1]);
% %delete(htext)    
% %patch(mesh.xyz(:,1),mesh.xyz(:,2),mesh.xyz(:,3),upstroke_heart')
 %patch('Faces',mesh.tri,'Vertices',mesh.xyz,'FaceVertexCData',transpose(upstroke_heart2),'FaceColor','interp','EdgeColor','none');
% patch('Faces',mesh.tri,'Vertices',mesh.xyz,'FaceVertexCData',transpose(output(permutation')),'FaceColor','interp','EdgeColor','none');
% %patch('Faces',mesh.tri,'Vertices',mesh.xyz,'FaceVertexCData',distances_to_min_index,'FaceColor','interp','EdgeColor','none');
% %patch('Faces',mesh.tri,'Vertices',mesh.xyz,'FaceVertexCData',CV','FaceColor','interp','EdgeColor','none');
% %caxis([0 100])
% title('conduction velocity [ms]')
% set(gca,'FontSize',21)
% %htext =text(0,-0.6,sprintf('time=%f',t));
% %patch('Faces',mesh.face,'Vertices',mesh.xyz,'FaceVertexCData',transpose(upstroke(P)),'FaceColor','interp','EdgeColor','none');
% %scatter3(mesh.xyz(ids_unstable,1),mesh.xyz(ids_unstable,2),mesh.xyz(ids_unstable,3))
% colorbar
% headlight
% %caxis([-95 -80])
% axis off
% end
% hist(output(P'))


mesh.xyzUpstrokeLGE1 = transpose(upstroke_heart2);

mesh.xyzUpstrokeControl = transpose(upstroke_control_heart2);
%upstroke_diff= transpose(upstroke_heart2- upstroke_control);
mesh.celltype =10;
write_VTK(mesh,'results_roots_base.vtk','binary');
