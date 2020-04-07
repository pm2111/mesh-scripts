addpath(genpath('C:\users\petnov\Dropbox\'));

cd D:\ARVC' meshing automatic'\patients\patient05\




%mesh_vtk = read_VTK('/data/sim3d/louie/HEART_with_UPS.vtk'); 
cd  D:\ARVC' meshing 'automatic\patients\patient05\results\TestHeart05smallpara
mesh = read_CHASTE('..\..\chaste_tiny\HEART');

fileinfo = hdf5info('results.h5');
permutation = load('D:\ARVC meshing automatic\patients\patient05\results\TestHeart05FullBeatandActivationScar\permutation.txt');
permutation=load('permutation.txt');
permutation = permutation+1;
upstroke = h5read('results.h5',fileinfo.GroupHierarchy.Datasets(3).Name);

upstroke_heart1 = upstroke(1,permutation,1); %select only the heart activation times!! 
upstroke_heart2 = upstroke(1,permutation,2); %select only the heart activation times!! 

N=size(mesh.xyz,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileinfo = hdf5info('results.h5');

close all
figure()
%camlight
hold on
%htext =text(0,-0.6,sprintf('time=%f',t));

 for t=1:7
output = h5read('results.h5','/Data',[1,1,t] ,[1,N,2]);
%delete(htext)    
%patch(mesh.xyz(:,1),mesh.xyz(:,2),mesh.xyz(:,3),upstroke_heart')
patch('Faces',mesh.tri,'Vertices',mesh.xyz,'FaceVertexCData',transpose(upstroke_heart1),'FaceColor','interp','EdgeColor','none');
patch('Faces',mesh.tri,'Vertices',mesh.xyz,'FaceVertexCData',transpose(output(permutation')),'FaceColor','interp','EdgeColor','none');
%patch('Faces',mesh.tri,'Vertices',mesh.xyz,'FaceVertexCData',distances_to_min_index,'FaceColor','interp','EdgeColor','none');
%patch('Faces',mesh.tri,'Vertices',mesh.xyz,'FaceVertexCData',CV','FaceColor','interp','EdgeColor','none');
%caxis([0 100])
title('conduction velocity [ms]')
set(gca,'FontSize',21)
%htext =text(0,-0.6,sprintf('time=%f',t));
%patch('Faces',mesh.face,'Vertices',mesh.xyz,'FaceVertexCData',transpose(upstroke(P)),'FaceColor','interp','EdgeColor','none');
%scatter3(mesh.xyz(ids_unstable,1),mesh.xyz(ids_unstable,2),mesh.xyz(ids_unstable,3))
colorbar
headlight
%caxis([-95 -80])
axis off
end
hisct(output(P'))


mesh.xyzUpstroke = transpose(upstroke_heart2);
%compare results with control values for this patient
% control_res = 'D:\ARVC meshing automatic\patients\patient05\results\TestHeart05FullBeatandActivation\';
% upstroke_control = h5read(strcat(control_res,'results.h5'),fileinfo.GroupHierarchy.Datasets(3).Name);
% upstroke_control = upstroke_control(1,permutation,1);
% mesh.xyzUpstrokeControl = transpose(upstroke_control);
% mesh.xyzUpstroke2 = transpose(upstroke_heart2);
% 
% mesh.xyzUpstrokeDifference = transpose(upstroke_heart- upstroke_control);
% upstroke_diff= transpose(upstroke_heart- upstroke_control);
mesh.celltype =10;
write_VTK(mesh,'results.vtk','binary');
