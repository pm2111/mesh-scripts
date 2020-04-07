addpath(genpath('C:\Users\peter\Dropbox\'));


%mesh_vtk = read_VTK('/data/sim3d/louie/HEART_with_UPS.vtk'); 
cd  D:\ARVC' meshing 'automatic\patients\patient05\results\TestHeart05Small_fib
mesh = read_CHASTE('/data/sim3d/ARVC_006/HEART');

%fileinfo = hdf5info('C:\Dphil\patient06\voltage.h5');
fid=fopen('C:\Dphil\patient06\voltages.bin805','r');
voltages=fread(fid,'double');

[EPI,LV,RV,~,~,~]=HEARTparts(mesh);

permutation = load('C:\Dphil\patient06\permutation.txt');
permutation = permutation+1;
voltages =voltages(permutation);


N=size(mesh.xyz,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileinfo = hdf5info('results.h5');

close all
figure()
%camlight
hold on
%htext =text(0,-0.6,sprintf('time=%f',t));

 %delete(htext)    
%patch(mesh.xyz(:,1),mesh.xyz(:,2),mesh.xyz(:,3),upstroke_heart')
patch('Faces',LV.tri,'Vertices',LV.xyz,'FaceVertexCData',voltages(1:87603),'FaceColor','interp','EdgeColor','none');
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
