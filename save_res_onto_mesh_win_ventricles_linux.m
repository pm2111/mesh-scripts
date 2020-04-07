addpath(genpath('/home/scratch/petnov/Dropbox/shared - copy'));cd  /data/sim3d/TestHeart09LGEStrong/







%mesh_vtk = read_VTK('/data/sim3d/louie/HEART_with_UPS.vtk'); 
cd  /data/sim3d/TestHeart09LGEStrong/
mesh = read_CHASTE('/data/meshes/chaste09/HEART');
permutation = load('permutation.txt');
permutation = permutation+1;

%[EPI,LV,RV,~,~,~] = HEARTparts(mesh);
fid=fopen('voltages/voltages.bin855','r');
voltages=fread(fid,'double');
fileinfo =h5info('results.h5');
upstroke = h5read('results.h5',strcat('/',fileinfo.Datasets(3).Name));
% upstroke_heart = upstroke(1,:,1); %select only the heart activation times!! 
upstroke_heart = upstroke(1,permutation,1);
t=2;
output = h5read('results.h5','/Data',[1,1,t] ,[1,N,1]);
voltages=output(permutation);
N=size(mesh.xyz,1);
fid2= fopen("pat05_small_par_o5081914.out",'r');
fid3=fopen("pat05_small_par_o5081990.out",'r');
acts= dlmread("pat05_small_par_o5081990.out"," ");

acts_full =zeros(size(mesh.xyz,1),1);
acts_full(acts(:,1)+1)=acts(:,2);

%acts_full= acts_full(permutation);
%upstroke_epi = upstroke_heart(1:size(EPI.nodes,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileinfo = hdf5info('results.h5');

close all
figure()
%camlight
hold on
%htext =text(0,-0.6,sprintf('time=%f',t));
% 
%  for t=1:10:200
% output = h5read('results.h5','/Data',[1,1,t] ,[1,N,1]);
% %delete(htext)    
% %patch(mesh.xyz(:,1),mesh.xyz(:,2),mesh.xyz(:,3),upstroke_heart')
patch('Faces',mesh.tri,'Vertices',mesh.xyz,'FaceVertexCData',transpose(upstroke_heart),'FaceColor','interp','EdgeColor','none');
% patch('Faces',EPI.tri,'Vertices',EPI.xyz,'FaceVertexCData',transpose(upstroke_epi),'FaceColor','interp','EdgeColor','none');
% %patch('Faces',mesh.tri,'Vertices',mesh.xyz,'FaceVertexCData',distances_to_min_index,'FaceColor','interp','EdgeColor','none');
% %patch('Faces',mesh.tri,'Vertices',mesh.xyz,'FaceVertexCData',CV','FaceColor','interp','EdgeColor','none');
% %caxis([0 100])
% title('conduction velocity [ms]')
% set(gca,'FontSize',21)
% %htext =text(0,-0.6,sprintf('time=%f',t));
% %patch('Faces',mesh.face,'Vertices',mesh.xyz,'FaceVertexCData',transpose(upstroke(P)),'FaceColor','interp','EdgeColor','none');
% %scatter3(mesh.xyz(ids_unstable,1),mesh.xyz(ids_unstable,2),mesh.xyz(ids_unstable,3))
% colorbar
% %caxis([-95 -80])
% axis off
% pause
% end
% hisct(output(P'))

mesh.xyzvoltage1= voltages';

mesh.xyzUpstroke = transpose(upstroke_heart);

mesh.xyzActs=acts_full;
%compare results with control values for this patient
% control_res = '/data/sim3d/ARVC_005/TestHeart05FullBeatAndActivation/';
% upstroke_control = h5read(strcat(control_res,'results.h5'),fileinfo.GroupHierarchy.Datasets(3).Name);
% upstroke_control = upstroke_control(1,permutation,1);
% mesh.xyzUpstrokeControl = transpose(upstroke_control);
% mesh.xyzUpstrokeDifference = transpose(upstroke_heart- upstroke_control);
% upstroke_diff= transpose(upstroke_heart- upstroke_control);
mesh.celltype =10;
write_VTK(mesh,'results.vtk','binary');
