addpath(genpath('C:\users\petnov\Dropbox\'));

cd D:\ARVC' meshing automatic'\patients\patient05\




%mesh_vtk = read_VTK('/data/sim3d/louie/HEART_with_UPS.vtk'); 
cd  D:\ARVC' meshing 'automatic\patients\patient09\results\TestHeart09Roots1111111GadLV
mesh = read_CHASTE('..\..\chaste_gmsh\HEART');

fileinfo = hdf5info('results_maps.h5');
permutation = load('D:\ARVC meshing automatic\patients\patient09\results\TestControlECGPat09\permutation.txt');
permutation = permutation+1;
upstroke_lge = h5read('D:\ARVC meshing automatic\patients\patient09\results\TestHeart09Roots1111111GadLV\results_maps.h5',fileinfo.GroupHierarchy.Datasets(1).Name);
upstroke= h5read('results.h5',fileinfo.GroupHierarchy.Datasets(1).Name);
upstroke_heart1 = upstroke(1,permutation,1); %select only the heart activation times!! 
upstroke_heart2 = upstroke(1,permutation,2); %select only the heart activation times!! 
% 
upstroke_lge2=upstroke_lge(1,permutation,2);
% ind_fake= find(upstroke_heart2<0);
% upstroke_heart2(ind_fake) =10000;

inds = load('D:\ARVC meshing automatic\patients\patient09\regions_pat09_in_mesh_vtk_ordering.txt');
epi_inds= load('D:\ARVC meshing automatic\patients\patient09\chaste_gmsh\HEART.epi')+1;

epi_cell_labels=load('D:\ARVC meshing automatic\patients\patient09\epi_cell_property.txt');

[EPI1,LV1,RV1,~,MAT1] = HEARTparts( mesh );


upstroke_epicardium = upstroke_heart1(epi_inds)-800;
upstroke_epicardium_lge=upstroke_lge2(epi_inds)-800;

inds_epi = inds(1:max(epi_inds));
coords=[]
for i=2:size(unique(inds),1)
    epi_and_region=find(epi_cell_labels==i);  %possibly the inds array is not in correct order
    [earliest_activation(i),ind]=min(upstroke_heart1(epi_and_region));
    local_coords=EPI1.xyz(epi_and_region,:);
    coords(i,1:3)=local_coords(ind,1:3);
end


%[reducedFaces,reducedVertices] = reducepatch(mesh.tri,mesh.xyz,300000);

dist = (coords(:,1).^2+coords(:,2).^2+coords(:,3).^2).^0.5;


figure()
p=patch( 'vertices', EPI1.xyz , 'faces', EPI1.tri  ,'facecolor','interp','cdata',upstroke_epicardium+800', 'edgecolor','none' ); 
hold on
[h,v] = IsoLine({EPI1.tri, EPI1.xyz},upstroke_epicardium,20,'k');
caxis([0 65])

axis off
axis equal; view(3); headlight
c=colorbar;
ylabel(c,'Upstroke Times [ms]','Fontsize',21);

set(gca,'Fontsize',21)



figure()
p=patch( 'vertices', EPI1.xyz , 'faces', EPI1.tri  ,'facecolor','interp','cdata',upstroke_epicardium_lge', 'edgecolor','none' ); 
[h,v] = IsoLine({EPI1.tri, EPI1.xyz},upstroke_epicardium_lge,20,'k');

axis off
axis equal; view(3); headlight
c=colorbar;
ylabel(c,'Upstroke Times [ms]','Fontsize',21);
set(gca,'Fontsize',21)
caxis([0 65])


csvwrite('D:\ARVC meshing automatic\patients\patient09\breakthroughs.csv',[coords(2:end,1) coords(2:end,2) coords(2:end,3)])


N=size(mesh.xyz,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %fileinfo = hdf5info('results.h5');
% t=3;
% close all
% figure()
% camlight
% % hold on
% % %htext =text(0,-0.6,sprintf('time=%f',t));
% % 
%  for t=1:7
t=10;
 output = h5read('D:\ARVC meshing automatic\patients\patient09\results\TestControlECGPat09\results.h5','/Data',[1,1,t] ,[1,N,1]);
 t=35;
  output1 = h5read('D:\ARVC meshing automatic\patients\patient09\results\TestControlECGPat09\results.h5','/Data',[1,1,t] ,[1,N,1]);

% % %delete(htext)    
% % %patch(mesh.xyz(:,1),mesh.xyz(:,2),mesh.xyz(:,3),upstroke_heart')
%  patch('Faces',mesh.tri,'Vertices',mesh.xyz,'FaceVertexCData',transpose(upstroke_heart2),'FaceColor','interp','EdgeColor','none');
% % patch('Faces',mesh.tri,'Vertices',mesh.xyz,'FaceVertexCData',transpose(output(permutation')),'FaceColor','interp','EdgeColor','none');
% % %patch('Faces',mesh.tri,'Vertices',mesh.xyz,'FaceVertexCData',distances_to_min_index,'FaceColor','interp','EdgeColor','none');
% % %patch('Faces',mesh.tri,'Vertices',mesh.xyz,'FaceVertexCData',CV','FaceColor','interp','EdgeColor','none');
% % %caxis([0 100])
% % title('conduction velocity [ms]')
% % set(gca,'FontSize',21)
% % %htext =text(0,-0.6,sprintf('time=%f',t));
% % %patch('Faces',mesh.face,'Vertices',mesh.xyz,'FaceVertexCData',transpose(upstroke(P)),'FaceColor','interp','EdgeColor','none');
% % %scatter3(mesh.xyz(ids_unstable,1),mesh.xyz(ids_unstable,2),mesh.xyz(ids_unstable,3))
% % colorbar
% % headlight
% % %caxis([-95 -80])
% % axis off
% % end
% hist(output(P'))


mesh.xyzUpstroke = transpose(upstroke_heart2)-800;
%compare results with control values for this patient
control_res = 'D:\ARVC meshing automatic\patients\patient09\results\TestControlECGPat09\';
hinfo = h5info(strcat(control_res,'results.h5'));
upstroke_control = h5read(strcat(control_res,'results.h5'),strcat('/',hinfo.Datasets(3).Name));
upstroke_control = upstroke_control(1,permutation,1);
mesh.xyzUpstrokeLGE1 = transpose(upstroke_control);

%mesh.xyzUpstrokeDifference = transpose(upstroke_heart2- upstroke_control);
%mesh.xyzControl=upstroke_control';
mesh.xyzV10=output(permutation)';
mesh.xyzV35=output1(permutation)';

upstroke_diff= transpose(upstroke_heart2- upstroke_control);
mesh.celltype =10;
write_VTK(mesh,'D:/ARVC meshing automatic/patients/patient09/results_roots_base.vtk','binary');


