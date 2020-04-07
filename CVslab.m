addpath(genpath('/home/scratch/petnov/Dropbox/shared - copy/IO/'));
addpath(genpath('/home/scratch/petnov/Dropbox/shared - copy/'));
addpath(genpath('/home/scratch/petnov/Dropbox/mesh scripts/'))


mesh =read_CHASTE('/data/meshes/mesh0.4mm/slab_0.4mm.ele');

cd  /data/sim3d/slab/TestSlabSize0.4NormalCV0.5/
permutation = load('permutation.txt');
permutation = permutation+1;

fileinfo =h5info('results.h5');
upstroke = h5read('results.h5',strcat('/',fileinfo.Datasets(3).Name));
upstroke_heart = upstroke(1,permutation,1);
N=size(mesh.xyz,1);
t=15;
voltage =h5read('results.h5','/Data',[1,1,t] ,[1,N,1]);
figure()
patch('Faces',mesh.tri,'Vertices',mesh.xyz,'FaceVertexCData',transpose(upstroke_heart),'FaceColor','interp','EdgeColor','none');
patch('Faces',mesh.tri,'Vertices',mesh.xyz,'FaceVertexCData',transpose(voltage(permutation)),'FaceColor','interp','EdgeColor','none');
colorbar
axis equal

names =dir('/data/sim3d/slab/');
for i=1:3
    [principal_cv(i),inter_fiber_cv(i),inter_sheet_cv(i+1)]=cv_calc(strcat('/data/sim3d/slab/',names(i+2).name,'/'),mesh);

end

mesh2 =read_CHASTE('/data/meshes/mesh0.2mm/slab_0.2mm.ele');



names2 =dir('/data/sim3d/slab0.2/');

for i=1:3
    [principal_cv2(i),inter_fiber_cv2(i),inter_sheet_cv2(i)]=cv_calc(strcat('/data/sim3d/slab0.2/',names2(i+2).name,'/'),mesh2);

end
mesh3 =read_CHASTE('/data/meshes/mesh0.1mm/slab_0.1mm.ele');

names3 =dir('/data/sim3d/slab0.1/');
for i=1:3
    [principal_cv3(i),inter_fiber_cv3(i),inter_sheet_cv3(i)]=cv_calc(strcat('/data/sim3d/slab0.1/',names3(i+2).name,'/'),mesh3);

end



figure()
hold on
bar([0.25 0.5 1],[principal_cv3; principal_cv2; principal_cv]','grouped')
legend('0.1mm','0.2mm','0.4mm','Location','southeast')
set(gca,'FontSize',21)
xlabel('conductivity scaling')
ylabel('CV along fiber direction [cm/s]')
xticks([0.25 0.5 1])




figure()
hold on
bar([0.25 0.5 1],[inter_fiber_cv3; inter_fiber_cv2; inter_fiber_cv]','grouped')
legend('0.1mm','0.2mm','0.4mm','Location','southeast')
set(gca,'FontSize',21)
xlabel('conductivity scaling')
ylabel('CV between fibers [cm/s]')
xticks([0.25 0.5 1])


figure()
hold on
bar([0.25 0.5 1],[inter_sheet_cv3; inter_sheet_cv2; inter_sheet_cv(2:4)]','grouped')
legend('0.1mm','0.2mm','0.4mm','Location','southeast')
set(gca,'FontSize',21)
xlabel('conductivity scaling')
ylabel('CV intersheet [cm/s]')
xticks([0.25 0.5 1])
%CV along principal diagonal


line = [linspace(0,1,100)*2;linspace(0,1,100)*0.8;linspace(0,1,100)*0.8];

idx=knnsearch(mesh.xyz,line);
idx2=knnsearch(mesh2.xyz,line);
idx3=knnsearch(mesh3.xyz,line);

line =line';
figure()
hold on
scatter3(mesh.xyz(idx,1),mesh.xyz(idx,2),mesh.xyz(idx,3));
scatter3(line(:,1)',line(:,2)',line(:,3)','r')

dists= @(x,z,y) sqrt(x.^2+y.^2+z.^2);


at=diag_AT('/data/meshes/mesh0.4mm/slab_0.4mm.ele','/data/sim3d/slab/TestSlabSize0.4NormalCV0.5/');
at2=diag_AT('/data/meshes/mesh0.2mm/slab_0.2mm.ele','/data/sim3d/slab0.2/TestSlabSize0.2NormalCV0.5/');
at3=diag_AT('/data/meshes/mesh0.1mm/slab_0.1mm.ele','/data/sim3d/slab0.1/TestSlabSize0.1NormalCV0.5/');


figure()
hold on
plot(dists(mesh3.xyz(idx3,1),mesh3.xyz(idx3,2),mesh3.xyz(idx3,3)),at3,'Linewidth',3)
plot(dists(mesh2.xyz(idx2,1),mesh2.xyz(idx2,2),mesh2.xyz(idx2,3)),at2,'Linewidth',3)
plot(dists(mesh.xyz(idx,1),mesh.xyz(idx,2),mesh.xyz(idx,3)),at,'Linewidth',3)
legend('edge length = 0.1mm', 'edge length = 0.2mm','edge length = 0.4mm')
set(gca,'FontSize',21)
ylabel('Activation time along main diagonal [ms]')
xlabel('distance from stimulation edge [cm]')
