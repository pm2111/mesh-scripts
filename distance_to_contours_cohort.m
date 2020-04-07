addpath(genpath('C:\users\petnov\Dropbox\'));
addpath('D:/peter tools/')
patients=dir('D:\ARVC meshing automatic\patients\');
enableVTK;
distances={};
for i=1:size(patients,1)
    try
        distances{i}=read_VTK(strcat('D:\ARVC meshing automatic\patients\',patients(i+2).name,'\mpp\d2HeartSurfaces.vtk'));
    catch
            fprintf(strcat('patient', num2str(  i), 'does not have a mesh'));
    end
end
% 
% meshes = {};
% 
% for i=1:size(patients,1)
%     try
%         meshes{i}=read_VTK(strcat('D:\ARVC meshing automatic\patients\',patients(i+2).name,'\HEART_2.00.vtk'));
%     catch 
%          fprintf(strcat('patient', num2str(  i), 'does not have a mesh'));
%     end
%     try
%         meshes{i}=read_VTK(strcat('D:\ARVC meshing automatic\patients\',patients(i+2).name,'\HEART_0.20.vtk'));
%     catch 
%          fprintf(strcat('patient', num2str(  i), 'does not have a mesh'));
% 
%     end
% end
% 
% meshes_full={};
% j=1;
% for i=1:size(meshes,2)
%     if isstruct(meshes{i})==1
%         meshes_full{j}=meshes{i};
%         j=j+1;
% 
%     end
% end
% 
% meshes_full{12}=meshes_full{13};
% figure()
% for i=4:12
%     subplot(3,3,i-3)
%     plotMESH(meshes_full(i),'ne','FaceColor',[1 1 0]*0.6)
%     headlight
%     axis off
% end
% results =[];
% for i=2:size(meshes_full,2)
%     [LVVendo, RVVendo, MuscleV] = volumef(meshes_full{i});
%     results(i-1,1)=LVVendo;
%     results(i-1,2)=RVVendo;
%     results(i-1,3)=MuscleV;
% end
% 
% control_names=dir('D:\meshes\control cohort\');
% meshes_control = {};
% 
% for i=1:4
%    
%         meshes_control{i}=read_VTK(strcat('D:\meshes\control cohort\',control_names(i+2).name));
% 
% end
% results_control =[];
% for i=1:size(meshes_control,2)
%     [LVVendo, RVVendo, MuscleV] = volumef(meshes_control{i});
%     results_control(i,1)=LVVendo;
%     results_control(i,2)=RVVendo;
%     results_control(i,3)=MuscleV;
% end

d2heart={};
for i=1:size(patients,1)
    try
        d2heart{i}=read_VTK(strcat('D:\ARVC meshing automatic\patients\',patients(i+2).name,'\mpp\d2HeartSurfaces.vtk'));
        mean_dist(i) = sum(d2heart{i}.triD)/size(d2heart{i}.triD,1);
    catch
        fprintf(strcat('no distances file found for pat',num2str(i)));
    end
    
end

patients1=dir('D:/ChrisReconstructions/');
for i=1:20
    ecgi=load(strcat('D:\ChrisReconstructions\',patients1(i+2).name,'\',patients1(i+2).name,'_Baseline.mat'));
    field=fieldnames(ecgi);
    mesh=struct();
    mesh.xyz=ecgi.(field{1}).nodes;
    mesh.tri=ecgi.(field{1}).mesh;
    mesh.celltype=5;
    write_VTK(mesh,strcat('D:\meshes\','ecgi',patients1(i+2).name,'.vtk'),'binary');
end
    

figure()
scatter(1:19,mean_dist([1:11,13:20]),100,'filled')
xlabel('patient number','FontSize',21)
ylabel('mean triangle centre to closest mesh node distance [mm]','FontSize',21)
set(gca,'FontSize',21)


figure()
hist(mean_dist,25)
xlabel('mean contour to closest mesh node distance [mm]','FontSize',21)
ylabel('frequency','FontSize',21)
set(gca,'FontSize',21)
xlim([0.25 3])



figure()
scatter(ones(size(results,1),1),results(:,1),80,'filled')
hold on
scatter(ones(size(results_control,1),1)*2,results_control(:,1),80,'r','filled')
title('LV Endocardial Volume [mm^3]','Fontsize',21)
set(gca,'FontSize',21)
xlim([0.5,2.5])
legend('ARVC','Control')


figure()
scatter(ones(size(results,1),1),results(:,2),80,'filled')
hold on
scatter(ones(size(results_control,1),1)*2,results_control(:,2),80,'r','filled')

title('RV Endocardial Volume [mm^3]','Fontsize',21)
set(gca,'FontSize',21)
legend('ARVC','Control')

xlim([0.5,2.5])

figure()
scatter(ones(size(results,1),1),results(:,3),80,'filled')
hold on
scatter(ones(size(results_control,1),1)*2,results_control(:,3),80,'r','filled')

title('Total Ventricular Muscle Volume [mm^3]','Fontsize',21)
set(gca,'FontSize',21)
xlim([0.5,2.5])
legend('ARVC','Control')



