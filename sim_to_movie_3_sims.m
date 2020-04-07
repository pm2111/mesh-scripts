close all
addpath(genpath('C:/users/petnov/Dropbox/shared - copy/'));
addpath(genpath('C:/Users/petnov/Dropbox/qrs/'));
addpath(genpath('C:/Users/petnov/Dropbox/mesh_scripts/'));
addpath(genpath('C:/users/petnov/Dropbox/shared - copy/IO/'));
addpath(genpath('C:\Users\petnov\Dropbox\sim3d scripts\'));
addpath(genpath('C:\Users\petnov\Dropbox\shared - copy\IO\'));



enableVTK;

%%
electrodes = load('D:\ARVC meshing automatic\patients\patient06\mpp\ECG_ELECTRODES.mat');

E = [ electrodes.ELECTRODES.LA ; electrodes.ELECTRODES.RA ; electrodes.ELECTRODES.LL ; electrodes.ELECTRODES.RL ; electrodes.ELECTRODES.V1 ; electrodes.ELECTRODES.V2 ; electrodes.ELECTRODES.V3 ; electrodes.ELECTRODES.V4 ; electrodes.ELECTRODES.V5 ; electrodes.ELECTRODES.V6 ];
E=E/10; %chaste is in cm
%plotMESH( MakeMesh( H , H.face ) , 'r' ); hplotMESH( T , 'nf' );
%hplot3d( E , '*m' ,'LineWidth',4,'MarkerSize',10 );
%axis(objbounds)

%[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e,Ie_f,IIe_f,IIIe_f,aVRe_f,aVLe_f,aVFe_f,V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f,ts_ecgi]=ecgi_ECG_pat09('D:\Data_ARVC_raw_electrodes_vest\Data_ARVC\ARVC9\ARVC_Baseline.mat','D:\ChrisReconstructions\ARVC9\ARVC9_Baseline.mat','D:\ARVC meshing automatic\patients\patient09\mpp\ECG_ELECTRODES.mat',6);


pop_path='D:\ARVC meshing automatic\patients\patient06\chaste06_coarse\';
plt_acts=0;
keep_models={};
patient_nr='06';
folder_name='Eikonal_coarse';
study='figure_acts_lat';% permutation = load('permutation.txt'); %the nodes are renumbered by chaste!!
% P=permutation+1; %matlab starts ordering from 1
mesh = read_VTK(strcat(pop_path,'HEART.vtk'));
mesh_large=read_CHASTE(strcat(pop_path,'HEART'));
nr_frames = 95;


%cd D:\ARVC' meshing automatic'\patients\patient05\results\TestHeart05Full
%fileinfo = hdf5info('results.h5');
for i=1:3
    sim_name=dir(strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\results\',folder_name,'\',study,'\',num2str(i)));
    acts(:,i)=dlmread(strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\results\',folder_name,'\',study,'\',num2str(i),'\',sim_name(3).name));
end
names={'Slow Purkinje', 'Intermediate Purkinje', 'Control Purkinje'};


v = VideoWriter(strcat(pop_path,'activation-pruk.avi','Motion JPEG AVI'));
open(v)
F=[];
Voltage=ones(size(acts,1),3)*-90;
for i=1:100
        h = figure('units','normalized','outerposition',[0 0 1 1]);

    for j=1:3
        hold on
    subplot(1,3,j)
    id=find(acts(:,j)==i);
    Voltage(id,j)=30;
    output_t=Voltage(:,j);
    text(0.8,-8,-10,names{j},'Fontsize',40)
    if j==3
    text(0.8,-8,10,strcat('time : ',num2str(i),'ms'),'Fontsize',45)
    end
    %set(h, 'Visible', 'off');
    %color=(output_epi+max(-output_epi))/(max(output_epi+max(-output_epi))+0.001).*ones(size(mesh.tri,1),1);
    patch('Faces',mesh_large.face,'Vertices',mesh_large.xyz,'FaceVertexCData',output_t,'Facecolor','flat','EdgeColor','none');
    view([-111,24])
    caxis([-90 30])
    colormap(flipud(jet))
        headlight
    %view([35 -91])
     % view([33 90])
    camroll(-180)
    zoom(1.5)
    axis off
    camlight
    axis equal
    end
       colorbar([0.93 0.2 0.025 0.7],'Fontsize',25)
    F=[F,getframe(h)];
    i
    if rem(i,10)==0
        writeVideo(v,F);
        clear F;
        F=[];
    end
    close(h);

end


%fig = figure;
%movie(fig,F(1:110),1)
close(v)
