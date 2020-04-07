close all
addpath(genpath('C:/users/petnov/Dropbox/shared - copy/'));
addpath(genpath('C:/Users/petnov/Dropbox/qrs/'));
addpath(genpath('C:/Users/petnov/Dropbox/mesh_scripts/'));
addpath(genpath('C:/users/petnov/Dropbox/shared - copy/IO/'));
addpath(genpath('C:\Users\petnov\Dropbox\sim3d scripts\'));
addpath(genpath('C:\Users\petnov\Dropbox\shared - copy\IO\'));



enableVTK;


 ecg=load('D:\ARVC meshing automatic\patients\patient06\results\Eikonal\exp13_ECGs\17.0x_0f_1.0fiber_roots_4_no_rv_pur_0_NEW_pred_ATMap.csvpseudo_ECG.mat');
    ecg=ecg.vars;
    Icon=ecg(:,1);
    IIcon=ecg(:,2);
    IIIcon=ecg(:,3);
    aVRcon=ecg(:,4);
    aVLcon=ecg(:,5);
    aVFcon=ecg(:,6);
    V1wcon=ecg(:,7);
    V2wcon=ecg(:,8);
    V3wcon=ecg(:,9);
    V4wcon=ecg(:,10);
    V5wcon=ecg(:,11);
    V6wcon=ecg(:,12);
    
%%
electrodes = load('D:\ARVC meshing automatic\patients\patient06\mpp\ECG_ELECTRODES.mat');

E = [ electrodes.ELECTRODES.LA ; electrodes.ELECTRODES.RA ; electrodes.ELECTRODES.LL ; electrodes.ELECTRODES.RL ; electrodes.ELECTRODES.V1 ; electrodes.ELECTRODES.V2 ; electrodes.ELECTRODES.V3 ; electrodes.ELECTRODES.V4 ; electrodes.ELECTRODES.V5 ; electrodes.ELECTRODES.V6 ];
E=E/10; %chaste is in cm
%plotMESH( MakeMesh( H , H.face ) , 'r' ); hplotMESH( T , 'nf' );
%hplot3d( E , '*m' ,'LineWidth',4,'MarkerSize',10 );
%axis(objbounds)

[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e,Ie_f,IIe_f,IIIe_f,aVRe_f,aVLe_f,aVFe_f,V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f,ts_ecgi]=ecgi_ECG_pat09('D:\Data_ARVC_raw_electrodes_vest\Data_ARVC\ARVC9\ARVC_Baseline.mat','D:\ChrisReconstructions\ARVC9\ARVC9_Baseline.mat','D:\ARVC meshing automatic\patients\patient09\mpp\ECG_ELECTRODES.mat',9,I,II,III,aVR,aVL,aVF,V1w,V2w,V3w,V4w,V5w,V6w);

nolge=horzcat(V1wcon,V2wcon,V3wcon,V4wcon,V5wcon,V6wcon);

pop_path='D:\ARVC meshing automatic\patients\patient06\';
plt_acts=0;
keep_models={};
patient_nr='06';
sim_name='Eikonal';
study='exp13';% permutation = load('permutation.txt'); %the nodes are renumbered by chaste!!
% P=permutation+1; %matlab starts ordering from 1
mesh = read_VTK(strcat(pop_path,'HEART.vtk'));
mesh_large=read_CHASTE(strcat(pop_path,'chaste',patient_nr,'\HEART'));
nr_frames = 1000;
v = VideoWriter(strcat(pop_path,'video1.avi','Motion JPEG AVI'));
open(v)
F=[];

%cd D:\ARVC' meshing automatic'\patients\patient05\results\TestHeart05Full
fileinfo = hdf5info('results.h5');
acts=dlmread(strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\results\',sim_name,'\',study,'\17.0x_0f_1.0fiber_roots_4_no_rv_pur_0_NEW_pred_ATMap','.csv'));
Voltage=ones(size(acts,1),1)*-90;
for i=1:max(acts)
    h = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on
    subplot(2,2,1)
    id=find(acts==i);
    Voltage(id)=30;
    output_t=Voltage;
    output_epi=output_t(1:size(mesh_large.face,1));
    %set(h, 'Visible', 'off');
    %color=(output_epi+max(-output_epi))/(max(output_epi+max(-output_epi))+0.001).*ones(size(mesh.tri,1),1);
    patch('Faces',mesh_large.face,'Vertices',mesh_large.xyz,'FaceVertexCData',output_t,'Facecolor','flat','EdgeColor','none');
    view([-148,16])
    %view([35 -91])
     % view([33 90])
    camroll(-180)
    axis off
    camlight
    axis equal
        subplot(2,2,2)
  patch('Faces',mesh_large.face,'Vertices',mesh_large.xyz,'FaceVertexCData',output_t,'Facecolor','flat','EdgeColor','none');
    hold on
    scatter3(electrodes.ELECTRODES.V1(1)/10,electrodes.ELECTRODES.V1(2)/10,electrodes.ELECTRODES.V1(3)/10,'m','filled')
        text(-3+electrodes.ELECTRODES.V1(1)/10,-3+electrodes.ELECTRODES.V1(2)/10,-3+electrodes.ELECTRODES.V1(3)/10,'Lead V1','FontSize',21)
           %text(-3+electrodes.ELECTRODES.V3(1)/10,-3+electrodes.ELECTRODES.V3(2)/10,-3+electrodes.ELECTRODES.V3(3)/10,'Lead V3','FontSize',21)

        zoom(1.4)

    %view([-13,66])
        view([33 90])

    axis off
    camlight
    axis equal

    subplot(2,2,[3 4])
    scatter(1:size(V1wcon,1),V1wcon,'b')
    hold on
    scatter(i,V1wcon(i),'r','filled')
    F=[F,getframe(h)];
    close(h);
    i
    if rem(i,10)==0
        writeVideo(v,F);
        clear F;
        F=[];
    end

end


%fig = figure;
%movie(fig,F(1:110),1)
close(v)
