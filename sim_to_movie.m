close all
addpath(genpath('C:/users/petnov/Dropbox/shared - copy/'));
addpath(genpath('C:/Users/petnov/Dropbox/qrs/'));
addpath(genpath('C:/Users/petnov/Dropbox/mesh_scripts/'));
addpath(genpath('C:/users/petnov/Dropbox/shared - copy/IO/'));
addpath(genpath('C:\Users\petnov\Dropbox\sim3d scripts\'));
addpath(genpath('C:\Users\petnov\Dropbox\shared - copy\IO\'));



enableVTK;
%[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e,Ie_f,IIe_f,IIIe_f,aVRe_f,aVLe_f,aVFe_f,V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f,ts_ecgi]=ecgi_ECG06('D:\Data_ARVC_raw_electrodes_vest\Data_ARVC\ARVC6\ARVC6_Baseline.mat','D:\ChrisReconstructions\ARVC6\ARVC6_Baseline.mat','D:\ARVC meshing automatic\patients\patient06\mpp\ECG_ELECTRODES.mat',6);
%ecg_ecgi=[V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f];

  %  offset=80;
   % offset2=size(Ie_f,1)-270;
   % full_ecg_ecgi={real(Ie_f(offset:2:offset2)),real(IIe_f(offset:2:offset2)),real(IIIe_f(offset:2:offset2)),real(aVRe_f(offset:2:offset2)),real(aVLe_f(offset:2:offset2)),real(aVFe_f(offset:2:offset2)),real(V1e_f(offset:2:offset2))-real(V1e_f(offset)),real(V2e_f(offset:2:offset2))-real(V2e_f(offset)),real(V3e_f(offset:2:offset2))-real(V3e_f(offset)),real(V4e_f(offset:2:offset2))-real(V4e_f(offset)),real(V5e_f(offset:2:offset2))-real(V5e_f(offset)),real(V6e_f(offset:2:offset2))-real(V6e_f(offset))};
 %plot clinical ecgs

 ecg=load('D:\ARVC meshing automatic\patients\patient06\results\Eikonal_coarse\scar_within_scar_non_uniform_ECGs\Purkinje_Speed_17.0_fib_0_scar_grad6scar_size_6epi_intact_1_NEW_pred_ATMap.csvpseudo_ECG.mat');
    ecg=ecg.vars;
  Ie=ecg(:,1)/10; 
        IIe=ecg(:,2)/10;
        IIIe=ecg(:,3)/10;
        aVRe=ecg(:,4)/10;
        aVLe=ecg(:,5)/10;
        aVFe=ecg(:,6)/10;
        V1e=ecg(:,7)/10;
        V2e=ecg(:,8)/10;
        V3e=ecg(:,9)/10;
        V4e=ecg(:,10)/10;
        V5e=ecg(:,11)/10;
        V6e=ecg(:,12)/10;
    full_ecg={};
        full_ecg={Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e};

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
sim_name='Eikonal_coarse';
study='scar_within_scar_non_uniform_AT';% permutation = load('permutation.txt'); %the nodes are renumbered by chaste!!
% P=permutation+1; %matlab starts ordering from 1
%mesh = read_VTK(strcat(pop_path,'HEART.vtk'));
mesh_large=read_CHASTE(strcat(pop_path,'HEART'));
nr_frames = 95;
v = VideoWriter(strcat(pop_path,'activation_epsilon_pos.avi','Motion JPEG AVI'));
open(v)
F=[];

%cd D:\ARVC' meshing automatic'\patients\patient05\results\TestHeart05Full
%fileinfo = hdf5info('results.h5');
acts=load(strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\results\',sim_name,'\','scar_within_scar_non_uniform_AT\Purkinje_Speed_17.0_fib_0_scar_grad6scar_size_6epi_intact_1\Purkinje_Speed_17.0_fib_0_scar_grad6scar_size_6epi_intact_1_NEW_pred_ATMap','.csv'));
Voltage=ones(size(acts,1),1)*-90;
for i=1:max(acts)
    h = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on
    subplot(2,2,1)
    id=find(acts==i);
    Voltage(id)=30;
    output_t=Voltage;
    %set(h, 'Visible', 'off');
    %color=(output_epi+max(-output_epi))/(max(output_epi+max(-output_epi))+0.001).*ones(size(mesh.tri,1),1);
    patch('Faces',mesh_large.face,'Vertices',mesh_large.xyz,'FaceVertexCData',output_t,'Facecolor','flat','EdgeColor','none');
    view([-154,31])
    caxis([0 100])
    colormap(flipud(jet))
        headlight
    %view([35 -91])
   % camroll(-180)
    axis off
    camlight
    axis equal
        %zoom(1.2);

        subplot(2,2,2)
  patch('Faces',mesh_large.face,'Vertices',mesh_large.xyz,'FaceVertexCData',output_t,'Facecolor','flat','EdgeColor','none');
    hold on
        caxis([0 100])
    colormap(flipud(jet))
        headlight
    scatter3(electrodes.ELECTRODES.V2(1)/10,electrodes.ELECTRODES.V2(2)/10,electrodes.ELECTRODES.V2(3)/10,'m','filled')
        text(-3+electrodes.ELECTRODES.V2(1)/10,-3+electrodes.ELECTRODES.V2(2)/10,-3+electrodes.ELECTRODES.V2(3)/10,'V2','FontSize',21,'color','w')
           %text(-3+electrodes.ELECTRODES.V3(1)/10,-3+electrodes.ELECTRODES.V3(2)/10,-3+electrodes.ELECTRODES.V3(3)/10,'Lead V3','FontSize',21)


    %view([-13,66])
      view([-50 -32])

    axis off
    camlight
    axis equal
      %  zoom(1.2)

    subplot(2,2,[3 4])
    scatter(1:size(ecg,1),ecg(:,8),'b')
    hold on
    scatter(i,ecg(i,8),'r','filled')
    grid on
    xlim([0 1000])
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
