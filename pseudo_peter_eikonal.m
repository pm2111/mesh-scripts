addpath(genpath('C:/Users/petnov/Dropbox/shared - copy'));
addpath(genpath('C:/Users/petnov/Dropbox/mesh scripts/'));

enableVTK;
sim_name='Eikonal_coarse';
pat_nr='06';

FBfine =transform(read_CHASTE(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\chaste',pat_nr,'_coarse\HEART')),'s',1/1);
H.face      = MeshBoundary( FBfine.tri );
H.triCENTER = meshFacesCenter( FBfine );
H.triVOL    = meshQuality( FBfine , 'volume' );
H.xyzIDS    = vtkClosestPoint( Mesh( FBfine ) , FBfine.xyz );
H.Gop       = meshGradient( FBfine );
H.tri =FBfine.tri;
H.xyz=FBfine.xyz;
N=size(FBfine.xyz,1);


%%
electrodes = load(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\mpp\ECG_ELECTRODES.mat'));
T=load(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\mpp\fittedVEST.mat'));

E = [ electrodes.ELECTRODES.LA ; electrodes.ELECTRODES.RA ; electrodes.ELECTRODES.LL ; electrodes.ELECTRODES.RL ; electrodes.ELECTRODES.V1 ; electrodes.ELECTRODES.V2 ; electrodes.ELECTRODES.V3 ; electrodes.ELECTRODES.V4 ; electrodes.ELECTRODES.V5 ; electrodes.ELECTRODES.V6 ];
E=E/10;
%body surface pot leads
T.VEST.xyz=T.VEST.xyz/10;
cols={'r','g','b','m','c','k'};
figure()
plotMESH( H , 'ne' ); 
hplotMESH( T.VEST , 'nf' );
axis equal
for i=5:10
    hold on
    scatter3( E(i,1), E(i,2), E(i,3),65,cols{i-4} ,'filled');
    axis(objbounds)
end
legend({'V1','V2','V3','V4','V5','V6'})
%%save the voltages into 1ms bins

%fileinfo = hdf5info('D:/ARVC meshing automatic/patients/patient18/results/TestHeart18ControlECG/results.h5');
%permutation = load(strcat('D:/ARVC meshing automatic/patients/patient',pat_nr,'/results/',sim_name,'/permutation.txt'));
%permutation = permutation+1;
voltage= dlmread('C:\Users\petnov\Dropbox\cell models\voltage_data_epi.csv',',',1,0);

    sim_type='eikonal06_scar\Purkinje_Speed_17.0_fib_1_scar_1.0\Purkinje_Speed_17.0_fib_1_scar_1.0_NEW_pred_ATMap';
    acts=dlmread(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',sim_name,'\',sim_type,'.csv'));
%save_dir=strcat('exp2_ECGs\17.0x_',num2str(i),'f_septal_',num2str(i+2),'f_1.0fiber_roots_0_no_rv_pur_0');
%acts=dlmread(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',sim_name,'\',sim_type,'.csv'));

%acts_slow=dlmread(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',sim_name,'\eikonal06_fine_NEW_pred_ATMap_250_2x_1f_NormalRoots.csv'));
%   figure()
%     headlight
%     hold on
%     patch('vertices',FBfine.xyz,'faces',FBfine.face,'edgecolor','none','FaceVertexCData',acts,'facecolor','interp')
%     colorbar
%     caxis([0 50])
%     axis off

ecg_full=[zeros(1,12)];
voltages_heart=voltage(9:end);
batch_size=5;
t_end=max(acts);
nr_batches=int16(t_end/batch_size)+1;
indices_activated=[];
time_since_activated=ones(size(acts));
for i=1:400
%           activated=find(acts==id);
%           activated1=find(acts==id+1);
%           activated2=find(acts==id+2);
%           activated3=find(acts==id+3);
%           activated4=find(acts==id+4);
%           time_since_act=id*ones(size(acts,1),1)-acts;
%           need_updating_AP=find(time_since_act>0);
%           
%           voltages_heart(activated)=apd_epi(1);
%           voltages_heart(activated1)= -89.0+ramp(1);
%           voltages_heart(activated2)= -89.0+ramp(2);
%           voltages_heart(activated3)= -89.0+ramp(3);
%           voltages_heart(activated4)= -89.0+ramp(4);
            
           indices_activated=[indices_activated;find(acts==i)];
           time_since_activated(indices_activated)=time_since_activated(indices_activated)+1;
           
          voltages_step=voltages_heart(time_since_activated);
           
           
          disp(i);
 
          %gradV(t}=reshape( H.Gop * getv( voltages_heart' , H.xyzIDS ) ,[],3);
          H.triG(:,:,1) = reshape( H.Gop * getv( voltages_step , H.xyzIDS ) ,[],3);
    
    [~,I,II,III,aVR,aVL,aVF,V1,V2,V3,V4,V5,V6] = computePECG(  H , E );
    int=[I(end),II(end),III(end),aVR(end),aVL(end),aVF(end),V1(end),V2(end),V3(end),V4(end),V5(end),V6(end)];
    ecg_full = vertcat(ecg_full,int);
    H=rmfield(H,'triG');
end
save(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\','\','ecgs\thin_RV_slow_prukinje_repol','.mat'),'ecg_full');

ecg=ecg_full/10;
  Ie=ecg(:,1); 
        IIe=ecg(:,2);
        IIIe=ecg(:,3);
        aVRe=ecg(:,4);
        aVLe=ecg(:,5);
        aVFe=ecg(:,6);
        V1e=ecg(:,7);
        V2e=ecg(:,8);
        V3e=ecg(:,9);
        V4e=ecg(:,10);
        V5e=ecg(:,11);
        V6e=ecg(:,12);
      % [Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e]=ecgi_ECG_sim_plot_mat(ecg1,V3e_f,0.1);
    ecgs{i-2}=[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e];
    full_ecg={};
        full_ecg={Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e};
        figure()
        plot_ecgs_2x6_t_wave(full_ecg','1',0,1,1:size(full_ecg{1},1));
        
        
        
        
        
            hold on
 precordial=[V1e,V2e,V3e,V4e,V5e,V6e];
    for j=1:6
       % [QRSs(j,i-2)]=qrs_estimator(precordial(:,j));
        [~,ratio(j,i-2),~,~,~,~,~,R_peak(j,i-2)]=qrs_estimator(precordial(:,j));
        [ ~,~,duration(j,i-2),swave(j,i-2),QRSs(j,i-2),S_start(j,i-2)]=s_wave_area(precordial(:,j));
    end
   
%    
%     
%     else
%         dummy_name=model_names{k-1};
%         dummy_name(7)='2';
%         A=strfind({names.name},dummy_name);
%         idx=find(not(cellfun('isempty',A)));




%save(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',sim_name,'\',save_dir,'.mat'),'ecg_full');

roots=dlmread(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\start_pos.txt'),',');
figure()
    patch('Faces',H.face,'Vertices',H.xyz,'FaceVertexCData',acts,'Facecolor','interp','EdgeColor','none');
headlight
hold on
scatter3(roots(:,1),roots(:,2),roots(:,3),135,'filled')
colorbar




[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e,Ie_f,IIe_f,IIIe_f,aVRe_f,aVLe_f,aVFe_f,V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f,ts_ecgi]=ecgi_ECG('D:\Data_ARVC_raw_electrodes_vest\Data_ARVC\ARVC18\ARVC18_Baseline.mat','D:\ChrisReconstructions\ARVC18\ARVC18_Baseline.mat','D:\ARVC meshing automatic\patients\patient18\mpp\ECG_ELECTRODES.mat',I_real,II_real,III_real,aVR_real,aVL_real,aVF_real,V1_real,V2_real,V3_real,V4_real,V5_real,V6_real);
ts=1:200;
offset=40;%idx
ts_ecgi=ts_ecgi-ts_ecgi(40);
figure('units','normalized','outerposition',[ 0 0 1 1])
grid on
i=1;
subplot('Position',[(.07+mod(i-1,3))/3 1.1-(ceil(i/3))/3 1/4 1/5])
%plot(ts(1:size(I_real,1)),I_real,'LineWidth',3)
hold all
%plot(ts_,'--r','LineWidth',3)
plot(ts_ecgi(offset:end),Ie(offset:end),'g','LineWidth',3)
plot(ts_ecgi(offset:end),Ie_f(offset:end),'k','LineWidth',3)
xlabel('ms')
ylabel('mV')
vec_pos=get(get(gca,'XLabel'),'Position');
set(get(gca,'XLabel'),'Position',vec_pos+[40 0 0])
title('Lead I')
grid on
i=2;
subplot('Position',[(.07+mod(i-1,3))/3 1.1-(ceil(i/3))/3 1/4 1/5])
%plot(ts(1:size(II_real,1)),II_real,'LineWidth',3)
hold on
%plot(III_real,'--r','LineWidth',3)
plot(ts_ecgi(offset:end),IIe(offset:end),'g','LineWidth',3)
plot(ts_ecgi(offset:end),IIe_f(offset:end),'k','LineWidth',3)
xlabel('ms')
ylabel('mV')
vec_pos=get(get(gca,'XLabel'),'Position');
set(get(gca,'XLabel'),'Position',vec_pos+[40 0 0])
title('Lead II')
grid on

i=3;
subplot('Position',[(.07+mod(i-1,3))/3 1.1-(ceil(i/3))/3 1/4 1/5])
%plot(ts,III_real,'LineWidth',3)
hold on
%plot(LA{1,1},IIIc,'--r','LineWidth',3)
plot(ts_ecgi(offset:end),IIIe(offset:end),'g','LineWidth',3)
plot(ts_ecgi(offset:end),IIIe_f(offset:end),'k','LineWidth',3)
xlabel('ms')
ylabel('mV')
vec_pos=get(get(gca,'XLabel'),'Position');
set(get(gca,'XLabel'),'Position',vec_pos+[40 0 0])
title('Lead III')
grid on

i=4;
subplot('Position',[(.07+mod(i-1,3))/3 1.19-(ceil(i/3))/3 1/4 1/5])
%plot(ts,aVR_real,'LineWidth',3,'LineWidth',3)
hold on
%plot(LA{1,1},aVRc,'--r','LineWidth',3)
plot(ts_ecgi(offset:end),aVRe(offset:end),'g','LineWidth',3)
plot(ts_ecgi(offset:end),aVRe_f(offset:end),'k','LineWidth',3)
xlabel('ms')
ylabel('mV')
vec_pos=get(get(gca,'XLabel'),'Position');

set(get(gca,'XLabel'),'Position',vec_pos+[40 0 0])

vec_pos=get(get(gca,'XLabel'),'Position');
title('Lead aVR')
title('Lead aVR')
grid on

i=5;
subplot('Position',[(.07+mod(i-1,3))/3 1.19-(ceil(i/3))/3 1/4 1/5])
%plot(ts,aVL_real,'LineWidth',3)
hold on
%plot(LA{1,1},aVLc,'--r','LineWidth',3)
plot(ts_ecgi(offset:end),aVLe(offset:end),'g','LineWidth',3)
plot(ts_ecgi(offset:end),aVLe_f(offset:end),'k','LineWidth',3)
xlabel('ms')
ylabel('mV')
vec_pos=get(get(gca,'XLabel'),'Position');

set(get(gca,'XLabel'),'Position',vec_pos+[40 10 0])
title('Lead aVL')
title('Lead aVL')
grid on

i=6;
subplot('Position',[(.07+mod(i-1,3))/3 1.19-(ceil(i/3))/3 1/4 1/5])
%plot(ts,aVF_real,'LineWidth',3)
hold on
%plot(LA{1,1},aVFc,'--r','LineWidth',3)
plot(ts_ecgi(offset:end),aVFe(offset:end),'g','LineWidth',3)
plot(ts_ecgi(offset:end),aVFe_f(offset:end),'k','LineWidth',3)


xlabel('ms')
ylabel('mV')
vec_pos=get(get(gca,'XLabel'),'Position');

set(get(gca,'XLabel'),'Position',vec_pos+[40 0 0])

title('Lead aVF')
grid on

i=7;
subplot('Position',[(.07+mod(i-1,3))/3 1.27-(ceil(i/3))/3 1/4 1/5])
%plot(ts,V1_real,'LineWidth',3)
hold on
%plot(LA{1,1},V1wc,'--r','LineWidth',3)
plot(ts_ecgi(offset:end),V1e(offset:end),'g','LineWidth',3)
plot(ts_ecgi(offset:end),V1e_f(offset:end),'k','LineWidth',3)

xlabel('ms')
ylabel('mV')
vec_pos=get(get(gca,'XLabel'),'Position');

set(get(gca,'XLabel'),'Position',vec_pos+[40 0 0])
title('Lead V1')
grid on

i=8;
subplot('Position',[(.07+mod(i-1,3))/3 1.27-(ceil(i/3))/3 1/4 1/5])
%plot(ts,V2_real,'LineWidth',3)
hold on
%plot(LA{1,1},V2wc,'--r','LineWidth',3)
plot(ts_ecgi(offset:end),V2e(offset:end),'g','LineWidth',3)
plot(ts_ecgi(offset:end),V2e_f(offset:end),'k','LineWidth',3)

title('Lead V2')
grid on
xlabel('ms')
ylabel('mV')
vec_pos=get(get(gca,'XLabel'),'Position');

set(get(gca,'XLabel'),'Position',vec_pos+[40 0 0])
i=9;
subplot('Position',[(.07+mod(i-1,3))/3 1.27-(ceil(i/3))/3 1/4 1/5])
%plot(ts,V3_real,'LineWidth',3)
hold on
%plot(LA{1,1},V3wc,'--r','LineWidth',3)
plot(ts_ecgi(offset:end),V3e(offset:end),'g','LineWidth',3)
plot(ts_ecgi(offset:end),V3e_f(offset:end),'k','LineWidth',3)
xlabel('ms')
ylabel('mV')
vec_pos=get(get(gca,'XLabel'),'Position');

set(get(gca,'XLabel'),'Position',vec_pos+[40 0 0])
title('Lead V3')
grid on

i=10;
subplot('Position',[(.07+mod(i-1,3))/3 1.36-(ceil(i/3))/3 1/4 1/5])
%plot(ts,V4_real,'LineWidth',3)
hold on
%plot(LA{1,1},V4wc,'--r','LineWidth',3)
plot(ts_ecgi(offset:end),V4e(offset:end),'g','LineWidth',3)
plot(ts_ecgi(offset:end),V4e_f(offset:end),'k','LineWidth',3)
xlabel('ms')
ylabel('mV')
vec_pos=get(get(gca,'XLabel'),'Position');

set(get(gca,'XLabel'),'Position',vec_pos+[40 0 0])
title('Lead V4')
grid on

i=11;
subplot('Position',[(.07+mod(i-1,3))/3 1.36-(ceil(i/3))/3 1/4 1/5])
%plot(ts,V5_real,'LineWidth',3)
hold on
%plot(LA{1,1},V5wc,'--r','LineWidth',3)
plot(ts_ecgi(offset:end),V5e(offset:end),'g','LineWidth',3)
plot(ts_ecgi(offset:end),V5e_f(offset:end),'k','LineWidth',3)
xlabel('ms')
ylabel('mV')
vec_pos=get(get(gca,'XLabel'),'Position');

set(get(gca,'XLabel'),'Position',vec_pos+[40 0 0])
title('Lead V5')
grid on

i=12;
subplot('Position',[(.07+mod(i-1,3))/3 1.36-(ceil(i/3))/3 1/4 1/5])
%plot(ts,V6_real,'LineWidth',3)
hold on
%plot(LA{1,1},V6wc,'--r','LineWidth',3)
plot(ts_ecgi(offset:end),V6e(offset:end),'g','LineWidth',3)
plot(ts_ecgi(offset:end),V6e_f(offset:end),'k','LineWidth',3)
xlabel('ms')
ylabel('mV')
vec_pos=get(get(gca,'XLabel'),'Position');

set(get(gca,'XLabel'),'Position',vec_pos+[40 0 0])
title('Lead V6')
grid on

legend('body surface pots', 'filtered body surface pots');
%%
mesh.xyzVoltages=voltages_heart;
mesh.celltype=10;
write_VTK(mesh,'D:/ARVC meshing automatic/patients/patient18/results/TestHeart18ControlECG/voltages200.vtk','binary');

figure;
subplot(1,8,1);  plot( ecg_full(1,:)  )
subplot(1,8,2); hplot( ecg_full(2,:) )
subplot(1,8,3); hplot( ecg_full(3,4) )
subplot(1,8,4); hplot( ecg_full(4,:) )
subplot(1,8,5); hplot( ecg_full(5,:) )
subplot(1,8,6); hplot( ecg_full(6,:) )
subplot(1,8,7); hplot( ecg_full(7,:) )
subplot(1,8,8); hplot( ecg_full(8,:) )
set( gaa(gcf) , 'YLim',[-90 50] )


%%
