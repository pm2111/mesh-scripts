addpath(genpath('C:/Users/petnov/Dropbox/shared - copy'));
enableVTK;
sim_name='Eikonal/exp1/17.0x_0f_1.0fiber_roots_0_no_rv_pur_0/';
pat_nr='09';

FBfine =transform(read_CHASTE(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\chaste',pat_nr,'\HEART')),'s',1/1);
H.face      = MeshBoundary( FBfine.tri );
H.triCENTER = meshFacesCenter( FBfine );
H.triVOL    = meshQuality( FBfine , 'volume' );
H.xyzIDS    = vtkClosestPoint( Mesh( FBfine ) , FBfine.xyz );
H.Gop       = meshGradient( FBfine );
H.tri =FBfine.tri;
H.xyz=FBfine.xyz;

mesh=struct();
mesh.xyz=H.xyz;
mesh.tri=H.tri;

%%
electrodes = load(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\mpp\ECG_ELECTRODES.mat'));


E = [ electrodes.ELECTRODES.LA ; electrodes.ELECTRODES.RA ; electrodes.ELECTRODES.LL ; electrodes.ELECTRODES.RL ; electrodes.ELECTRODES.V1 ; electrodes.ELECTRODES.V2 ; electrodes.ELECTRODES.V3 ; electrodes.ELECTRODES.V4 ; electrodes.ELECTRODES.V5 ; electrodes.ELECTRODES.V6 ];
E=E/10; %chaste is in cm
%plotMESH( MakeMesh( H , H.face ) , 'r' ); hplotMESH( T , 'nf' );
%hplot3d( E , '*m' ,'LineWidth',4,'MarkerSize',10 );
%axis(objbounds)

%%save the voltages into 1ms bins

%fileinfo = hdf5info('D:/ARVC meshing automatic/patients/patient18/results/TestHeart18ControlECG/results.h5');
permutation = load(strcat('D:/ARVC meshing automatic/patients/patient',pat_nr,'/results/',sim_name,'/permutation.txt'));
permutation = permutation+1;
N=size(FBfine.xyz,1);

ecg_full=[];
    fileinfo = hdf5info(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',sim_name,'\results.h5'));
for n=1:19

    for i=1:10
           % voltages = h5read('TestHeart05FullBeatandActivation\results.h5',fileinfo.GroupHierarchy.Datasets(1).Name,[1,1,t] ,[1,N,1]);
        voltages = h5read(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',sim_name,'\results.h5'),'/Data',[1,1,+i+800+(n-1)*10] ,[1,N,1]);

            voltages_heart = voltages(permutation); %select only the heart activation times!! 
          disp(i);

          %gradV(t}=reshape( H.Gop * getv( voltages_heart' , H.xyzIDS ) ,[],3);
          H.triG(:,:,i) = reshape( H.Gop * getv( voltages_heart' , H.xyzIDS ) ,[],3);
    end
    [~,I,II,III,aVR,aVL,aVF,V1,V2,V3,V4,V5,V6] = computePECG(  H , E );
    ecg_full = [ecg_full,I,II,III,aVR,aVL,aVF,V1,V2,V3,V4,V5,V6];
    H=rmfield(H,'triG');
end
I_real =[];
II_real= [];
III_real= [];
aVR_real= [];
aVL_real= [];
aVF_real=[];
V1_real= [];
V2_real= [];
V3_real=[];
V4_real=[];
V5_real=[];
V6_real= [];
for i=1:19
    I_real= vertcat(I_real,ecg_full(:,1+(i-1)*12));
    II_real= vertcat(II_real,ecg_full(:,2+(i-1)*12));
    III_real= vertcat(III_real,ecg_full(:,3+(i-1)*12));
    aVR_real= vertcat(aVR_real,ecg_full(:,4+(i-1)*12));
    aVL_real= vertcat(aVL_real,ecg_full(:,5+(i-1)*12));
    aVF_real= vertcat(aVF_real,ecg_full(:,6+(i-1)*12));
    V1_real= vertcat(V1_real,ecg_full(:,7+(i-1)*12));
    V2_real= vertcat(V2_real,ecg_full(:,8+(i-1)*12));
    V3_real= vertcat(V3_real,ecg_full(:,9+(i-1)*12));
    V4_real= vertcat(V4_real,ecg_full(:,10+(i-1)*12));
    V5_real= vertcat(V5_real,ecg_full(:,11+(i-1)*12));
    V6_real= vertcat(V6_real,ecg_full(:,12+(i-1)*12));
end


save(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',sim_name,'\ecg.mat'),'I_real','II_real','III_real','aVR_real','aVL_real','aVF_real','V1_real','V2_real','V3_real','V4_real','V5_real','V6_real');

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
