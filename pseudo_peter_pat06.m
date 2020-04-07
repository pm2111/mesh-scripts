addpath(genpath('C:/Users/petnov/Dropbox/'));
addpath(genpath('C:/Users/petnov/Dropbox/shared - copy\IO\'));
pat_nr='06';
sim_name='Eikonal';
enableVTK;

load('D:\ARVC meshing automatic\patients\patient06\results\TestHeart06Node2\ecg.mat');
Ic=I_real;
IIc=II_real;
IIIc=III_real;
aVRc=aVR_real;
aVLc=aVL_real;
aVFc=aVF_real;
V1wc=V1_real;
V2wc=V2_real;
V3wc=V3_real;
V4wc=V4_real;
V5wc=V5_real;
V6wc=V6_real;
load('D:\ARVC meshing automatic\patients\patient06\results\TestHeart06Node2LGE\ecg.mat');
I=I_real;
II=II_real;
III=III_real;
aVR=aVR_real;
aVL=aVL_real;
aVF=aVF_real;
V1w=V1_real;
V2w=V2_real;
V3w=V3_real;
V4w=V4_real;
V5w=V5_real;
V6w=V6_real;
%[Iei,IIei,IIIei,aVRei,aVLei,aVFei,V1wei,V2wei,V3wei,V4wei,V5wei,V6wei,LAei]=my_ecg_from_sim_eikonal('06','D:\ARVC meshing automatic\patients\patient06\results\Eikonal\',82);
ecg_eik=load('D:\ARVC meshing automatic\patients\patient06\results\Eikonal\ecg.mat');
ecg_eik=ecg_eik.ecg_full;
[Iei,IIei,IIIei,aVRei,aVLei,aVFei,V1wei,V2wei,V3wei,V4wei,V5wei,V6wei]=deal(ecg_eik(:,1) , ecg_eik(:,2),ecg_eik(:,3),ecg_eik(:,4),ecg_eik(:,5),ecg_eik(:,6),ecg_eik(:,7),ecg_eik(:,8),ecg_eik(:,9),ecg_eik(:,10),ecg_eik(:,11),ecg_eik(:,12));

%%
electrodes = load('D:\ARVC meshing automatic\patients\patient06\mpp\ECG_ELECTRODES.mat');

ecg=[ Ic IIc IIIc aVRc aVLc aVFc V1wc V2wc V3wc V4wc V5wc V6wc];
csvwrite(  strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',sim_name,'\eikonal06_fine_true_PseudoECG_150_1x.csv'),ecg);


E = [ electrodes.ELECTRODES.LA ; electrodes.ELECTRODES.RA ; electrodes.ELECTRODES.LL ; electrodes.ELECTRODES.RL ; electrodes.ELECTRODES.V1 ; electrodes.ELECTRODES.V2 ; electrodes.ELECTRODES.V3 ; electrodes.ELECTRODES.V4 ; electrodes.ELECTRODES.V5 ; electrodes.ELECTRODES.V6 ];
E=E/10; %chaste is in cm
% H=read_CHASTE('D:\ARVC meshing automatic\patients\patient06\chaste06\HEART');
% T=read_VTK('D:\ARVC meshing automatic\patients\patient06\mpp\BODY0.vtk');
% T.xyz=T.xyz/10;
% figure()
% plotMESH( MakeMesh( H , H.face ) , 'r' ); 
% hplotMESH( T , 'nf' );
% hplot3d( E , '*m' ,'LineWidth',4,'MarkerSize',5 );
% axis(objbounds)

[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e,Ie_f,IIe_f,IIIe_f,aVRe_f,aVLe_f,aVFe_f,V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f,ts_ecgi]=ecgi_ECG06('D:\Data_ARVC_raw_electrodes_vest\Data_ARVC\ARVC6\ARVC6_Baseline.mat','D:\ChrisReconstructions\ARVC6\ARVC6_Baseline.mat','D:\ARVC meshing automatic\patients\patient06\mpp\ECG_ELECTRODES.mat',6,I,II,III,aVR,aVL,aVF,V1w,V2w,V3w,V4w,V5w,V6w);


figure()
hold on
for i=1:12
    scatter(i,max(ecg_eik(:,i)));
      set(gca,'XTick',[1,2,3,4,5,6,7,8,9,10,11,12])
  set(gca,'XTickLabel',{'I','II','III','aVR','aVL','aVF','V1','V2','V3','V4','V5','V6'});
end


ampl=max(Ie_f);
Ic=Ic/max(Ic)*max(Ie_f);
IIc =IIc/max(IIc) * max(-IIe_f);
IIIc =IIIc/max(IIIc) * max(IIIe_f);
aVRc =aVRc/max(-aVRc) * max(-(aVRe_f));
aVLc = aVLc/max(aVLc) * max(aVLe_f);
aVFc = aVFc/max(aVFc) * max(-aVFe_f);
V1wc= V1wc/max(V1wc) * max(V1e_f);
V2wc = V2wc/max(V2wc) * max(V2e_f);
V3wc = V3wc/max(real(V3wc)) * max(real(V3e_f));
V4wc = V4wc/max(V4wc) * max(-V4e_f);
V5wc = V5wc/max(V5wc) * max(V5e_f);
V6wc = V6wc/max(V6wc) * max(V6e_f);


ampl=max(Ie_f);
I=I/max(I)*max(Ie_f);
II =II/max(II) * max(-IIe_f);
III =III/max(III) * max(IIIe_f);
aVR =aVR/max(-aVR) * max(-aVRe_f);
aVL = aVL/max(aVL) * max(aVLe_f);
aVF = aVF/max(aVF) * max(-aVFe_f);
V1w= V1w/max(V1w) * max(V1e_f);
V2w = V2w/max(V2w) * max(V2e_f);
V3w = V3w/max(V3w) * max(real(V3e_f));
V4w = V4w/max(V4w) * max(-V4e_f);
V5w = V5w/max(V5w) * max(V5e_f);
V6w = V6w/max(V6w) * max(V6e_f);


ampl=max(Ie_f);
Iei=Iei/max(Iei)*max(Ie_f);
IIei =IIei/max(IIei) * max(-IIe_f);
IIIei =IIIei/max(-IIIei) * max(-IIIe_f);
aVRei=aVRei/max(-aVRei) * max(-aVRe_f);
aVLei = aVLei/max(aVLei) * max(aVLe_f);
aVFei = aVFei/max(aVFei) * max(-aVFe_f);
V1wei= V1wei/max(V1wei) * max(V1e_f);
V2wei = V2wei/max(V2wei) * max(V2e_f);
V3wei = V3wei/max(V3wei) * max((V3e_f));
V4wei = V4wei/max(V4wei) * max(-V4e_f);
V5wei = V5wei/max(V5wei) * max(V5e_f);
V6wei = V6wei/max(V6wei) * max(V6e_f);


lge=1;
eikonal=1;

ts=1:208;
offset=1;%idx
ts_ecgi=ts_ecgi-ts_ecgi(offset);
figure()
grid on
i=1;
subplot('Position',[(.15+mod(i-1,3))/3 1.125-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),Ic,'LineWidth',3)
hold all
%plot(ts_,'--r','LineWidth',3)
if lge==1
    plot(ts(1:size(I,1)),I,'g','LineWidth',3)
end
if eikonal==1
    plot(ts(1:size(Iei,1)),Iei,'m','LineWidth',3)
end
plot(ts_ecgi(offset:end),Ie_f(offset:end),'k','LineWidth',3)
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
text(180,15,'Lead I','Fontsize',16)
grid on
set(gca,'Fontsize',16)

i=2;
subplot('Position',[(0.05+mod(i-1,3))/3 1.125-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),IIc,'LineWidth',3)
hold on
%plot(III_real,'--r','LineWidth',3)
if lge==1
    plot(ts(1:size(I,1)),II,'g','LineWidth',3)
end
if eikonal==1
    plot(ts(1:size(Iei,1)),IIei,'m','LineWidth',3)
end
plot(ts_ecgi(offset:end),IIe_f(offset:end),'k','LineWidth',3)
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
text(180,5,'Lead II','Fontsize',16)
grid on
set(gca,'Fontsize',16)

i=3;
subplot('Position',[(.05+mod(i-1,3))/3 1.125-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),IIIc,'LineWidth',3)
hold on
%plot(LA{1,1},IIIc,'--r','LineWidth',3)
if lge==1
    plot(ts(1:size(I,1)),III,'g','LineWidth',3)
end
if eikonal==1
    plot(ts(1:size(Iei,1)),IIIei,'m','LineWidth',3)
end
plot(ts_ecgi(offset:end),IIIe_f(offset:end),'k','LineWidth',3)
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)
text(180,30,'Lead III','Fontsize',16)
i=4;
subplot('Position',[(.15+mod(i-1,3))/3 1.195-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),aVRc,'LineWidth',3,'LineWidth',3)
hold on
%plot(LA{1,1},aVRc,'--r','LineWidth',3)
if lge==1
    plot(ts(1:size(I,1)),aVR,'g','LineWidth',3)
end
if eikonal==1
    plot(ts(1:size(Iei,1)),aVRei,'m','LineWidth',3)
end
plot(ts_ecgi(offset:end),aVRe_f(offset:end),'k','LineWidth',3)
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
text(180,7,'Lead aVR','Fontsize',16)
grid on
set(gca,'Fontsize',16)
i=5;
subplot('Position',[(.05+mod(i-1,3))/3 1.195-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),aVLc,'LineWidth',3)
hold on
%plot(LA{1,1},aVLc,'--r','LineWidth',3)
if lge==1
    plot(ts(1:size(I,1)),aVL,'g','LineWidth',3)
end
if eikonal==1
    plot(ts(1:size(Iei,1)),aVLei,'m','LineWidth',3)
end
plot(ts_ecgi(offset:end),aVLe_f(offset:end),'k','LineWidth',3)
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
text(180,25,'Lead aVL','Fontsize',16)


grid on
set(gca,'Fontsize',16)
i=6;
subplot('Position',[(.05+mod(i-1,3))/3 1.195-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),aVFc,'LineWidth',3)
hold on
%plot(LA{1,1},aVFc,'--r','LineWidth',3)
if lge==1

    plot(ts(1:size(I,1)),aVF,'g','LineWidth',3)
end
if eikonal==1
    plot(ts(1:size(Iei,1)),aVFei,'m','LineWidth',3)
end
plot(ts_ecgi(offset:end),aVFe_f(offset:end),'k','LineWidth',3)

xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)
text(180,15,'Lead aVF','Fontsize',16)


grid on

i=7;
subplot('Position',[(.15+mod(i-1,3))/3 1.279-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),V1wc,'LineWidth',3)
hold on
%plot(LA{1,1},V1wc,'--r','LineWidth',3)
if lge==1

    plot(ts(1:size(I,1)),V1w,'g','LineWidth',3)
end
if eikonal==1
    plot(ts(1:size(Iei,1)),V1wei,'m','LineWidth',3)
end
plot(ts_ecgi(offset:end),V1e_f(offset:end),'k','LineWidth',3)

xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)
text(180,30,'Lead V1','Fontsize',16)


grid on

i=8;
subplot('Position',[(.05+mod(i-1,3))/3 1.279-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),V2wc,'LineWidth',3)
hold on
%plot(LA{1,1},V2wc,'--r','LineWidth',3)
if lge==1

    plot(ts(1:size(I,1)),V2w,'g','LineWidth',3)
end
if eikonal==1
    plot(ts(1:size(Iei,1)),V2wei,'m','LineWidth',3)
end
 plot(ts_ecgi(offset:end),V2e_f(offset:end),'k','LineWidth',3)
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)

text(180,10,'Lead V2','Fontsize',16)


grid on

i=9;
subplot('Position',[(0.05+mod(i-1,3))/3 1.279-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),V3wc,'LineWidth',3)
hold on
%plot(LA{1,1},V3wc,'--r','LineWidth',3)
if lge==1
    plot(ts(1:size(I,1)),V3w,'g','LineWidth',3)
end
if eikonal==1
    plot(ts(1:size(Iei,1)),V3wei,'m','LineWidth',3)
end
plot(ts_ecgi(offset:end),V3e_f(offset:end),'k','LineWidth',3)
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)

text(180,10,'Lead V3','Fontsize',16)


grid on

i=10;
subplot('Position',[(.15+mod(i-1,3))/3 1.36-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),V4wc,'LineWidth',3)
hold on
%plot(LA{1,1},V4wc,'--r','LineWidth',3)
if lge==1

    plot(ts(1:size(I,1)),V4w,'g','LineWidth',3)
end
if eikonal==1
    plot(ts(1:size(Iei,1)),V4wei,'m','LineWidth',3)
end
plot(ts_ecgi(offset:end),V4e_f(offset:end),'k','LineWidth',3)
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)

text(180,15,'Lead V4','Fontsize',16)

grid on

i=11;
subplot('Position',[(.05+mod(i-1,3))/3 1.36-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),V5wc,'LineWidth',3)
hold on
%plot(LA{1,1},V5wc,'--r','LineWidth',3)
if lge==1

    plot(ts(1:size(I,1)),V5w,'g','LineWidth',3)
end
if eikonal==1
    plot(ts(1:size(Iei,1)),V5wei,'m','LineWidth',3)
end
plot(ts_ecgi(offset:end),V5e_f(offset:end),'k','LineWidth',3)
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');

set(gca,'Fontsize',16)
text(180,15,'Lead V5','Fontsize',16)

grid on

i=12;
subplot('Position',[(.05+mod(i-1,3))/3 1.36-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),V6wc,'LineWidth',3)
hold on
%plot(LA{1,1},V6wc,'--r','LineWidth',3)
if lge==1

    plot(ts(1:size(I,1)),V6w,'g','LineWidth',3)
end
if eikonal==1
    plot(ts(1:size(Iei,1)),V6wei,'m','LineWidth',3)
end
    plot(ts_ecgi(offset:end),V6e_f(offset:end),'k','LineWidth',3)
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
grid on
set(gca,'Fontsize',16)
text(180,3,'Lead V6','Fontsize',16)

grid on
if lge==1
    legend('control','LGE','eikonal', 'clinical recordings');
else
     legend('control', 'clinical recordings');

end
%Compute a goodness of fit to the ECGi data
% 
%  set(0,'RecursionLimit',1000);
% %[~,idx]=min(abs(V6w(i)*ones(size(V6e_f,1),1)-V6e_f));
%     y=[Ie_f,IIe_f,IIIe_f,aVRe_f,aVLe_f,aVFe_f,V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f];
%     fit=[];
%     j=1;
%     for x=[I,II,III,aVR,aVL,aVF,V1w,V2w,V3w,V4w,V5w,V6w]
%         sum_squared=0;
%         for i=1:size(x,1)
%             sum_squared = power(real(x(i))-real(y(i,j)),2)+sum_squared;
%         end
%         j
%         fit(j)=DiscreteFrechetDist(x,y(:,j));
%         j=j+1;
%     end
%     
%         y=[Ie_f,IIe_f,IIIe_f,aVRe_f,aVLe_f,aVFe_f,V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f];
%         j=1;
%     for x=[Ic,IIc,IIIc,aVRc,aVLc,aVFc,V1wc,V2wc,V3wc,V4wc,V5wc,V6wc]
% 
%         fit_control(j)=DiscreteFrechetDist(x,y(:,j));
%         j=j+1;
%     end
%     
  %try a Discrete Frechet distance as a measure of goodness of fit test to see if the values in the
  %simulation really come from a function which looks like body surface
  %potentials. This fuction measures the minimum distance required for a
  %dog and a walker to traverse their separate paths when a dog is on the
  %leash.



  figure()
  bar(1:12,vertcat(real(fit./fit_control), ones(1,12))','grouped')
  hold on
  set(gca,'XTick',[1,2,3,4,5,6,7,8,9,10,11,12])
  set(gca,'XTickLabel',{'I','II','III','aVR','aVL','aVF','V1','V2','V3','V4','V5','V6'});
  legend('LGE model', 'control model')
  set(gca,'Fontsize',21)
  ylabel('Normalised Frechet Distance')
  xlim([0 13])
%%
  figure()
  bar(1:12,fit,100,'filled')
  hold on
  scatter(1:12,fit_control,100,'filled')
  set(gca,'XTick',[1,2,3,4,5,6,7,8,9,10,11,12])
  set(gca,'XTickLabel',{'I','II','III','aVR','aVL','aVF','V1','V2','V3','V4','V5','V6'});
  legend('LGE model', 'control model')
  set(gca,'Fontsize',21)
  ylabel('Discrete Frechet Distance')
%%
