addpath(genpath('C:/Users/petnov/Dropbox/'));
addpath(genpath('C:/Users/petnov/Dropbox/qrs/'));

try
    load('D:\ARVC meshing automatic\patients\patient18\results\TestHeart18ControlECG\ecg.mat');
catch
    printf('no ecg file was found! Start computing the ecg...')
end

load('D:\ARVC meshing automatic\patients\patient09\results\TestHeart09Control\ecg.mat');
Icon=I_real;
IIcon=II_real;
IIIcon=III_real;
aVRcon=aVR_real;
aVLcon=aVL_real;
aVFcon=aVF_real;
V1wcon=V1_real;
V2wcon=V2_real;
V3wcon=V3_real;
V4wcon=V4_real;
V5wcon=V5_real;
V6wcon=V6_real;

Eikonal_light=1;
Eikonal_very_strong=0;

%load('D:\ARVC meshingautomatic\patients\patient09\results\TestHeart09LGEStrongCorrectRoots\ecg.mat');
if Eikonal_light==1
    ecg=load('D:\ARVC meshing automatic\patients\patient09\results\Eikonal\ecg150_NormalRootsLightFibLeftSeptumFrontStrong.mat');
    ecg=ecg.ecg_full;
    I=ecg(:,1);
    II=ecg(:,2);
    III=ecg(:,3);
    aVR=ecg(:,4);
    aVL=ecg(:,5);
    aVF=ecg(:,6);
    V1w=ecg(:,7);
    V2w=ecg(:,8);
    V3w=ecg(:,9);
    V4w=ecg(:,10);
    V5w=ecg(:,11);
    V6w=ecg(:,12);
else
load('D:\ARVC meshing automatic\patients\patient09\results\TestHeart09LGELight\ecg.mat');
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
end

if Eikonal_very_strong==1
    ecg=load('D:\ARVC meshing automatic\patients\patient09\results\Eikonal\ecg150_noseptlv.mat');
    ecg=ecg.ecg_full;
    Is=ecg(:,1);
    IIs=ecg(:,2);
    IIIs=ecg(:,3);
    aVRs=ecg(:,4);
    aVLs=ecg(:,5);
    aVFs=ecg(:,6);
    V1ws=ecg(:,7);
    V2ws=ecg(:,8);
    V3ws=ecg(:,9);
    V4ws=ecg(:,10);
    V5ws=ecg(:,11);
    V6ws=ecg(:,12);
else
    load('D:\ARVC meshing automatic\patients\patient09\results\TestHeart09LGEVeryStrong\ecg.mat');
    Iss=I_real;
    IIss=II_real;
    IIIss=III_real;
    aVRss=aVR_real;
    aVFss=aVF_real;
    aVLss=aVL_real;
    V1wss=V1_real;
    V2wss=V2_real;
    V3wss=V3_real;
    V4wss=V4_real;
    V5wss=V5_real;
    V6wss=V6_real;
end


% load('D:\ARVC meshing automatic\patients\patient09\results\TestHeart09LGELight\ecg.mat');
% I=I_real;
% II=II_real;
% III=III_real;
% aVR=aVR_real;
% aVL=aVL_real;
% aVF=aVF_real;
% V1w=V1_real;
% V2w=V2_real;
% V3w=V3_real;
% V4w=V4_real;
% V5w=V5_real;
% V6w=V6_real;


load('D:\ARVC meshing automatic\patients\patient09\results\TestHeart09LGEIntermediateCorrectRoots\ecg.mat');
Iint=I_real;
IIint=II_real;
IIIint=III_real;
aVRint=aVR_real;
aVLint=aVL_real;
aVFint=aVF_real;
V1wint=V1_real;
V2wint=V2_real;
V3wint=V3_real;
V4wint=V4_real;
V5wint=V5_real;
V6wint=V6_real;



%%
electrodes = load('D:\ARVC meshing automatic\patients\patient09\mpp\ECG_ELECTRODES.mat');


E = [ electrodes.ELECTRODES.LA ; electrodes.ELECTRODES.RA ; electrodes.ELECTRODES.LL ; electrodes.ELECTRODES.RL ; electrodes.ELECTRODES.V1 ; electrodes.ELECTRODES.V2 ; electrodes.ELECTRODES.V3 ; electrodes.ELECTRODES.V4 ; electrodes.ELECTRODES.V5 ; electrodes.ELECTRODES.V6 ];
E=E/10; %chaste is in cm
%plotMESH( MakeMesh( H , H.face ) , 'g' ); hplotMESH( T , 'nf' );
%hplot3d( E , '*m' ,'LineWidth',4,'MarkerSize',10 );
%axis(objbounds)

[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e,Ie_f,IIe_f,IIIe_f,aVRe_f,aVLe_f,aVFe_f,V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f,ts_ecgi]=ecgi_ECG_pat09('D:\Data_ARVC_raw_electrodes_vest\Data_ARVC\ARVC9\ARVC9_Baseline.mat','D:\ChrisReconstructions\ARVC9\ARVC9_Baseline.mat','D:\ARVC meshing automatic\patients\patient09\mpp\ECG_ELECTRODES.mat',9);

% ampl=max(real(V3e_f));
% Is=Is/max(Is)*ampl;
% IIs =IIs/max(IIs) * ampl;
% IIIs =IIIs/max(IIIs) *ampl;
% aVRs =aVRs/min(aVRs) *ampl;
% aVLs = aVLs/max(aVLs) *ampl;
% aVFs = aVFs/max(aVFs) * ampl;
% V1ws= V1ws/max(V1ws) * ampl;
% V2ws = V2ws/max(V2ws) * ampl;
% V3ws = V3ws/max(V3ws) * ampl;
% V4ws = V4ws/max(V4ws) *ampl;
% V5ws = V5ws/max(V5ws) * ampl;
% V6ws = V6ws/max(V6ws) * ampl;
% 


Iss=Iss/max(Iss)*ampl;
IIss =IIss/max(IIss) * ampl;
IIIss =IIIss/max(IIIss) *ampl;
aVRss =aVRss/min(aVRss) * ampl;
aVLss = aVLss/max(aVLss) * ampl;
aVFss = aVFss/max(aVFss) * ampl;
V1wss= V1wss/max(V1wss) *ampl;
V2wss = V2wss/max(V2wss) *ampl;
V3wss = V3wss/max(V3wss) * ampl;
V4wss = V4wss/max(V4wss) *ampl;
V5wss = V5wss/max(V5wss) * ampl;
V6wss = V6wss/max(V6wss) * ampl;

Icon=Icon/max(Icon)*ampl;
IIcon =IIcon/max(IIcon) *ampl;
IIIcon =IIIcon/max(IIIcon) * ampl;
aVRcon =aVRcon/max(-aVRcon) * ampl;
aVLcon = aVLcon/max(aVLcon) *ampl;
aVFcon = aVFcon/max(aVFcon) * ampl;
V1wcon= V1wcon/max(V1wcon) * ampl;
V2wcon = V2wcon/max(V2wcon) * ampl;
V3wcon = V3wcon/max(V3wcon) *ampl;
V4wcon = V4wcon/max(V4wcon) * ampl;
V5wcon = V5wcon/max(V5wcon) * ampl;
V6wcon = V6wcon/max(V6wcon) * ampl;



I=I/max(I)*ampl;
II =II/max(II) * ampl;
III =III/max(III) * ampl;
aVR =aVR/max(-aVR) * ampl;
aVL = aVL/max(aVL) *ampl;
aVF = aVF/max(aVF) * ampl;
V1w= V1w/max(V1w) * ampl;
V2w = V2w/max(V2w) * ampl;
V3w = V3w/max(V3w) *ampl;
V4w = V4w/max(V4w) *ampl;
V5w = V5w/max(V5w) * ampl;
V6w = V6w/max(V6w) * ampl;

% 
% Iint=Iint/max(Iint)*ampl;
% IIint =IIint/max(IIint) * ampl;
% IIIint =IIIint/max(IIIint) *ampl;
% aVRint =aVRint/max(-aVRint) * ampl;
% aVLint = aVLint/max(aVLint) * ampl;
% aVFint = aVFint/max(aVFint) *ampl;
% V1wint= V1wint/max(V1wint) * ampl;
% V2wint = V2wint/max(V2wint) *ampl;
% V3wint = V3wint/max(V3wint) * ampl;
% V4wint = V4wint/max(V4wint) * ampl;
% V5wint = V5wint/max(V5wint) * ampl;
% V6wint = V6wint/max(V6wint) * ampl;


%extract QRS biomarkers

medium = horzcat(V1wint,V2wint,V3wint,V4wint,V5wint,V6wint);
%strong = horzcat(V1ws,V2ws,V3ws,V4ws,V5ws,V6ws);
light=horzcat(V1w,V2w,V3w,V4w,V5w,V6w);
nolge=horzcat(V1wcon,V2wcon,V3wcon,V4wcon,V5wcon,V6wcon);

ecg=horzcat(V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f);
dt_ecgi=ts_ecgi(8)-ts_ecgi(7);
for i=1:size(medium,2)
            QRSs(i,1)=qrs_estimator(nolge(:,i));

    QRSs(i,2)=qrs_estimator(light(:,i));
 %   QRSs(i,3)=qrs_estimator(medium(:,i));
%    QRSs(i,4)=qrs_estimator(strong(:,i));

        QRSs(i,5)=qrs_estimator_bsp(ecg(:,i))*dt_ecgi;

end

%no significant difference in QRS durations
%now look at S wave area/peak


Swave=zeros(size(medium,2),10);
duration =zeros(size(medium,2),4);
for i=1:4
      [ Swave(i,1),Swave(i,2),duration(i,1)]=s_wave_area(nolge(:,i));

    [Swave(i,3),Swave(i,4),duration(i,2)]=s_wave_area(light(:,i));
  %[Swave(i,5),Swave(i,6),duration(i,3)]=s_wave_area(medium(:,i));
 %[ Swave(i,7),Swave(i,8),duration(i,4)]=s_wave_area(strong(:,i));

       [Swave(i,9),Swave(i,10),duration(i,5)]=s_wave_area_bsp(ecg(:,i),ts_ecgi);

end

ratio=zeros(size(medium,2),5);
for i=1:size(medium,2)
            [~,ratio(i,1)]=qrs_estimator(nolge(:,i));

    [~,ratio(i,2)]=qrs_estimator(light(:,i));
   %[~, ratio(i,3)]=qrs_estimator(medium(:,i));
   % [~,ratio(i,4)]=qrs_estimator(strong(:,i));

       [~, ratio(i,5)]=qrs_estimator_bsp(ecg(:,i));

end
% 
% 
% figure()
%   bar(1:4,duration(1:4,:)./QRSs(1:4,:),'grouped')
%   hold on
%   set(gca,'XTick',[1,2,3])
%   set(gca,'XTickLabel',{'V1','V2','V3','V4'});
%   legend( 'no lge','light LGE','medium LGE','strong LGE','recordings')
%   set(gca,'Fontsize',21)
%   ylabel('S wave as fraction of QRS Duration ')
% 
%   
% 
% figure()
%   bar(1:3,Swave(1:3,[2 4 6 8 10 ])./duration(1:3,:),'grouped')
%   hold on
%   set(gca,'XTick',[1,2,3])
%   set(gca,'XTickLabel',{'V1','V2','V3'});
%   legend('no lge','light LGE','medium LGE','strong LGE','recordings')
%   set(gca,'Fontsize',21)
%   ylabel('S wave area/(s wave duration x s wave amplitude) ')
% 
% 
% 
% figure()
%   bar(1:6,QRSs,'grouped')
%   hold on
%   set(gca,'XTick',[1,2,3,4,5,6])
%   set(gca,'XTickLabel',{'V1','V2','V3','V4','V5','V6'});
%   legend('no lge','light LGE','medium LGE','strong LGE','recordings')
%   set(gca,'Fontsize',21)
%   ylabel('QRS Duration (ms)')
% 
%   
%   
% figure()
%   bar(1:4,horzcat(Swave(1:4,1),Swave(1:4,3),Swave(1:4,5),Swave(1:4,7),Swave(1:4,9)),'grouped')
%   hold on
%   set(gca,'XTick',[1,2,3,4])
%   set(gca,'XTickLabel',{'V1','V2','V3','V4'});
%   legend('no LGE','light LGE','medium LGE','strong LGE + septal LGE eikonal','recordings')
%   set(gca,'Fontsize',21)
%   ylabel('Area under S-wave (Vms)')
%   
%   
%     figure()
%   bar(1:4,horzcat(ratio(1:4,1),ratio(1:4,2),ratio(1:4,3),ratio(1:4,4),ratio(1:4,5)),'grouped')
%   hold on
%   set(gca,'XTick',[1,2,3,4])
%   set(gca,'XTickLabel',{'V1','V2','V3','V4'});
%   legend('no LGE','light LGE','medium LGE','strong LGE', 'recordings')
%   set(gca,'Fontsize',21)
%   ylabel('R peak to S peak ratio')
%   
%   
%   

  close all


lge=0;
lge2=0;
lge_strong=0;
lge_very_strong=0;
control=0;
data=1;
offset=100;%idx
ts=1:size(Icon,1)+offset;

ts_ecgi=ts_ecgi-ts_ecgi(offset);
figure()
grid on
i=1;
subplot('Position',[(.15+mod(i-1,3))/3 1.125-(ceil(i/3))/3 1/4 1/5])
hold all

if control==1
        plot(ts(offset:size(Icon,1)+offset-1),Icon,'b','LineWidth',3)
end
if lge==1
    plot(ts(offset:size(I,1)+offset-1),I,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(offset:size(I,1)),Iint,'g','LineWidth',3)
end
if lge_strong==1
    plot(ts(offset:size(I,1)),Is,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if lge_very_strong==1
    plot(ts(1:end),size(I,1),'Color','m','LineWidth',3)
end
%plot(ts_,'--r','LineWidth',3)
if data==1
    plot(ts_ecgi(offset:end),Ie_f(offset:end),'Color',[0.91 0.41 0.17],'LineWidth',3)
end


plot(ts_ecgi(1:end),Ie_f(1:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
text(200,25,'Lead I','Fontsize',16)
grid on
set(gca,'Fontsize',16)
ylim([-70 35])
i=2;
subplot('Position',[(0.05+mod(i-1,3))/3 1.125-(ceil(i/3))/3 1/4 1/5])
hold all

if control==1
        plot(ts(offset:size(Icon,1)+offset-1),IIcon,'b','LineWidth',3)
end
if lge==1
    plot(ts(offset:size(I,1)+offset-1),II,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(offset:end),IIint,'g','LineWidth',3)
end
if lge_strong==1

    plot(ts(offset:end),IIs,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if lge_very_strong==1
    plot(ts(offset:end),IIss,'Color','m','LineWidth',3)
end
if data==1

    plot(ts_ecgi(1:end),IIe_f(1:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end
ylim([-70 40])

xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
text(200,25,'Lead II','Fontsize',16)
grid on
set(gca,'Fontsize',16)

i=3;
subplot('Position',[(.05+mod(i-1,3))/3 1.125-(ceil(i/3))/3 1/4 1/5])
hold all

if control==1
        plot(ts(offset:size(Icon,1)+offset-1),IIIcon,'b','LineWidth',3)
end
if lge==1
    plot(ts(offset:size(I,1)+offset-1),III,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(offset:size(IIIint,1)),IIIint,'g','LineWidth',3)
end
if lge_strong==1

    plot(ts(offset:end),IIIs,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if lge_very_strong==1
    plot(ts(offset:end),IIIss,'Color','m','LineWidth',3)
end
if data==1

    plot(ts_ecgi(1:end),IIIe_f(1:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)
ylim([-70 40])

text(200,7,'Lead III','Fontsize',16)
i=4;
subplot('Position',[(.15+mod(i-1,3))/3 1.195-(ceil(i/3))/3 1/4 1/5])
hold all
ylim([-70 40])
xlim([100 250])

if control==1
        plot(ts(offset:size(Icon,1)+offset-1),aVRcon,'b','LineWidth',3)
end
if lge==1
    plot(ts(offset:size(I,1)+offset-1),aVR,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(offset:size(aVRint,1)),aVRint,'g','LineWidth',3)
end
if lge_strong==1

plot(ts(offset:end),aVRs,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if lge_very_strong==1
    plot(ts(offset:end),aVRss,'Color','m','LineWidth',3)
end
plot(ts_ecgi(1:end),aVRe_f(1:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)

xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
text(200,5,'Lead aVR','Fontsize',16)
grid on
set(gca,'Fontsize',16)
i=5;
subplot('Position',[(.05+mod(i-1,3))/3 1.195-(ceil(i/3))/3 1/4 1/5])
hold all
ylim([-70 40])
xlim([100 250])

if control==1
        plot(ts(offset:size(Icon,1)+offset-1),aVLcon,'b','LineWidth',3)
end
if lge==1
    plot(ts(offset:size(I,1)+offset-1),aVL,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(offset:end),aVLint,'g','LineWidth',3)
end
if lge_strong==1
plot(ts(offset:end),aVLs,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if lge_very_strong==1
    plot(ts(offset:end),aVLss,'Color','m','LineWidth',3)
end
if data==1

    plot(ts_ecgi(1:end),aVLe_f(1:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
text(200,10,'Lead aVL','Fontsize',16)
xlim([100 250])


grid on
set(gca,'Fontsize',16)
i=6;
subplot('Position',[(.05+mod(i-1,3))/3 1.195-(ceil(i/3))/3 1/4 1/5])
hold all

if control==1
        plot(ts(offset:size(Icon,1)+offset-1),aVFcon,'b','LineWidth',3)
end
if lge==1
    plot(ts(offset:size(I,1)+offset-1),aVF,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(offset:end),aVFint,'g','LineWidth',3)
end
if lge_strong==1

plot(ts(offset:size(aVFs,1)),aVFs,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if lge_very_strong==1
    plot(ts(offset:end),aVFss,'Color','m','LineWidth',3)
end
if data==1

    plot(ts_ecgi(1:end),aVFe_f(1:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end
ylim([-70 40])

xlim([100 250])

xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)
text(200,15,'Lead aVF','Fontsize',16)


grid on

i=7;
subplot('Position',[(.15+mod(i-1,3))/3 1.279-(ceil(i/3))/3 1/4 1/5])
hold all

if control==1
        plot(ts(offset:size(Icon,1)+offset-1),V1wcon,'b','LineWidth',3)
end
if lge==1
    plot(ts(offset:size(I,1)+offset-1),V1w,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(offset:end),V1wint,'g','LineWidth',3)
end
if lge_strong==1

plot(ts(offset:size(V1ws,1)),V1ws,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if lge_very_strong==1
    plot(ts(offset:end),V1wss,'Color','m','LineWidth',3)
end
if data==1

    plot(ts_ecgi(1:end),V1e_f(1:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
ylim([-70 40])
xlim([100 250])

vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)
text(165,22,'Lead V1','Fontsize',16)


grid on

i=8;
subplot('Position',[(.05+mod(i-1,3))/3 1.279-(ceil(i/3))/3 1/4 1/5])
hold all

if control==1
        plot(ts(offset:size(Icon,1)+offset-1),V2wcon,'b','LineWidth',3)
end
if lge==1
    plot(ts(offset:size(I,1)+offset-1),V2w,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(offset:end),V2wint,'g','LineWidth',3)
end
if lge_strong==1

plot(ts(offset:end),V2ws,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if lge_very_strong==1
    plot(ts(offset:end),V2wss,'Color','m','LineWidth',3)
end
if data==1

     plot(ts_ecgi(1:end),V2e_f(1:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end    
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)
ylim([-70 40])
xlim([100 250])

text(165,22,'Lead V2','Fontsize',16)


grid on

i=9;
subplot('Position',[(0.05+mod(i-1,3))/3 1.279-(ceil(i/3))/3 1/4 1/5])
hold all

if control==1
        plot(ts(offset:size(Icon,1)+offset-1),V3wcon,'b','LineWidth',3)
end
if lge==1
    plot(ts(offset:size(I,1)+offset-1),V3w,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(offset:end),V3wint,'g','LineWidth',3)
end
if lge_strong==1

plot(ts(offset:end),V3ws,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if lge_very_strong==1
    plot(ts(offset:end),V3wss,'Color','m','LineWidth',3)
end
if data==1

    plot(ts_ecgi(1:end),V3e_f(1:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)
ylim([-70 40])
xlim([100 250])

text(165,22,'Lead V3','Fontsize',16)


grid on

i=10;
subplot('Position',[(.15+mod(i-1,3))/3 1.36-(ceil(i/3))/3 1/4 1/5])
hold all

if control==1
        plot(ts(offset:size(Icon,1)+offset-1),V4wcon,'b','LineWidth',3)
end
if lge==1
    plot(ts(offset:size(I,1)+offset-1),V4w,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(offset:end),V4wint,'g','LineWidth',3)
end
if lge_strong==1

plot(ts(1:end),V4ws,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if lge_very_strong==1
    plot(ts(offset:end),V4wss,'Color','m','LineWidth',3)
end
if data==1

    plot(ts_ecgi(1:end),V4e_f(1:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)

text(165,22,'Lead V4','Fontsize',16)
ylim([-70 40])
xlim([100 250])

grid on

i=11;
subplot('Position',[(.05+mod(i-1,3))/3 1.36-(ceil(i/3))/3 1/4 1/5])
hold all

if control==1
        plot(ts(offset:size(Icon,1)+offset-1),V5wcon,'b','LineWidth',3)
end
if lge==1
    plot(ts(offset:size(I,1)+offset-1),V5w,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(offset:end),V5wint,'g','LineWidth',3)
end
if lge_strong==1

plot(ts(1:end),V5ws,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if lge_very_strong==1
    plot(ts(offset:end),V5wss,'Color','m','LineWidth',3)
end
if data==1

    plot(ts_ecgi(1:end),V5e_f(1:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');

set(gca,'Fontsize',16)
text(165,22,'Lead V5','Fontsize',16)
ylim([-70 40])
xlim([100 250])

grid on

i=12;
subplot('Position',[(.05+mod(i-1,3))/3 1.36-(ceil(i/3))/3 1/4 1/5])
hold all

if control==1
        plot(ts(offset:size(Icon,1)+offset-1),V6wcon,'b','LineWidth',3)
end
if lge==1
    plot(ts(offset:size(I,1)+offset-1),V6w,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(offset:end),V6wint,'g','LineWidth',3)
end
if lge_strong==1

plot(ts(offset:end),V6ws,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if lge_very_strong==1
    plot(ts(offset:end),V6wss,'Color','m','LineWidth',3)
end
if data==1

    plot(ts_ecgi(1:end),V6e_f(1:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
grid on
set(gca,'Fontsize',16)
text(165,22,'Lead V6','Fontsize',16)
ylim([-70 40])
xlim([100 250])

grid on
if lge==1 && lge2==1
    legend('sim no LGE','sim lge light', 'sim lge intermediate','sim lge light LV anterior + intermediate lge septal + eikonal','sim lge very strong','recorded ECG');
else
     legend('control sim', 'septal and lv anterior LGE','recorded ECG');

end
    %legend( 'sim lge intermediate','sim lge strong','sim lge very strong','recorded ECG');




%only precordial leads below




lge=0;
lge2=0;
lge_strong=1;

control=0;
data=1;
ts=1:208;
offset=1;%idx
ts_ecgi=ts_ecgi-ts_ecgi(offset);
figure()
grid on


i=1;
subplot('Position',[(.05+mod(i-1,3))/3 1.125-(ceil(i/3))/3 1/4 1/5])
hold all

if control==1
        plot(ts(1:size(V1wcon,1)),V1wcon,'b','LineWidth',3)
end
if lge==1
    plot(ts(1:size(V1w,1)),V1w,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(1:size(V1wint,1)),V1wint,'g','LineWidth',3)
end
if lge_strong==1

plot(ts(1:size(V1ws,1)),V1ws,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if data==1

    plot(ts_ecgi(offset:end),V1e_f(offset:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
ylim([-70 40])

vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)
text(85,71,'Lead V1','Fontsize',16)


grid on

i=2;
subplot('Position',[(.05+mod(i-1,3))/3 1.125-(ceil(i/3))/3 1/4 1/5])
hold all

if control==1
        plot(ts(1:size(V2wcon,1)),V2wcon,'b','LineWidth',3)
end
if lge==1
    plot(ts(1:size(V2w,1)),V2w,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(1:size(V2wint,1)),V2wint,'g','LineWidth',3)
end
if lge_strong==1

plot(ts(1:size(V2ws,1)),V2ws,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if data==1

     plot(ts_ecgi(offset:end),V2e_f(offset:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end    
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)
ylim([-70 40])

text(85,70,'Lead V2','Fontsize',16)


grid on

i=3;
subplot('Position',[(.05+mod(i-1,3))/3 1.125-(ceil(i/3))/3 1/4 1/5])
hold all

if control==1
        plot(ts(1:size(V3wcon,1)),V3wcon,'b','LineWidth',3)
end
if lge==1
    plot(ts(1:size(V3w,1)),V3w,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(1:size(V3wint,1)),V3wint,'g','LineWidth',3)
end
if lge_strong==1

plot(ts(1:size(V3ws,1)),V3ws,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if data==1

    plot(ts_ecgi(offset:end),V3e_f(offset:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)
ylim([-70 40])

text(85,70,'Lead V3','Fontsize',16)


grid on

i=4;
subplot('Position',[(.05+mod(i-1,3))/3 1.195-(ceil(i/3))/3 1/4 1/5])
hold all

if control==1
        plot(ts(1:size(V4wcon,1)),V4wcon,'b','LineWidth',3)
end
if lge==1
    plot(ts(1:size(V4w,1)),V4w,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(1:size(V4wint,1)),V4wint,'g','LineWidth',3)
end
if lge_strong==1

plot(ts(1:size(V4ws,1)),V4ws,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if data==1

    plot(ts_ecgi(offset:end),V4e_f(offset:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)

text(85,70,'Lead V4','Fontsize',16)
ylim([-70 40])

grid on

i=5;
subplot('Position',[(.05+mod(i-1,3))/3 1.195-(ceil(i/3))/3 1/4 1/5])
hold all

if control==1
        plot(ts(1:size(V5wcon,1)),V5wcon,'b','LineWidth',3)
end
if lge==1
    plot(ts(1:size(V5w,1)),V5w,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(1:size(V5wint,1)),V5wint,'g','LineWidth',3)
end
if lge_strong==1

plot(ts(1:size(V5ws,1)),V5ws,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if data==1

    plot(ts_ecgi(offset:end),V5e_f(offset:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');

set(gca,'Fontsize',16)
text(85,70,'Lead V5','Fontsize',16)
ylim([-70 40])

grid on

i=6;
subplot('Position',[(.05+mod(i-1,3))/3 1.195-(ceil(i/3))/3 1/4 1/5])
hold all

if control==1
        plot(ts(1:size(V6wcon,1)),V6wcon,'b','LineWidth',3)
end
if lge==1
    plot(ts(1:size(V6w,1)),V6w,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(1:size(V6wint,1)),V6wint,'g','LineWidth',3)
end
if lge_strong==1

plot(ts(1:size(V6ws,1)),V6ws,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if data==1
    plot(ts_ecgi(offset:end),V6e_f(offset:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
grid on
set(gca,'Fontsize',16)
text(85,70,'Lead V6','Fontsize',16)
ylim([-70 40])

grid on
if lge==1 && lge2==1
    legend('sim no LGE','sim lge light', 'sim lge intermediate','sim lge strong','recorded ECG');
else
     legend('control sim', 'recorded ECG');

end
%Compute a goodness of fit to the ECGi data


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
%     for x=[Is,IIs,IIIs,aVRs,aVLs,aVFs,V1ws,V2ws,V3ws,V4ws,V5ws,V6ws]
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
  
  
  %quantify the biomarkers on the QRS
  
  
  
  
%   
%   
% a=real(fit./fit_control);
%   figure()
%   bar(1:3,vertcat([a(1),a(7),a(8)], ones(1,3))','grouped')
%   hold on
%   set(gca,'XTick',[1,2,3,4,5,6,7,8,9,10,11,12])
%   set(gca,'XTickLabel',{'I','V1','V2'});
%   legend('LGE model', 'control model')
%   set(gca,'Fontsize',21)
%   ylabel('Normalised Frechet Distance')
%   xlim([0 13])
%   
%     figure()
%   bar(1:12,vertcat(real(fit(1./fit_control), ones(1,12))','grouped'))
%   hold on
%   set(gca,'XTick',[1,2,3,4,5,6,7,8,9,10,11,12])
%   set(gca,'XTickLabel',{'I','II','III','aVR','aVL','aVF','V1','V2','V3','V4','V5','V6'});
%   legend('LGE model', 'control model')
%   set(gca,'Fontsize',21)
%   ylabel('Normalised Frechet Distance')
%   xlim([0 13])
%%
  %make a movie of lead V4 in order to investigate epsilon wave

% 
% for i=1:71
%     h = figure();
% 
%     hold on
%        plot(ts(1:size(I,1)),V4w,'Color',[0.54 0.82 0.99],'LineWidth',3)
%        scatter(ts(i),V4w(i),40,'filled')
% 
%     
%     F(i) = getframe(h);
%     
% end
% 
% 
% fig = figure;
% movie(fig,F(1:71),1)
% 
% v = VideoWriter('D:\ARVC meshing automatic\patients\patient09\results\TestHeart09LGE2\ecg_vid','Motion JPEG AVI');
% 
% open(v)
% writeVideo(v,F(1:71));
% 
% close(v)

%%
