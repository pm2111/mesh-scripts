addpath(genpath('C:/Users/petnov/Dropbox/'));
addpath(genpath('C:/Users/petnov/Dropbox/qrs/'));
addpath(genpath('C:/Users/petnov/Dropbox/sim 3d scripts/'));
addpath('C:\Users\peter\Dropbox\sim3d scripts\');

try
    load('D:\ARVC meshing automatic\patients\patient18\results\TestHeart18ControlECG\ecg.mat');
catch
    printf('no ecg file was found! Start computing the ecg...')
end
Eikonal_light=1;
Eikonal_int=1;
Eikonal_strong=1;
Eikonal_very_strong=0;
Eikonal_control=1;
if Eikonal_control==1
     ecg=load('D:\ARVC meshing automatic\patients\patient06\results\Eikonal\exp12\17.0x_0f_1.0fiber_roots_3_no_rv_pur_1\17.0x_0f_1.0fiber_roots_3_no_rv_pur_1_NEW_pred_ATMap.mat');
    ecg=ecg.ecg_full;
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
else
load('D:\ARVC meshing automatic\patients\patient06\results\TestHeart06Node2\ecg.mat');
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
end

%load('D:\ARVC meshingautomatic\patients\patient09\results\TestHeart09LGEStrongCorrectRoots\ecg.mat');
if Eikonal_strong==1
    ecg=load('D:\ARVC meshing automatic\patients\patient06\results\Eikonal\exp8_ECGs\17.0x_3f_1.0fiber_NEW_pred_ATMap.csvpseudo_ECG.mat');
    ecg=ecg.vars;
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
load('D:\ARVC meshing automatic\patients\patient06\results\TestHeart06Node2LGE\ecg.mat');
    Is=I_real;
    IIs=II_real;
    IIIs=III_real;
    aVRs=aVR_real;
    aVFs=aVF_real;
    aVLs=aVL_real;
    V1ws=V1_real;
    V2ws=V2_real;
    V3ws=V3_real;
    V4ws=V4_real;
    V5ws=V5_real;
    V6ws=V6_real;
end

if Eikonal_very_strong==1
    ecg=load('D:\ARVC meshing automatic\patients\patient09\results\Eikonal\ecg150SmallFibStrong.mat');
    ecg=ecg.vars;
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
load('D:\ARVC meshing automatic\patients\patient06\results\TestHeart06Node2LGE\ecg.mat');
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

[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e,Ie_f,IIe_f,IIIe_f,aVRe_f,aVLe_f,aVFe_f,V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f,ts_ecgi]=ecgi_ECG06('D:\Data_ARVC_raw_electrodes_vest\Data_ARVC\ARVC6\ARVC6_Baseline.mat','D:\ChrisReconstructions\ARVC6\ARVC6_Baseline.mat','D:\ARVC meshing automatic\patients\patient06\mpp\ECG_ELECTRODES.mat',6);

if Eikonal_light==1
    try
    ecg= csvread('C:\Users\peter\Dropbox\sim3D\eikonal06_fine_pred_PseudoECG_200_2x_1f.csv')';
    catch
   ecg= load('D:\ARVC meshing automatic\patients\patient06\results\Eikonal\exp8_ECGs\17.0x_1f_1.0fiber_NEW_pred_ATMap.csvpseudo_ECG.mat')';
        ecg=ecg.vars;
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
    end 
   %ecg=load('D:\ARVC meshing automatic\patients\patient06\results\Eikonal\ecg100SmallFibLight.mat');
    %ecg=ecg.ecg_full;
    else
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
end



if Eikonal_int==1
    ecg=load('D:\ARVC meshing automatic\patients\patient06\results\Eikonal\exp8_ECGs\17.0x_2f_1.0fiber_NEW_pred_ATMap.csvpseudo_ECG.mat');
        ecg=ecg.vars;
    Iint=ecg(:,1);
    IIint=ecg(:,2);
    IIIint=ecg(:,3);
    aVRint=ecg(:,4);
    aVLint=ecg(:,5);
    aVFint=ecg(:,6);
    V1wint=ecg(:,7);
    V2wint=ecg(:,8);
    V3wint=ecg(:,9);
    V4wint=ecg(:,10);
    V5wint=ecg(:,11);
    V6wint=ecg(:,12);
else
load('D:\ARVC meshing automatic\patients\patient06\results\TestHeart06Node2LGE\ecg.mat');
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

end



%
electrodes = load('D:\ARVC meshing automatic\patients\patient06\mpp\ECG_ELECTRODES.mat');


E = [ electrodes.ELECTRODES.LA ; electrodes.ELECTRODES.RA ; electrodes.ELECTRODES.LL ; electrodes.ELECTRODES.RL ; electrodes.ELECTRODES.V1 ; electrodes.ELECTRODES.V2 ; electrodes.ELECTRODES.V3 ; electrodes.ELECTRODES.V4 ; electrodes.ELECTRODES.V5 ; electrodes.ELECTRODES.V6 ];
E=E/10; %chaste is in cm




offset=200;%idx

ampl=max(real(V5e_f(offset:offset+200)));
ampl_min=min(real(V5e_f(1:600)));

ampl_prev= max(V5ws);
ampl_prev_min=min(V5ws);
Is=Is/ampl_prev*ampl;
IIs =IIs/ampl_prev * ampl;
IIIs =IIIs/ampl_prev *ampl;
aVRs =aVRs/ampl_prev *ampl;
aVLs = aVLs/ampl_prev *ampl;
aVFs = aVFs/ampl_prev * ampl;
if max(V1ws)<=0
    V1ws= V1ws/ampl_prev_min * ampl_min;
else
        V1ws= V1ws/ampl_prev * ampl;
end
V2ws = V2ws/max(V2ws) * ampl;
if max(V3ws)<=0
    V3ws= V3ws/ampl_prev_min * ampl_min;
else
    V3ws = V3ws/ampl_prev * ampl;
end
V4ws = V4ws/ampl_prev *ampl;
V5ws = V5ws/ampl_prev * ampl;
V6ws = V6ws/ampl_prev * ampl;


ampl_prev= max(V5wss);
ampl_prev_min=min(V5wss);
Iss=Iss/ampl_prev*ampl;
IIss =IIss/ampl_prev * ampl;
IIIss =IIIss/ampl_prev *ampl;
aVRss =aVRss/ampl_prev * ampl;
aVLss = aVLss/ampl_prev * ampl;
aVFss = aVFss/ampl_prev * ampl;
V1wss= V1wss/ampl_prev *ampl;
V2wss = V2wss/ampl_prev *ampl;
V3wss = V3wss/ampl_prev * ampl;
V4wss = V4wss/ampl_prev *ampl;
V5wss = V5wss/ampl_prev * ampl;
V6wss = V6wss/ampl_prev * ampl;

ampl_prev= max(V5wcon);

Icon=Icon/ampl_prev*ampl;
IIcon =IIcon/ampl_prev *ampl;
IIIcon =IIIcon/ampl_prev * ampl;
aVRcon =aVRcon/ampl_prev * ampl;
aVLcon = aVLcon/ampl_prev *ampl;
aVFcon = aVFcon/ampl_prev * ampl;
V1wcon= V1wcon/ampl_prev * ampl;
V2wcon = V2wcon/ampl_prev * ampl;
V3wcon = V3wcon/ampl_prev *ampl;
V4wcon = V4wcon/ampl_prev * ampl;
V5wcon = V5wcon/ampl_prev * ampl;
V6wcon = V6wcon/ampl_prev * ampl;


ampl_prev= max(V5w);
ampl_prev_min= abs(min(V5w));

I=I/ampl_prev*ampl;
II =II/ampl_prev * ampl;
III =III/ampl_prev* ampl;
aVR =aVR/ampl_prev * ampl;
aVL = aVL/ampl_prev*ampl;
aVF = aVF/ampl_prev * ampl;
if max(V1w) <=0
    V1w= V1w/ampl_prev_min * ampl_min;
else
     V1w= V1w/ampl_prev * ampl;
end
V2w = V2w/ampl_prev * ampl;
if max(V3w) <=0
    V3w= V3w/ampl_prev_min * ampl_min;
else
     V3w= V3w/ampl_prev * ampl;
end
endV4w = V4w/ampl_prev *ampl;
V5w = V5w/ampl_prev * ampl;
V6w = V6w/ampl_prev * ampl;

ampl_prev=max(V5wint);
Iint=Iint/ampl_prev*ampl;
IIint =IIint/ampl_prev * ampl;
IIIint =IIIint/ampl_prev *ampl;
aVRint =aVRint/ampl_prev * ampl;
aVLint = aVLint/ampl_prev * ampl;
aVFint = aVFint/ampl_prev *ampl;
V1wint= V1wint/ampl_prev * ampl;
V2wint = V2wint/ampl_prev *ampl;
V3wint = V3wint/ampl_prev * ampl;
V4wint = V4wint/ampl_prev * ampl;
V5wint = V5wint/ampl_prev * ampl;
V6wint = V6wint/ampl_prev * ampl;


%extract QRS biomarkers

medium = horzcat(V1wint,V2wint,V3wint,V4wint,V5wint,V6wint);
strong = horzcat(V1ws,V2ws,V3ws,V4ws,V5ws,V6ws);
light=horzcat(V1w,V2w,V3w,V4w,V5w,V6w);
nolge=horzcat(V1wcon,V2wcon,V3wcon,V4wcon,V5wcon,V6wcon);

ecg=horzcat(V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f);
dt_ecgi=ts_ecgi(8)-ts_ecgi(7);
R_peak=[];
Rs=[];
for i=1:size(medium,2)
%             [QRSs(i,1)]=qrs_estimator(nolge(:,i));
% 
%     %QRSs(i,2)=qrs_estimator(light(:,i));
%     QRSs(i,3)=qrs_estimator(medium(:,i));
%     QRSs(i,4)=qrs_estimator(strong(:,i));
% 
       [ QRSs(i,5),~,~,~,~,~,~,~,R_peak(i),ST_el(i)]=qrs_estimator_bsp(ecg(:,i));

end
QRSs(:,5)=QRSs(:,5)*dt_ecgi;
%no significant difference in QRS durations
%now look at S wave area/peak


Swave=zeros(size(medium,2),10);
duration =zeros(size(medium,2),4);
for i=1:4
      [ Swave(i,1),Swave(i,2),duration(i,1)]=s_wave_area(nolge(:,i));

    %[Swave(i,3),Swave(i,4),duration(i,2)]=s_wave_area(light(:,i));
  %[Swave(i,5),Swave(i,6),duration(i,3)]=s_wave_area(medium(:,i));
% [ Swave(i,7),Swave(i,8),duration(i,4)]=s_wave_area(strong(:,i));

       [Swave(i,9),Swave(i,10),duration(i,5)]=s_wave_area_bsp(ecg(:,i),ts_ecgi);

end

ratio=zeros(size(medium,2),5);
for i=1:size(medium,2)
           [~,ratio(i,1),~,~,~,~,~,R_peaks(i)]=qrs_estimator(nolge(:,i));

    [~,ratio(i,2),~,~,~,~,~,R_peaks2(i)]=qrs_estimator(light(:,i));
   [~, ratio(i,3),~,~,~,~,~,R_peaks3(i)]=qrs_estimator(medium(:,i));
  %  [~,ratio(i,4)]=qrs_estimator(strong(:,i));

      [~, ratio(i,5)]=qrs_estimator_bsp(ecg(:,i));

end
Rs=vertcat(Rs,R_peak,R_peaks,R_peaks2,R_peaks3);


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
  idx=find(ratio>10);
  for i = 1:size(idx)
      j=idx(i);
   if (j<7)
        ratio(j)=R_peaks(j);
   elseif (j>6) & (j<13)
         ratio(j)=R_peaks2(j-6);
   elseif (j>12) & (j<19)
         ratio(j)=R_peaks3(j-12);
       
   else
       
       ratio(j)=R_peak(j-24);
   end
  end
  
    figure()
  bar(1:6,horzcat(ratio(1:6,1)/max(ratio(1:6,1)),ratio(1:6,2)/max(ratio(1:6,2)),ratio(1:6,3)/max(ratio(1:6,3)),ratio(1:6,5)/max(ratio(1:6,5))),'grouped')
  hold on
  set(gca,'XTick',[1,2,3,4])
  set(gca,'XTickLabel',{'V1','V2','V3','V4','V5','V6'});
  legend('no LGE','light LGE','intermediate LGE','recordings')
  set(gca,'Fontsize',21)
  ylabel('R peak to S peak ratio')
  
%   
  



lge=1;
lge2=1;
lge_strong=1;
lge_very_strong=0;
control=1;
data=1;
ts=1:size(V1w,1)+900;
ts_ecgi=ts_ecgi-ts_ecgi(offset);
figure()
grid on
i=1;
subplot('Position',[(.15+mod(i-1,3))/3 1.125-(ceil(i/3))/3 1/4 1/5])
hold all

if control==1
        plot(ts(1:size(Icon,1)),Icon,'b','LineWidth',3)
end
if lge==1
    plot(ts(1:size(I,1)),I,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(1:size(Iint,1)),Iint,'g','LineWidth',3)
end
if lge_strong==1
    plot(ts(1:size(Is,1)),Is,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if lge_very_strong==1
    plot(ts(1:size(Iss,1)),Iss,'Color','m','LineWidth',3)
end
%plot(ts_,'--r','LineWidth',3)
if data==1
    plot(ts_ecgi(offset:end),Ie_f(offset:end),'Color',[0.91 0.41 0.17],'LineWidth',3)
end


plot(ts_ecgi(offset:end),Ie_f(offset:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
text(200,25,'Lead I','Fontsize',16)
grid on
ylim([min(real(ecg(:,3)))-5   max(real(ecg(:,3)))+5])

set(gca,'Fontsize',16)
i=2;
subplot('Position',[(0.05+mod(i-1,3))/3 1.125-(ceil(i/3))/3 1/4 1/5])
hold all

if control==1
        plot(ts(1:size(IIcon,1)),IIcon,'b','LineWidth',3)
end
if lge==1
    plot(ts(1:size(II,1)),II,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(1:size(IIint,1)),IIint,'g','LineWidth',3)
end
if lge_strong==1

    plot(ts(1:size(IIs,1)),IIs,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if lge_very_strong==1
    plot(ts(1:size(Iss,1)),IIss,'Color','m','LineWidth',3)
end
if data==1

    plot(ts_ecgi(offset:end),IIe_f(offset:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end

xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
text(200,25,'Lead II','Fontsize',16)
grid on
set(gca,'Fontsize',16)
ylim([min(real(ecg(:,3)))-5   max(real(ecg(:,3)))+5])

i=3;
subplot('Position',[(.05+mod(i-1,3))/3 1.125-(ceil(i/3))/3 1/4 1/5])
hold all

if control==1
        plot(ts(1:size(IIIcon,1)),IIIcon,'b','LineWidth',3)
end
if lge==1
    plot(ts(1:size(III,1)),III,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(1:size(IIIint,1)),IIIint,'g','LineWidth',3)
end
if lge_strong==1

    plot(ts(1:size(IIIs,1)),IIIs,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if lge_very_strong==1
    plot(ts(1:size(Iss,1)),IIIss,'Color','m','LineWidth',3)
end
if data==1

    plot(ts_ecgi(offset:end),IIIe_f(offset:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)
ylim([min(real(ecg(:,3)))-5   max(real(ecg(:,3)))+5])

text(200,7,'Lead III','Fontsize',16)
i=4;
subplot('Position',[(.15+mod(i-1,3))/3 1.195-(ceil(i/3))/3 1/4 1/5])
hold all

if control==1
        plot(ts(1:size(aVRcon,1)),aVRcon,'b','LineWidth',3)
end
if lge==1
    plot(ts(1:size(aVR,1)),aVR,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(1:size(aVRint,1)),aVRint,'g','LineWidth',3)
end
if lge_strong==1

plot(ts(1:size(aVRs,1)),aVRs,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if lge_very_strong==1
    plot(ts(1:size(Iss,1)),aVRss,'Color','m','LineWidth',3)
end
plot(ts_ecgi(offset:end),aVRe_f(offset:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)

xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
text(200,5,'Lead aVR','Fontsize',16)
grid on
set(gca,'Fontsize',16)
i=5;
subplot('Position',[(.05+mod(i-1,3))/3 1.195-(ceil(i/3))/3 1/4 1/5])
hold all
ylim([min(real(ecg(:,3)))-5   max(real(ecg(:,3)))+5])

if control==1
        plot(ts(1:size(aVLcon,1)),aVLcon,'b','LineWidth',3)
end
if lge==1
    plot(ts(1:size(aVL,1)),aVL,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(1:size(aVLint,1)),aVLint,'g','LineWidth',3)
end
if lge_strong==1
plot(ts(1:size(aVLs,1)),aVLs,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if lge_very_strong==1
    plot(ts(1:size(Iss,1)),aVLss,'Color','m','LineWidth',3)
end
if data==1

    plot(ts_ecgi(offset:end),aVLe_f(offset:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
text(200,10,'Lead aVL','Fontsize',16)
ylim([min(real(ecg(:,3)))-5   max(real(ecg(:,3)))+5])


grid on
set(gca,'Fontsize',16)
i=6;
subplot('Position',[(.05+mod(i-1,3))/3 1.195-(ceil(i/3))/3 1/4 1/5])
hold all

if control==1
        plot(ts(1:size(aVFcon,1)),aVFcon,'b','LineWidth',3)
end
if lge==1
    plot(ts(1:size(aVF,1)),aVF,'Color',[0.54 0.82 0.99],'LineWidth',3)
end
if lge2==1
    plot(ts(1:size(aVFint,1)),aVFint,'g','LineWidth',3)
end
if lge_strong==1

plot(ts(1:size(aVFs,1)),aVFs,'Color',[0.91 0.41 0.17],'LineWidth',3)
end
if lge_very_strong==1
    plot(ts(1:size(Iss,1)),aVFss,'Color','m','LineWidth',3)
end
if data==1

    plot(ts_ecgi(offset:end),aVFe_f(offset:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end


xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)
text(200,15,'Lead aVF','Fontsize',16)

ylim([min(real(ecg(:,3)))-5   max(real(ecg(:,3)))+5])

grid on

i=7;
subplot('Position',[(.15+mod(i-1,3))/3 1.279-(ceil(i/3))/3 1/4 1/5])
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
if lge_very_strong==1
    plot(ts(1:size(Iss,1)),V1wss,'Color','m','LineWidth',3)
end
if data==1

    plot(ts_ecgi(offset:end),V1e_f(offset:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
ylim([min(real(ecg(:,3)))-5   max(real(ecg(:,3)))+5])

vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)
text(85,71,'Lead V1','Fontsize',16)
ylim([min(real(ecg(:,3)))-5   max(real(ecg(:,3)))+5])


grid on

i=8;
subplot('Position',[(.05+mod(i-1,3))/3 1.279-(ceil(i/3))/3 1/4 1/5])
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
if lge_very_strong==1
    plot(ts(1:size(Iss,1)),V2wss,'Color','m','LineWidth',3)
end
if data==1

     plot(ts_ecgi(offset:end),V2e_f(offset:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end    
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)
ylim([min(real(ecg(:,3)))-5   max(real(ecg(:,3)))+5])

text(85,70,'Lead V2','Fontsize',16)


grid on

i=9;
subplot('Position',[(0.05+mod(i-1,3))/3 1.279-(ceil(i/3))/3 1/4 1/5])
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
if lge_very_strong==1
    plot(ts(1:size(Iss,1)),V3wss,'Color','m','LineWidth',3)
end
if data==1

    plot(ts_ecgi(offset:end),V3e_f(offset:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)
ylim([min(real(ecg(:,3)))-5   max(real(ecg(:,3)))+5])

text(85,70,'Lead V3','Fontsize',16)


grid on

i=10;
subplot('Position',[(.15+mod(i-1,3))/3 1.36-(ceil(i/3))/3 1/4 1/5])
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
if lge_very_strong==1
    plot(ts(1:size(Iss,1)),V4wss,'Color','m','LineWidth',3)
end
if data==1

    plot(ts_ecgi(offset:end),V4e_f(offset:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)
ylim([min(real(ecg(:,3)))-5   max(real(ecg(:,3)))+5])

text(85,70,'Lead V4','Fontsize',16)

grid on

i=11;
subplot('Position',[(.05+mod(i-1,3))/3 1.36-(ceil(i/3))/3 1/4 1/5])
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
if lge_very_strong==1
    plot(ts(1:size(Iss,1)),V5wss,'Color','m','LineWidth',3)
end
if data==1

    plot(ts_ecgi(offset:end),V5e_f(offset:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');

set(gca,'Fontsize',16)
text(85,70,'Lead V5','Fontsize',16)
ylim([min(real(ecg(:,3)))-5   max(real(ecg(:,3)))+5])

grid on

i=12;
subplot('Position',[(.05+mod(i-1,3))/3 1.36-(ceil(i/3))/3 1/4 1/5])
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
if lge_very_strong==1
    plot(ts(1:size(Iss,1)),V6wss,'Color','m','LineWidth',3)
end
if data==1

    plot(ts_ecgi(offset:end),V6e_f(offset:end),'Color',[0.55, 0.0 , 0.0],'LineWidth',3)
end
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
grid on
set(gca,'Fontsize',16)
text(85,70,'Lead V6','Fontsize',16)
ylim([min(real(ecg(:,3)))-5   max(real(ecg(:,3)))+5])

grid on
if lge==1 && lge2==1
    legend('sim no LGE','sim lge light', 'sim lge intermediate','sim lge strong','recorded ECG');
else
     legend('control sim', 'LGE light','recorded ECG');

end
    %legend( 'sim lge intermediate','sim lge strong','sim lge very strong','recorded ECG');




%only precordial leads below



for i=1:size(ecg,2)
    score(i,1)=QRSs(i,5)>110;
    score(i,2)=duration(i,5)>55;
    if i<6
        score(i,3)=R_peak(i+1)>R_peak(i);
    end
    score(i,4)=ST_el(i);
end
    

