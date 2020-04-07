addpath(genpath('C:/Users/petnov/Dropbox/'));
addpath(genpath('C:/Users/petnov/Dropbox/qrs/'));

try
    load('D:\ARVC meshing automatic\patients\patient18\results\TestHeart18ControlECG\ecg.mat');
catch
    printf('no ecg file was found! Start computing the ecg...')
end

[Icon,IIcon,IIIcon,aVRcon,aVLcon,aVFcon,V1wcon,V2wcon,V3wcon,V4wcon,V5wcon,V6wcon,LAcon] =my_ecg_from_sim('09','D:\ARVC meshing automatic\patients\patient09\results\TestControlECGPat09\',202);
load('D:\ARVC meshing automatic\patients\patient09\results\TestHeart09LGELight\ecg_old.mat');
Is=I_real;
IIs=II_real;
IIIs=III_real;
aVRs=aVR_real;
aVLs=aVL_real;
aVFs=aVF_real;
V1ws=V1_real;
V2ws=V2_real;
V3ws=V3_real;
V4ws=V4_real;
V5ws=V5_real;
V6ws=V6_real;


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


load('D:\ARVC meshing automatic\patients\patient09\results\TestHeart09LGELight\ecg_old.mat');
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
locs=load('D:\Data_ARVC_raw_electrodes_vest\Data_ARVC\ARVC_electrode_locations');
locs=locs.arvc9_electrodes/10;
E_left=E;
E_left(:,1)=E(:,1)-2.5;
idx=[];
for i=1:10
    idx(i)=knnsearch(locs,E(i,:));
end

FBfine =transform(read_CHASTE(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\chaste',pat_nr,'\HEART')),'s',1/1);
H.face      = MeshBoundary( FBfine.tri );
H.triCENTER = meshFacesCenter( FBfine );
H.xyz=FBfine.xyz;

T= read_VTK(strcat('D:\ARVC meshing automatic\patients\patient','09','\mpp\','VEST0.vtk'));
T.xyz=T.xyz/10;
% 
figure()
plotMESH( MakeMesh( H , H.face ) , 'r' ); hplotMESH( T , 'nf' );
hplot3d( E , '*m' ,'LineWidth',4,'MarkerSize',10 );
hplot3d( locs(idx,:) , '*b' ,'LineWidth',4,'MarkerSize',10 );
axis(objbounds)
% 


[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e,Ie_f,IIe_f,IIIe_f,aVRe_f,aVLe_f,aVFe_f,V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f,ts_ecgi]=ecgi_ECG_pat09('D:\Data_ARVC_raw_electrodes_vest\Data_ARVC\ARVC9\ARVC_Baseline.mat','D:\ChrisReconstructions\ARVC9\ARVC9_Baseline.mat','D:\ARVC meshing automatic\patients\patient09\mpp\ECG_ELECTRODES.mat',9,I,II,III,aVR,aVL,aVF,V1w,V2w,V3w,V4w,V5w,V6w);

ampl=max(Ie_f);
Is=Is/max(Is)*max(Ie_f);
IIs =IIs/max(IIs) * max(IIe_f);
IIIs =IIIs/max(IIIs) * max(IIIe_f);
aVRs =aVRs/max(-aVRs) * max(-aVRe_f);
aVLs = aVLs/max(aVLs) * max(aVLe_f);
aVFs = aVFs/max(aVFs) * max(aVFe_f);
V1ws= V1ws/max(-V1ws) * max(-V1e_f);
V2ws = V2ws/max(-V2ws) * max(-V2e_f);
V3ws = V3ws/max(-V3ws) * max(-V3e_f);
V4ws = V4ws/max(V4ws) * max(V4e_f);
V5ws = V5ws/max(V5ws) * max(V5e_f);
V6ws = V6ws/max(V6ws) * max(V6e_f);

Icon=Icon/max(Icon)*max(Ie_f);
IIcon =IIcon/max(IIcon) * max(IIe_f);
IIIcon =IIIcon/max(IIIcon) * max(IIIe_f);
aVRcon =aVRcon/max(-aVRcon) * max(-aVRe_f);
aVLcon = aVLcon/max(aVLcon) * max(aVLe_f);
aVFcon = aVFcon/max(aVFcon) * max(aVFe_f);
V1wcon= V1wcon/max(-V1wcon) * max(-V1e_f);
V2wcon = V2wcon/max(-V2wcon) * max(-V2e_f);
V3wcon = V3wcon/max(-V3wcon) * max(-V3e_f);
V4wcon = V4wcon/max(V4wcon) * max(V4e_f);
V5wcon = V5wcon/max(V5wcon) * max(V5e_f);
V6wcon = V6wcon/max(V6wcon) * max(V6e_f);



ampl=max(Ie_f);
I=I/max(I)*max(Ie_f);
II =II/max(II) * max(IIe_f);
III =III/max(III) * max(IIIe_f);
aVR =aVR/max(-aVR) * max(-aVRe_f);
aVL = aVL/max(aVL) * max(aVLe_f);
aVF = aVF/max(aVF) * max(aVFe_f);
V1w= V1w/max(-V1w) * max((-V1e_f));
V2w = V2w/max(-V2w) * max((-V2e_f));
V3w = V3w/max(-V3w) * max((-V3e_f));
V4w = V4w/max(V4w) * max(V4e_f);
V5w = V5w/max(V5w) * max(V5e_f);
V6w = V6w/max(V6w) * max(V6e_f);


ampl=max(Ie_f);
% Iint=Iint/max(Iint)*max(Ie_f);
% IIint =IIint/max(IIint) * max(IIe_f);
% IIIint =IIIint/max(IIIint) * max(IIIe_f);
% aVRint =aVRint/max(-aVRint) * max(-aVRe_f);
% aVLint = aVLint/max(aVLint) * max(aVLe_f);
% aVFint = aVFint/max(aVFint) * max(aVFe_f);
% V1wint= V1wint/max(-V1wint) * max(-(V1e_f));
% V2wint = V2wint/max(-V2wint) * max(-(V2e_f));
% V3wint = V3wint/max(-V3wint) * max(-(V3e_f));
% V4wint = V4wint/max(V4wint) * max(V4e_f);
% V5wint = V5wint/max(V5wint) * max(V5e_f);
% V6wint = V6wint/max(V6wint) * max(V6e_f);


%extract QRS biomarkers

%medium = horzcat(V1wint,V2wint,V3wint,V4wint,V5wint,V6wint);
strong = horzcat(V1ws,V2ws,V3ws,V4ws,V5ws,V6ws);
light=horzcat(V1w,V2w,V3w,V4w,V5w,V6w);
nolge=horzcat(V1wcon,V2wcon,V3wcon,V4wcon,V5wcon,V6wcon);

ecg=horzcat(V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f);
dt_ecgi=ts_ecgi(8)-ts_ecgi(7);
for i=1:size(medium,2)
    QRSs(i,1)=qrs_estimator(light(:,i));
  %  QRSs(i,2)=qrs_estimator(medium(:,i));
    QRSs(i,3)=qrs_estimator(strong(:,i));
        QRSs(i,4)=qrs_estimator(nolge(:,i));

        QRSs(i,5)=qrs_estimator_bsp(ecg(:,i))*dt_ecgi;

end

%no significant difference in QRS durations
%now look at S wave area/peak


Swave=zeros(size(medium,2),10);
duration =zeros(size(medium,2),4);
for i=1:3
    [Swave(i,1),Swave(i,2),duration(i,1)]=s_wave_area(light(:,i));
  %[Swave(i,3),Swave(i,4),duration(i,2)]=s_wave_area(medium(:,i));
 [ Swave(i,5),Swave(i,6),duration(i,3)]=s_wave_area(strong(:,i));
  [ Swave(i,7),Swave(i,8),duration(i,4)]=s_wave_area(nolge(:,i));

       [Swave(i,9),Swave(i,10),duration(i,5)]=s_wave_area_bsp(ecg(:,i),ts_ecgi);

end

ratio=zeros(size(medium,2),5);
for i=1:size(medium,2)
    [~,ratio(i,1)]=qrs_estimator(light(:,i));
 %  [~, ratio(i,2)]=qrs_estimator(medium(:,i));
    [~,ratio(i,3)]=qrs_estimator(strong(:,i));
        [~,ratio(i,4)]=qrs_estimator(nolge(:,i));

       [~, ratio(i,5)]=qrs_estimator_bsp(ecg(:,i));

end


figure()
  bar(1:3,duration(1:3,:)./QRSs(1:3,:),'grouped')
  hold on
  set(gca,'XTick',[1,2,3])
  set(gca,'XTickLabel',{'V1','V2','V3'});
  legend('light LGE','medium LGE','strong LGE', 'no lge','recordings')
  set(gca,'Fontsize',21)
  ylabel('S wave as fraction of QRS Duration ')

  

figure()
  bar(1:3,Swave(1:3,[2 4 6 8 10 ])./duration(1:3,:),'grouped')
  hold on
  set(gca,'XTick',[1,2,3])
  set(gca,'XTickLabel',{'V1','V2','V3'});
  legend('light LGE','medium LGE','strong LGE', 'no lge','recordings')
  set(gca,'Fontsize',21)
  ylabel('S wave area/(s wave duration x s wave amplitude) ')



figure()
  bar(1:6,QRSs,'grouped')
  hold on
  set(gca,'XTick',[1,2,3,4,5,6])
  set(gca,'XTickLabel',{'V1','V2','V3','V4','V5','V6'});
  legend('light LGE','medium LGE','strong LGE', 'no lge','recordings')
  set(gca,'Fontsize',21)
  ylabel('QRS Duration (ms)')

  
  
figure()
  bar(1:3,horzcat(Swave(1:3,1),Swave(1:3,3),Swave(1:3,5),Swave(1:3,7)),'grouped')
  hold on
  set(gca,'XTick',[1,2,3,4])
  set(gca,'XTickLabel',{'V1','V2','V3'});
  legend('light LGE','medium LGE','strong LGE', 'recordings')
  set(gca,'Fontsize',21)
  ylabel('Area under S-wave (Vms)')
  
  
  figure()
  bar(1:3,horzcat(Swave(1:3,2),Swave(1:3,4),Swave(1:3,6),Swave(1:3,8)),'grouped')
  hold on
  set(gca,'XTick',[1,2,3,4])
  set(gca,'XTickLabel',{'V1','V2','V3'});
  legend('light LGE','medium LGE','strong LGE', 'recordings')
  set(gca,'Fontsize',21)
  ylabel('Area under S-wave/Peak S-wave (ms)')
  
    figure()
  bar(1:3,horzcat(ratio(1:3,1),ratio(1:3,2),ratio(1:3,3),ratio(1:3,4)),'grouped')
  hold on
  set(gca,'XTick',[1,2,3,4])
  set(gca,'XTickLabel',{'V1','V2','V3'});
  legend('light LGE','medium LGE','strong LGE', 'recordings')
  set(gca,'Fontsize',21)
  ylabel('R peak to S peak ratio')
  
  
  

  


lge=1;
lge2=0;
control=0;
ts=1:208;
offset=1;%idx
ts_ecgi=ts_ecgi-ts_ecgi(offset);
figure()
grid on
i=1;
subplot('Position',[(.15+mod(i-1,3))/3 1.125-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Is,1)),Is,'LineWidth',3)
hold all
%plot(ts_,'--r','LineWidth',3)
if lge==1
    plot(ts(1:size(I,1)),I,'g','LineWidth',3)
end
if lge2==1
    plot(ts(1:size(Iint,1)),Iint,'r','LineWidth',3)
end
if control==1
        plot(ts(1:size(Icon,1)),Icon,'m','LineWidth',3)
end
plot(ts_ecgi(offset:end),Ie_f(offset:end),'k','LineWidth',3)
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
text(200,25,'Lead I','Fontsize',16)
grid on
set(gca,'Fontsize',16)

i=2;
subplot('Position',[(0.05+mod(i-1,3))/3 1.125-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Is,1)),IIs,'LineWidth',3)
hold on
%plot(III_real,'--r','LineWidth',3)
if lge==1
    plot(ts(1:size(I,1)),II,'g','LineWidth',3)
end
if lge2==1
    plot(ts(1:size(IIint,1)),IIint,'r','LineWidth',3)
end
if control==1
        plot(ts(1:size(IIcon,1)),IIcon,'m','LineWidth',3)
end
plot(ts_ecgi(offset:end),IIe_f(offset:end),'k','LineWidth',3)
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
text(200,25,'Lead II','Fontsize',16)
grid on
set(gca,'Fontsize',16)

i=3;
subplot('Position',[(.05+mod(i-1,3))/3 1.125-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Is,1)),IIIs,'LineWidth',3)
hold on
%plot(LA{1,1},IIIc,'--r','LineWidth',3)
if lge==1
    plot(ts(1:size(I,1)),III,'g','LineWidth',3)
end
if lge2==1
    plot(ts(1:size(IIIint,1)),IIIint,'r','LineWidth',3)
end
if control==1
        plot(ts(1:size(IIIcon,1)),IIIcon,'m','LineWidth',3)
end
plot(ts_ecgi(offset:end),IIIe_f(offset:end),'k','LineWidth',3)
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)
text(200,7,'Lead III','Fontsize',16)
i=4;
subplot('Position',[(.15+mod(i-1,3))/3 1.195-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Is,1)),aVRs,'LineWidth',3,'LineWidth',3)
hold on
%plot(LA{1,1},aVRc,'--r','LineWidth',3)
if lge==1
    plot(ts(1:size(I,1)),aVR,'g','LineWidth',3)
end
if lge2==1
    plot(ts(1:size(Iint,1)),aVRint,'r','LineWidth',3)
end
if control==1
        plot(ts(1:size(aVRcon,1)),aVRcon,'m','LineWidth',3)
end
plot(ts_ecgi(offset:end),aVRe_f(offset:end),'k','LineWidth',3)
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
text(200,5,'Lead aVR','Fontsize',16)
grid on
set(gca,'Fontsize',16)
i=5;
subplot('Position',[(.05+mod(i-1,3))/3 1.195-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Is,1)),aVLs,'LineWidth',3)
hold on
%plot(LA{1,1},aVLc,'--r','LineWidth',3)
if lge==1
    plot(ts(1:size(I,1)),aVL,'g','LineWidth',3)
end
if lge2==1
    plot(ts(1:size(Iint,1)),aVLint,'r','LineWidth',3)
end
if control==1
        plot(ts(1:size(aVLcon,1)),aVLcon,'m','LineWidth',3)
end
plot(ts_ecgi(offset:end),aVLe_f(offset:end),'k','LineWidth',3)
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
text(200,10,'Lead aVL','Fontsize',16)


grid on
set(gca,'Fontsize',16)
i=6;
subplot('Position',[(.05+mod(i-1,3))/3 1.195-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Is,1)),aVFs,'LineWidth',3)
hold on
%plot(LA{1,1},aVFc,'--r','LineWidth',3)
if lge==1

    plot(ts(1:size(I,1)),aVF,'g','LineWidth',3)
end
if lge2==1
    plot(ts(1:size(Iint,1)),aVFint,'r','LineWidth',3)
end
if control==1
        plot(ts(1:size(aVFcon,1)),aVFcon,'m','LineWidth',3)
end
plot(ts_ecgi(offset:end),aVFe_f(offset:end),'k','LineWidth',3)

xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)
text(200,15,'Lead aVF','Fontsize',16)


grid on

i=7;
subplot('Position',[(.15+mod(i-1,3))/3 1.279-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Is,1)),V1ws,'LineWidth',3)
hold on
%plot(LA{1,1},V1wc,'--r','LineWidth',3)
if lge==1

    plot(ts(1:size(I,1)),V1w,'g','LineWidth',3)
end
if lge2==1
    plot(ts(1:size(Iint,1)),V1wint,'r','LineWidth',3)
end
if control==1
        plot(ts(1:size(V1wcon,1)),V1wcon,'m','LineWidth',3)
end
plot(ts_ecgi(offset:end),V1e_f(offset:end),'k','LineWidth',3)

xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)
text(200,35,'Lead V1','Fontsize',16)


grid on

i=8;
subplot('Position',[(.05+mod(i-1,3))/3 1.279-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Is,1)),V2ws,'LineWidth',3)
hold on
%plot(LA{1,1},V2wc,'--r','LineWidth',3)
if lge==1

    plot(ts(1:size(I,1)),V2w,'g','LineWidth',3)
end
if lge2==1
    plot(ts(1:size(Iint,1)),V2wint,'r','LineWidth',3)
end
if control==1
        plot(ts(1:size(V2wcon,1)),V2wcon,'m','LineWidth',3)
end
 plot(ts_ecgi(offset:end),V2e_f(offset:end),'k','LineWidth',3)
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)

text(200,30,'Lead V2','Fontsize',16)


grid on

i=9;
subplot('Position',[(0.05+mod(i-1,3))/3 1.279-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Is,1)),V3ws,'LineWidth',3)
hold on
%plot(LA{1,1},V3wc,'--r','LineWidth',3)
if lge==1
    plot(ts(1:size(I,1)),V3w,'g','LineWidth',3)
end
if lge2==1
    plot(ts(1:size(Iint,1)),V3wint,'r','LineWidth',3)
end
if control==1
        plot(ts(1:size(V3wcon,1)),V3wcon,'m','LineWidth',3)
end
plot(ts_ecgi(offset:end),V3e_f(offset:end),'k','LineWidth',3)
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)

text(200,35,'Lead V3','Fontsize',16)


grid on

i=10;
subplot('Position',[(.15+mod(i-1,3))/3 1.36-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Is,1)),V4ws,'LineWidth',3)
hold on
%plot(LA{1,1},V4wc,'--r','LineWidth',3)
if lge==1

    plot(ts(1:size(I,1)),V4w,'g','LineWidth',3)
end
if lge2==1
    plot(ts(1:size(Iint,1)),V4wint,'r','LineWidth',3)
end
if control==1
        plot(ts(1:size(V4wcon,1)),V4wcon,'m','LineWidth',3)
end
plot(ts_ecgi(offset:end),V4e_f(offset:end),'k','LineWidth',3)
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');
grid on
set(gca,'Fontsize',16)

text(200,250,'Lead V4','Fontsize',16)

grid on

i=11;
subplot('Position',[(.05+mod(i-1,3))/3 1.36-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Is,1)),V5ws,'LineWidth',3)
hold on
%plot(LA{1,1},V5wc,'--r','LineWidth',3)
if lge==1

    plot(ts(1:size(I,1)),V5w,'g','LineWidth',3)
end
if lge2==1
    plot(ts(1:size(Iint,1)),V5wint,'r','LineWidth',3)
end
if control==1
        plot(ts(1:size(V5wcon,1)),V5wcon,'m','LineWidth',3)
end
plot(ts_ecgi(offset:end),V5e_f(offset:end),'k','LineWidth',3)
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
vec_pos=get(get(gca,'XLabel'),'Position');

set(gca,'Fontsize',16)
text(200,80,'Lead V5','Fontsize',16)

grid on

i=12;
subplot('Position',[(.05+mod(i-1,3))/3 1.36-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Is,1)),V6ws,'LineWidth',3)
hold on
%plot(LA{1,1},V6wc,'--r','LineWidth',3)
if lge==1

    plot(ts(1:size(I,1)),V6w,'g','LineWidth',3)
end
if lge2==1
    plot(ts(1:size(Iint,1)),V6wint,'r','LineWidth',3)
end
if control==1
        plot(ts(1:size(V6wcon,1)),V6wcon,'m','LineWidth',3)
end
    plot(ts_ecgi(offset:end),V6e_f(offset:end),'k','LineWidth',3)
xlabel('ms','Fontsize',16)
ylabel('mV','Fontsize',16)
grid on
set(gca,'Fontsize',16)
text(200,10,'Lead V6','Fontsize',16)

grid on
if lge==1 && lge2==1
    legend('sim lge statistical torso positions','sim lge ecgi vest locations','recorded ECG');
elseif lge==1
        
    legend('sim lge statistical torso positions','sim lge ecgi vest locations','recorded ECG');

end
%Compute a goodness of fit to the ECGi data


%[~,idx]=min(abs(V6w(i)*ones(size(V6e_f,1),1)-V6e_f));
    y=[Ie_f,IIe_f,IIIe_f,aVRe_f,aVLe_f,aVFe_f,V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f];
    fit=[];
    j=1;
    for x=[I,II,III,aVR,aVL,aVF,V1w,V2w,V3w,V4w,V5w,V6w]
        sum_squared=0;
        for i=1:size(x,1)
            sum_squared = power(real(x(i))-real(y(i,j)),2)+sum_squared;
        end
        j
        fit(j)=DiscreteFrechetDist(x,y(:,j));
        j=j+1;
    end
    
        y=[Ie_f,IIe_f,IIIe_f,aVRe_f,aVLe_f,aVFe_f,V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f];
        j=1;
    for x=[Is,IIs,IIIs,aVRs,aVLs,aVFs,V1ws,V2ws,V3ws,V4ws,V5ws,V6ws]

        fit_control(j)=DiscreteFrechetDist(x,y(:,j));
        j=j+1;
    end
    
  %try a Discrete Frechet distance as a measure of goodness of fit test to see if the values in the
  %simulation really come from a function which looks like body surface
  %potentials. This fuction measures the minimum distance required for a
  %dog and a walker to traverse their separate paths when a dog is on the
  %leash.
  
  
  %quantify the biomarkers on the QRS
  
  
  
  
  
  
a=real(fit./fit_control);
  figure()
  bar(1:3,vertcat([a(1),a(7),a(8)], ones(1,3))','grouped')
  hold on
  set(gca,'XTick',[1,2,3,4,5,6,7,8,9,10,11,12])
  set(gca,'XTickLabel',{'I','V1','V2'});
  legend('LGE model', 'control model')
  set(gca,'Fontsize',21)
  ylabel('Normalised Frechet Distance')
  xlim([0 13])
  
    figure()
  bar(1:12,vertcat(real(fit(1./fit_control), ones(1,12))','grouped'))
  hold on
  set(gca,'XTick',[1,2,3,4,5,6,7,8,9,10,11,12])
  set(gca,'XTickLabel',{'I','II','III','aVR','aVL','aVF','V1','V2','V3','V4','V5','V6'});
  legend('LGE model', 'control model')
  set(gca,'Fontsize',21)
  ylabel('Normalised Frechet Distance')
  xlim([0 13])
%%
  %make a movie of lead V4 in order to investigate epsilon wave

% 
% for i=1:71
%     h = figure();
% 
%     hold on
%        plot(ts(1:size(I,1)),V4w,'g','LineWidth',3)
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
