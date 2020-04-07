addpath(genpath('C:/Users/petnov/Dropbox/'));
enableVTK;

try
    load('D:\ARVC meshing automatic\patients\patient18\results\TestHeart18ControlECG\ecg.mat');
catch
    printf('no ecg file was found! Start computing the ecg...')
end

load('D:\ARVC meshing automatic\patients\patient18\results\TestHeart18ControlECG\ecg.mat');
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
load('D:\ARVC meshing automatic\patients\patient18\results\TestHeart18LVRVLGEBidomain\ecg.mat');
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

%%
electrodes = load('D:\ARVC meshing automatic\patients\patient18\mpp\ECG_ELECTRODES.mat');


E = [ electrodes.ELECTRODES.LA ; electrodes.ELECTRODES.RA ; electrodes.ELECTRODES.LL ; electrodes.ELECTRODES.RL ; electrodes.ELECTRODES.V1 ; electrodes.ELECTRODES.V2 ; electrodes.ELECTRODES.V3 ; electrodes.ELECTRODES.V4 ; electrodes.ELECTRODES.V5 ; electrodes.ELECTRODES.V6 ];
E=E/10; %chaste is in cm
%plotMESH( MakeMesh( H , H.face ) , 'r' ); hplotMESH( T , 'nf' );
%hplot3d( E , '*m' ,'LineWidth',4,'MarkerSize',10 );
%axis(objbounds)

[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e,Ie_f,IIe_f,IIIe_f,aVRe_f,aVLe_f,aVFe_f,V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f,ts_ecgi]=ecgi_ECG('D:\Data_ARVC_raw_electrodes_vest\Data_ARVC\ARVC18\ARVC18_Baseline.mat','D:\ChrisReconstructions\ARVC18\ARVC18_Baseline.mat','D:\ARVC meshing automatic\patients\patient18\mpp\ECG_ELECTRODES.mat',I,II,III,aVR,aVL,aVF,V1w,V2w,V3w,V4w,V5w,V6w);


ampl=max(Ie_f);
Ic=Ic/max(Ic)*max(Ie_f);
IIc =IIc/max(IIc) * max(IIe_f);
IIIc =IIIc/max(IIIc) * max(IIIe_f);
aVRc =aVRc/max(-aVRc) * max(-aVRe_f);
aVLc = aVLc/max(aVLc) * max(-aVLe_f);
aVFc = aVFc/max(aVFc) * max(aVFe_f);
V1wc= V1wc/max(V1wc) * max(V1e_f);
V2wc = V2wc/max(V2wc) * max(V2e_f);
V3wc = V3wc/max(V3wc) * max(V3e_f);
V4wc = V4wc/max(V4wc) * max(V4e_f);
V5wc = V5wc/max(V5wc) * max(V5e_f);
V6wc = V6wc/max(V6wc) * max(V6e_f);


ampl=max(Ie_f);
I=I/max(I)*max(Ie_f);
II =II/max(II) * max(IIe_f);
III =III/max(III) * max(IIIe_f);
aVR =aVR/max(-aVR) * max(-aVRe_f);
aVL = aVL/max(aVL) * max(-aVLe_f);
aVF = aVF/max(aVF) * max(aVFe_f);
V1w= V1w/max(V1w) * max(V1e_f);
V2w = V2w/max(V2w) * max(V2e_f);
V3w = V3w/max(V3w) * max(V3e_f);
V4w = V4w/max(V4w) * max(V4e_f);
V5w = V5w/max(V5w) * max(V5e_f);
V6w = V6w/max(V6w) * max(V6e_f);



ts=1:size(aVRc,1);
offset=1;%idx
ts_ecgi=ts_ecgi-ts_ecgi(offset);
figure()
grid on
i=1;
subplot('Position',[(.05+mod(i-1,3))/3 1.1-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),Ic,'LineWidth',3)
hold all
%plot(ts_,'--r','LineWidth',3)
plot(ts(1:size(I,1)),I,'g','LineWidth',3)
plot(ts_ecgi(offset:end),Ie_f(offset:end),'k','LineWidth',3)

title('Lead I')
grid on
i=2;
subplot('Position',[(.05+mod(i-1,3))/3 1.1-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),IIc,'LineWidth',3)
hold on
%plot(III_real,'--r','LineWidth',3)
plot(ts(1:size(I,1)),II,'g','LineWidth',3)
plot(ts_ecgi(offset:end),IIe_f(offset:end),'k','LineWidth',3)

title('Lead II')
grid on

i=3;
subplot('Position',[(.05+mod(i-1,3))/3 1.1-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),IIIc,'LineWidth',3)
hold on
%plot(LA{1,1},IIIc,'--r','LineWidth',3)
plot(ts(1:size(I,1)),III,'g','LineWidth',3)
plot(ts_ecgi(offset:end),IIIe_f(offset:end),'k','LineWidth',3)

title('Lead III')
grid on

i=4;
subplot('Position',[(.05+mod(i-1,3))/3 1.19-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),aVRc,'LineWidth',3,'LineWidth',3)
hold on
%plot(LA{1,1},aVRc,'--r','LineWidth',3)
plot(ts(1:size(I,1)),aVR,'g','LineWidth',3)
plot(ts_ecgi(offset:end),aVRe_f(offset:end),'k','LineWidth',3)

title('Lead aVR')
grid on

i=5;
subplot('Position',[(.05+mod(i-1,3))/3 1.19-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),aVLc,'LineWidth',3)
hold on
%plot(LA{1,1},aVLc,'--r','LineWidth',3)
plot(ts(1:size(I,1)),aVL,'g','LineWidth',3)
plot(ts_ecgi(offset:end),aVLe_f(offset:end),'k','LineWidth',3)

title('Lead aVL')
grid on

i=6;
subplot('Position',[(.05+mod(i-1,3))/3 1.19-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),aVFc,'LineWidth',3)
hold on
%plot(LA{1,1},aVFc,'--r','LineWidth',3)
plot(ts(1:size(I,1)),aVF,'g','LineWidth',3)
plot(ts_ecgi(offset:end),aVFe_f(offset:end),'k','LineWidth',3)


title('Lead aVF')
grid on

i=7;
subplot('Position',[(.05+mod(i-1,3))/3 1.27-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),V1wc,'LineWidth',3)
hold on
%plot(LA{1,1},V1wc,'--r','LineWidth',3)
plot(ts(1:size(I,1)),V1w,'g','LineWidth',3)
plot(ts_ecgi(offset:end),V1e_f(offset:end),'k','LineWidth',3)


title('Lead V1')
grid on

i=8;
subplot('Position',[(.05+mod(i-1,3))/3 1.27-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),V2wc,'LineWidth',3)
hold on
%plot(LA{1,1},V2wc,'--r','LineWidth',3)
plot(ts(1:size(I,1)),V2w,'g','LineWidth',3)
plot(ts_ecgi(offset:end),V2e_f(offset:end),'k','LineWidth',3)

title('Lead V2')
grid on

i=9;
subplot('Position',[(.05+mod(i-1,3))/3 1.27-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),V3wc,'LineWidth',3)
hold on
%plot(LA{1,1},V3wc,'--r','LineWidth',3)
plot(ts(1:size(I,1)),V3w,'g','LineWidth',3)
plot(ts_ecgi(offset:end),V3e_f(offset:end),'k','LineWidth',3)

title('Lead V3')
grid on

i=10;
subplot('Position',[(.05+mod(i-1,3))/3 1.36-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),V4wc,'LineWidth',3)
hold on
%plot(LA{1,1},V4wc,'--r','LineWidth',3)
plot(ts(1:size(I,1)),V4w,'g','LineWidth',3)
plot(ts_ecgi(offset:end),V4e_f(offset:end),'k','LineWidth',3)

title('Lead V4')
grid on

i=11;
subplot('Position',[(.05+mod(i-1,3))/3 1.36-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),V5wc,'LineWidth',3)
hold on
%plot(LA{1,1},V5wc,'--r','LineWidth',3)
plot(ts(1:size(I,1)),V5w,'g','LineWidth',3)
plot(ts_ecgi(offset:end),V5e_f(offset:end),'k','LineWidth',3)

title('Lead V5')
grid on

i=12;
subplot('Position',[(.05+mod(i-1,3))/3 1.36-(ceil(i/3))/3 1/4 1/5])
plot(ts(1:size(Ic,1)),V6wc,'LineWidth',3)
hold on
%plot(LA{1,1},V6wc,'--r','LineWidth',3)
plot(ts(1:size(I,1)),V6w,'g','LineWidth',3)
plot(ts_ecgi(offset:end),V6e_f(offset:end),'k','LineWidth',3)

title('Lead V6')
grid on

legend('control sim','sim + fibrosis', 'body surface potentials');

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
    for x=[Ic,IIc,IIIc,aVRc,aVLc,aVFc,V1wc,V2wc,V3wc,V4wc,V5wc,V6wc]

        fit_control(j)=DiscreteFrechetDist(x,y(:,j));
        j=j+1;
    end
    
  %try a Discrete Frechet distance as a measure of goodness of fit test to see if the values in the
  %simulation really come from a function which looks like body surface
  %potentials. This fuction measures the minimum distance required for a
  %dog and a walker to traverse their separate paths when a dog is on the
  %leash.


  figure()
  scatter(1:12,fit,100,'filled')
  hold on
  scatter(1:12,fit_control,100,'filled')
  set(gca,'XTick',[1,2,3,4,5,6,7,8,9,10,11,12])
  set(gca,'XTickLabel',{'I','II','III','aVR','aVL','aVF','V1','V2','V3','V4','V5','V6'});
  legend('LGE model', 'control model')
  set(gca,'Fontsize',21)
  ylabel('Discrete Frechet Distance')
%%
