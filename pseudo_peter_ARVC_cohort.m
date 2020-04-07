addpath(genpath('C:/Users/petnov/Dropbox/'));
enableVTK;

frechet={};
for k=1:3
if k==1
    [Ic,IIc,IIIc,aVRc,aVLc,aVFc,V1wc,V2wc,V3wc,V4wc,V5wc,V6wc,LAc] =my_ecg_from_sim('09','D:\ARVC meshing automatic\patients\patient09\results\TestControlECGPat09\',202);
    [I,II,III,aVR,aVL,aVF,V1w,V2w,V3w,V4w,V5w,V6w,LA] =my_ecg_from_sim('09','D:\ARVC meshing automatic\patients\patient09\results\TestHeart09Roots1111111GadLV\',210);
[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e,Ie_f,IIe_f,IIIe_f,aVRe_f,aVLe_f,aVFe_f,V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f,ts_ecgi]=ecgi_ECG_pat09('D:\Data_ARVC_raw_electrodes_vest\Data_ARVC\ARVC9\ARVC_Baseline.mat','D:\ChrisReconstructions\ARVC9\ARVC9_Baseline.mat','D:\ARVC meshing automatic\patients\patient09\mpp\ECG_ELECTRODES.mat',I,II,III,aVR,aVL,aVF,V1w,V2w,V3w,V4w,V5w,V6w);

elseif k==2
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
    [Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e,Ie_f,IIe_f,IIIe_f,aVRe_f,aVLe_f,aVFe_f,V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f,ts_ecgi]=ecgi_ECG_pat09('D:\Data_ARVC_raw_electrodes_vest\Data_ARVC\ARVC18\ARVC18_Baseline.mat','D:\ChrisReconstructions\ARVC18\ARVC18_Baseline.mat','D:\ARVC meshing automatic\patients\patient18\mpp\ECG_ELECTRODES.mat',I,II,III,aVR,aVL,aVF,V1w,V2w,V3w,V4w,V5w,V6w);

else
      [Ic,IIc,IIIc,aVRc,aVLc,aVFc,V1wc,V2wc,V3wc,V4wc,V5wc,V6wc,LAc] =my_ecg_from_sim('05','D:\ARVC meshing automatic\patients\patient05\results\TestHeart05FullBeatandActivation\',202);
    [I,II,III,aVR,aVL,aVF,V1w,V2w,V3w,V4w,V5w,V6w,LA] =my_ecg_from_sim('05','D:\ARVC meshing automatic\patients\patient05\results\TestHeart05Roots0111111GadLVRV\',210);
    [Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e,Ie_f,IIe_f,IIIe_f,aVRe_f,aVLe_f,aVFe_f,V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f,ts_ecgi]=ecgi_ECG_pat09('D:\Data_ARVC_raw_electrodes_vest\Data_ARVC\ARVC5\ARVC5_Baseline.mat','D:\ChrisReconstructions\ARVC5\ARVC5_Baseline.mat','D:\ARVC meshing automatic\patients\patient05\mpp\ECG_ELECTRODES.mat',I,II,III,aVR,aVL,aVF,V1w,V2w,V3w,V4w,V5w,V6w);

end
ampl=max(Ie_f);
Ic=Ic/max(Ic)*max(Ie_f);
IIc =IIc/max(IIc) * max(IIe_f);
IIIc =IIIc/max(IIIc) * max(IIIe_f);
aVRc =aVRc/max(-aVRc) * max(-aVRe_f);
aVLc = aVLc/max(aVLc) * max(aVLe_f);
aVFc = aVFc/max(aVFc) * max(aVFe_f);
V1wc= V1wc/max(-V1wc) * max(-V1e_f);
V2wc = V2wc/max(-V2wc) * max(-V2e_f);
V3wc = V3wc/max(-V3wc) * max(-V3e_f);
V4wc = V4wc/max(-V4wc) * max(-V4e_f);
V5wc = V5wc/max(V5wc) * max(V5e_f);
V6wc = V6wc/max(V6wc) * max(V6e_f);


ampl=max(Ie_f);
I=I/max(I)*max(Ie_f);
II =II/max(II) * max(IIe_f);
III =III/max(III) * max(IIIe_f);
aVR =aVR/max(-aVR) * max(-aVRe_f);
aVL = aVL/max(aVL) * max(aVLe_f);
aVF = aVF/max(aVF) * max(aVFe_f);
V1w= V1w/max(-V1w) * max(-V1e_f);
V2w = V2w/max(-V2w) * max(-V2e_f);
V3w = V3w/max(-V3w) * max(-V3e_f);
V4w = V4w/max(-V4w) * max(-V4e_f);
V5w = V5w/max(V5w) * max(V5e_f);
V6w = V6w/max(V6w) * max(V6e_f);
%[~,idx]=min(abs(V6w(i)*ones(size(V6e_f,1),1)-V6e_f));
    y=[Ie_f,IIe_f,IIIe_f,aVRe_f,aVLe_f,aVFe_f,V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f];
    fit=[];
        set(0,'RecursionLimit',1000)

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
     frechet{j,1}=fit;
    frechet{j,2}=fit_control;
    
        y=[Ie_f,IIe_f,IIIe_f,aVRe_f,aVLe_f,aVFe_f,V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f];
        j=1;
    for x=[Ic,IIc,IIIc,aVRc,aVLc,aVFc,V1wc,V2wc,V3wc,V4wc,V5wc,V6wc]

        fit_control(j)=DiscreteFrechetDist(x,y(:,j));
        j=j+1;
    end
   
    end
  %try a Discrete Frechet distance as a measure of goodness of fit test to see if the values in the
  %simulation really come from a function which looks like body surface
  %potentials. This fuction measures the minimum distance required for a
  %dog and a walker to traverse their separate paths when a dog is on the
  %leash.
frechet_lead_mean = (frechet{1,1}+frechet{2,1}+frechet{3,1})/3
frechet_lead_mean_con = (frechet{1,2}+frechet{2,2}+frechet{3,2})/3

for i=1:3
    frechet_mean_patient(i) = mean([frechet{i,1}]);
    frechet_mean_patient_con(i) = mean(frechet{i,2});
end

  figure()
  scatter(1:12,frechet_lead_mean,100,'filled')
  hold on
  scatter(1:12,frechet_lead_mean_con,100,'filled')
  set(gca,'XTick',[1,2,3,4,5,6,7,8,9,10,11,12])
  set(gca,'XTickLabel',{'I','II','III','aVR','aVL','aVF','V1','V2','V3','V4','V5','V6'});
  legend('LGE model', 'control model')
  set(gca,'Fontsize',21)
  ylabel('Discrete Frechet Distance')
%%

  figure()
  scatter(1:3,frechet_mean_patient,100,'filled')
  hold on
  scatter(1:3,frechet_mean_patient_con,100,'filled')
  set(gca,'XTick',[1,2,3])
  set(gca,'XTickLabel',{'01','02','03'});
  legend('LGE model', 'control model')
  set(gca,'Fontsize',21)
  ylabel('Discrete Frechet Distance')
xlim([0 4])