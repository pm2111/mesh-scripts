addpath(genpath('/data/DropPeter/Dropbox/qrs/'));
addpath(genpath('/data/DropPeter/Dropbox/tools/'));
addpath(genpath('/data/DropPeter/Dropbox/mesh scripts/'));
addpath(genpath('/data/DropPeter/Dropbox/sim3d scripts/'));


close all
                [ecg,ts_ecgi]=ecgi_ECG06(['/data/Data_ARVC_raw_electrodes_vest/Data_ARVC/ARVC' num2str(j) '/ARVC' num2str(j) '_Baseline.mat'],['/data/ChrisReconstructions/ARVC' num2str(j) '/ARVC' num2str(j) '_Baseline.mat'],['/run/media/petnov/Seagate Expansion Drive/peter data backup windows PC/ARVC meshing automatic/patients/patient' pat_nr '/mpp/ECG_ELECTRODES.mat'],6);
ecg_ecgi=[V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f];
sim_path='/run/media/petnov/Seagate Expansion Drive/peter data backup windows PC/patients/patient06/results/Eikonal/exp30/';
plt_acts=0;
keep_models={};
patient_nr='09';
sim_name='Eikonal';
study='exp32';
patient_sim_folder =  strcat('/run/media/petnov/Seagate Expansion Drive/peter data backup windows PC/patients/patient',patient_nr,'/');
%H_mesh = read_CHASTE( strcat(patient_sim_folder,    'chaste',patient_nr,'\HEART'));
for z=[1]
    if z==1
        pop_path='/run/media/petnov/Seagate Expansion Drive/peter data backup windows PC/ARVC meshing automatic/patients/patient06/results/Eikonal_coarse/eikonal06_coarse_large_population4_ECGs/';
    else
        pop_path='/run/media/petnov/Seagate Expansion Drive/peter data backup windows PC/ARVC meshing automatic/patients/patient09/results/Eikonal/exp7_ECGs/';
    end
names = dir(pop_path);
duration=[];
QRSs=[];
ratio=[];
R_peak=[];
nb=size(names,1);
ecgs={};
upstroke={};
% ecg1=load(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',sim_name,'\',sim_type,'.mat'));
% ecg1=ecg1.ecg_full;
sim_path='/run/media/petnov/Seagate Expansion Drive/peter data backup windows PC/patients/patient06/results/Eikonal/exp_roots_light_pheno/';

folders_acts=dir(sim_path);
for i = [1]%3:nb

    names_acts=dir(strcat(sim_path,folders_acts(i+2).name));
    upstroke{i}=load(strcat(sim_path,folders_acts(i+2).name,'/',names_acts(3).name));
    
end



    for j=1:6
      [QRSs_ecgi(j)]=abs(qrs_estimator_bsp(ecg_ecgi(:,j))*(ts_ecgi(2)-ts_ecgi(1)));
        [~,ratio_ecgi(j),~,~,~,~,~,~]=qrs_estimator_bsp(ecg_ecgi(:,j));
        [ ~,~,duration_ecgi(j),~,RS_grad(j)]=s_wave_area_bsp(ecg_ecgi(:,j),1:size(ecg_ecgi,1)); 
    end



labels= {  'no LGE and all roots present'   'LBBB control' 'no LGE global conduction scaling 0.7 + LBBB'    'LGE regional conduction scaling + LBBB'    };



k=1;
model_names ={};
%for o=2:nb/2
figure('Position',[0 0 250 500]); %left bottom width height
z=1;
for i =[3]%o+(1:2)%o+122*(0:1)%o+(1:4)%[-5+o+ 5*(1:5)]%[2+o+36*(0:3)]% [-3+o+ 6*(1:6)]%[-5+o+ 5*(1:5)] %[o+(0:3)*25]%[-2+ 5*(1:5)] %3:nb
    idx=[];
    ecg=load([pop_path names(i).name]);
     model_names{z}=names(i).name;

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
      % [Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e]=ecgi_ECG_sim_plot_mat(ecg1,V3e_f,0.1);
    ecgs{i-2}=[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e];
    full_ecg={};
        full_ecg={Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e};

        plot_ecgs_2x6_tight(full_ecg',names(3).name,0,i-2,1:size(V1e,1));
               % plot_ecgs_2x6_clinical(full_ecg',names(3).name,0,z);

            hold on

%    
%     
%     else
%         dummy_name=model_names{k-1};
%         dummy_name(7)='2';
%         A=strfind({names.name},dummy_name);
%         idx=find(not(cellfun('isempty',A)));
z=z+1;
end


word_dir='/data/Data_ARVC_raw_electrodes_vest/ecgs/';
scores={};
s_waves=[];
sts=[];
QRS_pop=[];
grads=[];
PRs=[];
for j=[1:20]
    if j<10
      pat_nr=[ '0' num2str(j)];
    else
       pat_nr=num2str(j);
    end
% electrodes = load(['D:\ARVC meshing automatic\patients\patient' pat_nr '\mpp\ECG_ELECTRODES.mat']);
% 
% 
% E = [ electrodes.ELECTRODES.LA ; electrodes.ELECTRODES.RA ; electrodes.ELECTRODES.LL ; electrodes.ELECTRODES.RL ; electrodes.ELECTRODES.V1 ; electrodes.ELECTRODES.V2 ; electrodes.ELECTRODES.V3 ; electrodes.ELECTRODES.V4 ; electrodes.ELECTRODES.V5 ; electrodes.ELECTRODES.V6 ];
% E=E/10; %chaste is in cm
% %plotMESH( MakeMesh( H , H.face ) , 'g' ); hplotMESH( T , 'nf' );
% %hplot3d( E , '*m' ,'LineWidth',4,'MarkerSize',10 );
% %axis(objbounds)
if   j==19  ||j==12 || j==6 
    try 
                [ecg,ts_ecgi]=ecgi_ECG06(['/data/Data_ARVC_raw_electrodes_vest/Data_ARVC/ARVC' num2str(j) '/ARVC' num2str(j) '_Baseline.mat'],['/data/ChrisReconstructions/ARVC' num2str(j) '/ARVC' num2str(j) '_Baseline.mat'],['/run/media/petnov/Seagate Expansion Drive/peter data backup windows PC/ARVC meshing automatic/patients/patient' pat_nr '/mpp/ECG_ELECTRODES.mat'],j);
    catch
        [ecg,ts_ecgi]=ecgi_ECG06(['/data/Data_ARVC_raw_electrodes_vest/Data_ARVC/ARVC' num2str(j) '/ARVC' num2str(j) '_Baseline2.mat'],['/data/ChrisReconstructions/ARVC' num2str(j) '/ARVC' num2str(j) '_Baseline.mat'],['/run/media/petnov/Seagate Expansion Drive/peter data backup windows PC/ARVC meshing automatic/patients/patient' pat_nr '/mpp/ECG_ELECTRODES.mat'],j);
    end
elseif  j==20|| j==8 || j==17 ||j==13 || j==11
            [ecg,ts_ecgi]=ecgi_ECG06(['/data/Data_ARVC_raw_electrodes_vest/Data_ARVC/ARVC' num2str(j) '/ARVC' num2str(j) '_Baseline2.mat'],['/data/ChrisReconstructions/ARVC' num2str(j) '/ARVC' num2str(j) '_Baseline.mat'],['/run/media/petnov/Seagate Expansion Drive/peter data backup windows PC/ARVC meshing automatic/patients/patient' pat_nr '/mpp/ECG_ELECTRODES.mat'],j);
elseif j==15 
        [ecg,ts_ecgi]=ecgi_ECG06(['/data/Data_ARVC_raw_electrodes_vest/Data_ARVC/ARVC' num2str(j) '/ARVC' num2str(j) '_Baseline3.mat'],['/data/ChrisReconstructions/ARVC' num2str(j) '/ARVC' num2str(j) '_Baseline.mat'],['/run/media/petnov/Seagate Expansion Drive/peter data backup windows PC/ARVC meshing automatic/patients/patient' pat_nr '/mpp/ECG_ELECTRODES.mat'],j);
elseif j==9 ||  j==14 ||j==18 ||j==19|| j==2 || j==4 || j==7 || j==19 ||j==5
    
           [ecg,ts_ecgi]=ecgi_ECG_pat09(['/data/Data_ARVC_raw_electrodes_vest/Data_ARVC/ARVC' num2str(j) '/ARVC' num2str(j) '_Baseline.mat'],['/data/ChrisReconstructions/ARVC' num2str(j) '/ARVC' num2str(j) '_Baseline.mat'],['/run/media/petnov/Seagate Expansion Drive/peter data backup windows PC/ARVC meshing automatic/patients/patient' pat_nr '/mpp/ECG_ELECTRODES.mat'],j);
    
elseif j==10 ||j==4  
          [ecg,ts_ecgi]=ecgi_ECG_pat09(['/data/Data_ARVC_raw_electrodes_vest/Data_ARVC/ARVC' num2str(j) '/ARVC' num2str(j) '_Baseline2.mat'],['/data/ChrisReconstructions/ARVC' num2str(j) '/ARVC' num2str(j) '_Baseline.mat'],['/run/media/petnov/Seagate Expansion Drive/peter data backup windows PC/ARVC meshing automatic/patients/patient' pat_nr '/mpp/ECG_ELECTRODES.mat'],j);

else
        [ecg,ts_ecgi]=ecgi_ECG06(['/data/Data_ARVC_raw_electrodes_vest/Data_ARVC/ARVC' num2str(j) '/ARVC' num2str(j) '_Baseline.mat'],['/data/ChrisReconstructions/ARVC' num2str(j) '/ARVC' num2str(j) '_Baseline.mat'],['/run/media/petnov/Seagate Expansion Drive/peter data backup windows PC/ARVC meshing automatic/patients/patient' pat_nr '/mpp/ECG_ELECTRODES.mat'],j);
    
end
   % ecg=horzcat(Ie_f,IIe_f,IIIe_f,aVRe_f,aVLe_f,aVFe_f,V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f);
    ECGs{j}=ecg;
dt_ecgi=ts_ecgi(8)-ts_ecgi(7);
R_peak=[];

end
r_idx=zeros(6,20);
s_idx=zeros(6,20);
s_nadir_idx=zeros(6,20);
for l=[1:20]
    ecg=ECGs{l};
    data=load(strcat('/data/ChrisReconstructions/ARVC',num2str(l),'/ARVC',num2str(l),'_Baseline.mat'));
    nm=fieldnames(data);
    ts_ecgi=data.(nm{1}).epipots.t;
    dt=1000/2048;
    ts_mod=[1:size(ecg,1)]*dt;
    dummy=ECGs{l};
  

 
    offset=1;
    for m=1:12
          ecg(:,m)=ecg(:,m)-ecg(offset,m)*ones(size(ecg,1),1);
    end
full_ecg={ecg(:,1),ecg(:,2),ecg(:,3),ecg(:,4),ecg(:,5),ecg(:,6),ecg(:,7),ecg(:,8),ecg(:,9),ecg(:,10),ecg(:,11),ecg(:,12)};%plot clinical ecgs

Fs=2048;
f_low=35;
f_high=0.5;

for z=1:12
    [D,C]=butter(7,f_low/(Fs/2),'low');
        [D1,C1]=butter(5,f_high/(Fs/2),'high');

    full_ecg_filtered_high{z}=baseline_ecg(filtfilt(D,C,full_ecg{z}),2048,0.5);
        %qrs_filtered_high{z}=baseline_ecg(filtfilt(D,C,full_ecg_ecgi{z}),2048,0.5);
   % full_ecg_filtered_high{z}=filtfilt(D1,C1,filtfilt(D,C,full_ecg{z}));

end
for j=1:6
      [QRS_ecgi(j,l),ratio_ecgi(j,l),S_wave(j,l),r_idx(j,l),s_idx(j,l),s_nadir_idx(j,l),duration_ecgi(j,l),r_peak_idx(j,l)]=qrs_estimator_bsp2(full_ecg_filtered_high{j+6});
end
 


%manual tuning of qrs params
if l==20
  
       s_idx(6,l)=970;
       s_idx(5,l)=991;
       s_idx(4,l)=963;
       s_idx(2,l)=933;
       s_idx(1,l)=991;

    duration_ecgi(:,l)=(s_idx(:,l)-s_nadir_idx(:,l))*1000/2048;
    QRS_ecgi(:,l)=abs(r_idx(:,l)-s_idx(:,l))*1000/2048;
end
if l==19
  
       s_idx(5,l)=450;
             r_idx(4,l)=245;

    duration_ecgi(:,l)=(s_idx(:,l)-s_nadir_idx(:,l))*1000/2048;
    QRS_ecgi(:,l)=abs(r_idx(:,l)-s_idx(:,l))*1000/2048;
end
if l==18
  
       s_idx(4,l)=543;
       r_idx(4,l)=345;
       s_idx(3,l)=527;
       r_idx(2,l)=319;
       s_idx(2,l)=528;
       s_idx(1,l)=550;

    duration_ecgi(:,l)=(s_idx(:,l)-s_nadir_idx(:,l))*1000/2048;
    QRS_ecgi(:,l)=abs(r_idx(:,l)-s_idx(:,l))*1000/2048;
end
if l==17
  
       s_idx(6,l)=571;
       r_idx(6,l)=388;
       r_idx(4,l)=367;
       s_idx(4,l)=578;


    duration_ecgi(:,l)=(s_idx(:,l)-s_nadir_idx(:,l))*1000/2048;
    QRS_ecgi(:,l)=abs(r_idx(:,l)-s_idx(:,l))*1000/2048;
end
if l==16
  
    r_idx(3,l)=587;
       r_idx(2,l)=587;
       s_idx(2,l)=832;


    duration_ecgi(:,l)=(s_idx(:,l)-s_nadir_idx(:,l))*1000/2048;
    QRS_ecgi(:,l)=abs(r_idx(:,l)-s_idx(:,l))*1000/2048;
end
if l==15
   s_idx(6,l)=552;
    r_idx(6,l)=314;
    r_idx(4,l)=263;
    r_idx(3,l)=247;
    r_idx(2,l)=257;
    r_idx(1,l)=267;

    duration_ecgi(:,l)=(s_idx(:,l)-s_nadir_idx(:,l))*1000/2048;
    QRS_ecgi(:,l)=abs(r_idx(:,l)-s_idx(:,l))*1000/2048;
end
if l==14
   s_idx(6,l)=573;
   s_idx(5,l)=545;
    r_idx(5,l)=334;

    duration_ecgi(:,l)=(s_idx(:,l)-s_nadir_idx(:,l))*1000/2048;
    QRS_ecgi(:,l)=abs(r_idx(:,l)-s_idx(:,l))*1000/2048;
end
if l==13
    r_idx(5,l)=445;

    duration_ecgi(:,l)=(s_idx(:,l)-s_nadir_idx(:,l))*1000/2048;
    QRS_ecgi(:,l)=abs(r_idx(:,l)-s_idx(:,l))*1000/2048;
end
if l==12
    r_idx(6,l)=70;
    r_idx(5,l)=65;
    s_idx(5,l)=258;
    r_idx(4,l)=43;
    r_idx(2,l)=49;
    r_idx(1,l)=57;

    duration_ecgi(:,l)=(s_idx(:,l)-s_nadir_idx(:,l))*1000/2048;
    QRS_ecgi(:,l)=abs(r_idx(:,l)-s_idx(:,l))*1000/2048;
end
if l==9
    r_idx(1,l)=307;
    r_idx(2,l)=290;
    r_idx(3,l)=299;
    r_idx(4,l)=313;
    r_idx(5,l)=313;
    r_idx(6,l)=310;
    s_idx(1,l)=490;
    s_idx(2,l)=495;
    s_idx(3,l)=507;
    s_idx(4,l)=518;

    duration_ecgi(:,l)=(s_idx(:,l)-s_nadir_idx(:,l))*1000/2048;
    QRS_ecgi(:,l)=abs(r_idx(:,l)-s_idx(:,l))*1000/2048;
end
if l==11
    r_idx(6,l)=218;
    r_idx(5,l)=172;
     s_idx(5,l)=401;
    r_idx(4,l)=209;
    r_idx(3,l)=203;
    r_idx(1,l)=184;

    duration_ecgi(:,l)=(s_idx(:,l)-s_nadir_idx(:,l))*1000/2048;
    QRS_ecgi(:,l)=abs(r_idx(:,l)-s_idx(:,l))*1000/2048;
end
if l==10
    r_idx(1,l)=255;
    s_idx(1,l)=493;
    s_idx(2,l)=491;
    s_idx(3,l)=516;
    s_idx(4,l)=511;
    r_idx(6,l)=293;

      duration_ecgi(:,l)=(s_idx(:,l)-s_nadir_idx(:,l))*1000/2048;
    QRS_ecgi(:,l)=abs(r_idx(:,l)-s_idx(:,l))*1000/2048;
end
if l==1
    r_idx(3,l)=341;

    duration_ecgi(:,l)=(s_idx(:,l)-s_nadir_idx(:,l))*1000/2048;
    QRS_ecgi(:,l)=abs(r_idx(:,l)-s_idx(:,l))*1000/2048;
end
if l==2
    s_idx(5,l)=536;
    r_idx(4,l)=325;
    r_idx(3,l)=323;
    s_idx(2,l)=561;
    r_idx(1,l)=320;
    duration_ecgi(:,l)=(s_idx(:,l)-s_nadir_idx(:,l))*1000/2048;
    QRS_ecgi(:,l)=abs(r_idx(:,l)-s_idx(:,l))*1000/2048;
end
if l==3
    r_idx(5,l)=185;
    s_idx(4,l)=422;
    r_idx(2,l)=208;
    r_idx(1,l)=217;

    duration_ecgi(:,l)=(s_idx(:,l)-s_nadir_idx(:,l))*1000/2048;
    QRS_ecgi(:,l)=abs(r_idx(:,l)-s_idx(:,l))*1000/2048;
end
if l==4
    s_idx(1,l)=332;
    r_idx(2,l)=96;
    s_idx(2,l)=370;
    r_idx(3,l)=106;
    r_idx(4,l)=114;
    r_idx(5,l)=122;
    r_idx(6,l)=133;
    duration_ecgi(:,l)=(s_idx(:,l)-s_nadir_idx(:,l))*1000/2048;
    QRS_ecgi(:,l)=abs(r_idx(:,l)-s_idx(:,l))*1000/2048;
end
if l==5
    r_idx(5,l)=353;
     r_idx(4,l)=296;
     s_idx(4,l)=578;
    r_idx(3,l)=300;
    s_idx(3,l)=561;
    s_idx(1,l)=495;
    duration_ecgi(:,l)=(s_idx(:,l)-s_nadir_idx(:,l))*1000/2048;
    QRS_ecgi(:,l)=abs(r_idx(:,l)-s_idx(:,l))*1000/2048;
end
if l==6
    r_idx(1,l)=274;
  
    duration_ecgi(:,l)=(s_idx(:,l)-s_nadir_idx(:,l))*1000/2048;
    QRS_ecgi(:,l)=abs(r_idx(:,l)-s_idx(:,l))*1000/2048;
end
if l==8
    s_idx(6,l)=481;
    r_idx(5,l)=281;
      s_idx(2,l)=438;
      s_idx(1,l)=430;

    duration_ecgi(:,l)=(s_idx(:,l)-s_nadir_idx(:,l))*1000/2048;
    QRS_ecgi(:,l)=abs(r_idx(:,l)-s_idx(:,l))*1000/2048;
end
for z=1:12

    dummy=full_ecg_filtered_high{z};
    qrs_filtered{z}=dummy(min(r_idx(:,l)-10):max(s_idx(:,l)));
   % qrs_filtered{z}(453)=0;
end

figure()
%plot_ecgs_2x6_bios(full_ecg_filtered_high,'1',0,5,ts_mod,[r_idx(:,l) ;r_idx(:,l)],[s_idx(:,l); s_idx(:,l) ],[ s_nadir_idx(:,l);s_nadir_idx(:,l) ]);
plot_ecgs_2x6_tight(qrs_filtered,'1',0,5,ts_mod(1:size(qrs_filtered{1},1)));

% hold on
% plot_ecgs_2x6_bios(full_ecg,'1',0,3,ts_mod,[r_idx(:,l) ;r_idx(:,l)],[s_idx(:,l); s_idx(:,l)]);

legend({strcat('Patient Nr',num2str(l))},'Location','eastoutside','Fontsize',30)
%saveas(gcf,strcat('C:\users\petnov\Dropbox\ecgs3\',num2str(l),'.png'));
close all
end

labs={'control and 50% global INa reduction','control', 'slow purkinje and 50% INa reduction','slow Purkinje'};
labs2 = {'control','RV slow', 'RV slower','RV slowest'};
labs3 ={'control','low density LGE','intermediate density LGE','high density LGE'};
legend({labs3{i-2}},'Interpreter', 'none','Fontsize',25);

tad_sim=[];
tad_sim={};
for j=3:6
     ecg=load([pop_path names(j).name]);
         max_ampl=max(max(ecg.vars));
    ecg=ecg.vars/max_ampl;
     ecg_cell={};
    qrs_sim= qrs_estimator(ecg(:,i));
             [~,s_ids_sim]=min(ecg);
        [~,r_ids_sim]=max(ecg);
        dummy=abs(size(ecg,1)-s_ids_sim(7:10))';
        tad_sim(:,j-2)= dummy;
%      for i=7:12
% %          qrs_sim(j)=size(ecg,1);
% %              qrs_sim(i-6,j-2)= qrs_estimator(ecg(:,i));
% %               [~,~,duration(i-6,j-2)]= s_wave_area(ecg(:,j));
% %          if r_ids_sim(i)>s_ids_sim(i)
% %              r_ids_sim(i)=1;
% %          end
%          [~,s_start_sim(i-6,j-2)]=min(abs(ecg(r_ids_sim(i):s_ids_sim(i),i)));
%         s_start_sim(i-6,j-2)=  s_start_sim(i-6,j-2)+r_ids_sim(i);
%      end
        %discrepancy metric
        discrepancy(j-2) = (sum(abs(QRS_ecgi(1:6,str2num(patient_nr))-qrs_sim)./QRS_ecgi(1:6,str2num(patient_nr)))+sum((abs(duration_ecgi(1:4,str2num(patient_nr))-tad_sim(:,j-2)))./duration_ecgi(1:4,str2num(patient_nr))))/10;
        
    for k=7:12
        mxs(k-6)=max(full_ecg_filtered_high{k});
    end
        for k=7:12
            
          full_ecg_filtered_high{k}=full_ecg_filtered_high{k}/max(mxs);
        end
        Ie=ecg(:,1)/1000; 
        IIe=ecg(:,2)/1000;
        IIIe=ecg(:,3)/1000;
        aVRe=ecg(:,4)/1000;
        aVLe=ecg(:,5)/1000;
        aVFe=ecg(:,6)/1000;
        V1e=ecg(:,7)/1000;
        V2e=ecg(:,8)/1000;
        V3e=ecg(:,9)/1000;
        V4e=ecg(:,10)/1000;
        V5e=ecg(:,11)/1000;
        V6e=ecg(:,12)/1000;
        comparison_points=[];
        
        [~,s_ids_sim]=min(ecg);
        duration(1:4,j)=size(ecg,1)-s_ids_sim(7:10);
        [~,r_ids_sim]=max(ecg);
        ref_points=[s_ids_sim(7:10),r_ids_sim(11:12)];
        ts_ecg=[ts_mod(s_nadir_idx(1:4,str2num(patient_nr))),ts_mod(r_peak_idx(5:6,str2num(patient_nr)))] ;%change the patient nr!!!  
        ids_ecg = [s_nadir_idx(1:4,str2num(patient_nr))',r_peak_idx(5:6,str2num(patient_nr))'];
        t_start_ecg= ts_ecg-ref_points;
        t_end_ecg = ts_ecg+(size(ecg,1)-ref_points);
    for i=7:12
        comparison_points=[];
        for z=1:size(ecg,1)
            
           n=t_start_ecg(i-6)+z;
    
          dummy= find(abs(ts_mod-n)<1);
          
           comparison_points=[comparison_points,dummy(1)];
        end
        [~,start_t_id]=min(abs(ts_mod-t_start_ecg(i-6)));
        dummy=ts_mod(start_t_id:s_idx(i-6,str2num(patient_nr)))-ts_mod(start_t_id);
        
                 %try finding the comparison points which start at the s-wave start
         %(V=0)
         [~,s_start]=min(abs(ecg(r_ids_sim(i):s_ids_sim(i),i)));
         s_start=s_start+r_ids_sim(i);
         if isempty(s_start)
             s_start=1;
         end
             
       % ecg_inter(:,i)=interp1(1:size(ecg,1),ecg(:,j),ts_mod(start_t_id)-ts_mod(start_t_id):dt:ts_mod(s_idx(i-6,str2num(patient_nr)))-ts_mod(start_t_id));
        figure()
      %  plot(full_ecg_filtered_high{i}(start_t_id:s_idx(i-6,str2num(patient_nr))))
        plot(full_ecg_filtered_high{i}(comparison_points(s_start:end)));
        hold on
         plot(ecg(s_start:end,i))
         close()
        
         %try finding the comparison points which start at the s-wave start
         %(V=0)
        abs_dist(i-6,j-2)= sum(abs(full_ecg_filtered_high{i}(comparison_points)-ecg(:,i)))/size(ecg,1);
        std_dist(i-6,j-2)=std(abs(full_ecg_filtered_high{i}(comparison_points)-ecg(:,i)));
    corr_coef(i-6,j-2)=sum(full_ecg_filtered_high{i}(comparison_points(s_start:end)).*ecg(s_start:end,i))/sum(sqrt(full_ecg_filtered_high{i}(comparison_points(s_start:end)).^2.*ecg(s_start:end,i).^2));
    end
end

corr_coef_mean = sum(corr_coef)/6;

end
if isempty(idx)==0
else
    legend(model_names{end},'Interpreter', 'none');
end
%legend({'Intermediate/Advanced ARVC phenotype','Mild ARVC Phenotype'},'Interpreter', 'none','Fontsize',21);




%ecg biomarkers from sims
model_names ={};
R_peak=[];
QRSs=[];
duration=[];
ratio=[];
S_start=[];
for i=3:6
    idx=[];
    ecg=load([pop_path names(i).name]);
     %model_names{i-2}=names(i).name;

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
      % [Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e]=ecgi_ECG_sim_plot_mat(ecg1,V3e_f,0.1);
    ecgs{i-2}=[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e];
    full_ecg={};
        full_ecg={Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e};

     %   plot_ecgs_2x6(full_ecg',names(3).name,0,z,1:size(V1e,1));
   amplitude_V1(i-2)=max(V1e);
    precordial=[V1e,V2e,V3e,V4e,V5e,V6e];
    for j=1:6
       % [QRSs(j,i-2)]=qrs_estimator(precordial(:,j));
        [~,ratio(j,i-2),~,~,~,~,~,R_peak(j,i-2)]=qrs_estimator(precordial(:,j));
        [ ~,~,duration(j,i-2),swave(j,i-2),QRSs(j,i-2),S_start(j,i-2)]=s_wave_area(precordial(:,j));
    end
   
end






        inf= h5info(strcat(pop_path,'results.h5'));

    upstroke_con= h5read(strcat(pop_path,'results.h5'),'/UpstrokeTimeMap_0',[1,1,1] ,[1,size(H_mesh.xyz,1),1]);
    permutation=load(strcat(pop_path,'permutation.txt'));
    upstroke_con=upstroke_con(permutation+1);
%ecg_diff=[ecgs{2}-ecgs{1}];
%figure()
%ecgi_ECG_sim_plot_mat(ecg_diff,V3e_f,i/4);
 offset=150;
    offset2=size(Ie_f,1)-300;
    full_ecg_ecgi={real(Ie_f(offset:2:offset2)),real(IIe_f(offset:2:offset2)),real(IIIe_f(offset:2:offset2)),real(aVRe_f(offset:2:offset2)),real(aVLe_f(offset:2:offset2)),real(aVFe_f(offset:2:offset2)),real(V1e_f(offset:2:offset2)),real(V2e_f(offset:2:offset2)),real(V3e_f(offset:2:offset2)),real(V4e_f(offset:2:offset2)),real(V5e_f(offset:2:offset2)),real(V6e_f(offset:2:offset2))};
 %plot clinical ecgs
    for m=1:12
          full_ecg_ecgi{m}=full_ecg_ecgi{m}-full_ecg_ecgi{m}(1)*ones(size(full_ecg_ecgi{m},1),1);
    end
    figure('Position',[0 0 500 1000]); %left bottom width height
    
        plot_ecgs_2x6(full_ecg_ecgi,names(3).name,0,5,ts_ecgi(offset:2:offset2)-ts_ecgi(offset));

patient_sim_folder='/run/media/petnov/Seagate Expansion Drive/peter data backup windows PC/patients/patient06/';
H_mesh = read_CHASTE( strcat(patient_sim_folder,    'chaste',patient_nr,'/HEART'));
H_mesh.xyzActs=upstroke{1};
H_mesh.xyzActsNew=upstroke{2};

write_VTK(H_mesh,strcat(sim_path,'control.vtk'),'binary');
sim_path='/run/media/petnov/Seagate Expansion Drive/peter data backup windows PC/patients/patient06/results/Eikonal/figure_dense_acts/';

          figure()
labs2={'control','control + LGE', 'slow purkinje','slow Purkinje + LGE'};

    %labels={'no LGE', 'LGE max conduction scaling factor 0.5','LGE max conduction scaling factor 0.4','LGE max conduction scaling factor 0.25'};
for i=1:4
subplot_tight(2,2,i)
    
    hold on
    patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',[upstroke{i}],'facecolor','interp')
   % patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',[upstroke_con]','facecolor','interp')

    axis off
    caxis([0 100])
    headlight
    title(['activation times [ms] ' ' ' labs2{i}],'FontSize',25)
    view ([ -6 58])
    colormap(flipud(jet))
end
%colorbar('southoutside')
%set(gca,'Fontsize',25)




ecg_ecgi=[full_ecg_ecgi{6},full_ecg_ecgi{7},full_ecg_ecgi{8},full_ecg_ecgi{9},full_ecg_ecgi{10},full_ecg_ecgi{11},full_ecg_ecgi{12}];
for j=1:6
      [QRSs_ecgi(j)]=qrs_estimator_bsp(ecg_ecgi(:,j));
        [~,ratio_ecgi(j),~,~,~,~,~,~]=qrs_estimator_bsp(ecg_ecgi(:,j));
        [ ~,~,duration_ecgi(j)]=s_wave_area_bsp_new(ecg_ecgi(:,j),1:size(ecg_ecgi,1)); 
end

% difference_QRS=[];
% difference_S_wave_dur=[];
% difference_RS_ratio=[];
% for i=1:size(QRSs,2)
%     difference_QRS(:,i)= QRSs(:,i)'-QRSs_ecgi;
%     difference_S_wave_dur(:,i)=duration(1:4,i)-duration_ecgi(1:4)';
%     difference_RS_ratio(:,i)=ratio(:,i)-ratio_ecgi';
% end

%cum_RS_ratio= sum(abs(difference_RS_ratio),1);

averaged_QRSs=mean(QRSs,1);
averaged_QRSs_ecgi=mean(QRSs_ecgi);
averaged_S_wave=mean(duration(1:4,:),1);



%diff_averaged_qrs=averaged_QRSs-averaged_QRSs_ecgi;

%averaged_Swave_diff=mean(abs(difference_S_wave_dur),1);
%averaged_QRS_diff=mean(abs(difference_QRS),1);
%good_models=find(averaged_Swave_diff<25);
%good_QRS_models=find(abs(diff_averaged_qrs)<15);
avg_s_wave=mean(swave(1:4,:),1);
%v_good_models=find(abs(diff_averaged_qrs)<25 & abs(averaged_Swave_diff)<15 & avg_s_wave>0.2);
deep_s_waves=find(avg_s_wave>0.5);
RS_ratio_mean=mean(ratio(1:4,:),1);
%score the models

qrs_score= abs(averaged_QRSs)>130;
s_wave_score= abs(averaged_S_wave)>80;
%s_wave_ampl_score=avg_s_wave>0.5;
ratio_score=RS_ratio_mean<0.7;
R_prog_tot=zeros(5,size(QRSs,2));
for i=1:5
        Rprog{i}=R_peak(i+1,:)>R_peak(i,:);
          R_prog_tot(i,:)=Rprog{i};

 
end
R_wave_prog=[];
for i=1:size(R_prog_tot,2)
    ids=strfind(R_prog_tot(:,i)',[1 1 1]);
    if isempty(ids)
        R_wave_prog(i)=0;
    
    else
        R_wave_prog(i)=1;
    end

end

%s wave start 
avg_s_wave_start=sum(S_start(1:4,:),1)/4;

s_start_score=avg_s_wave_start>80;

r_score=sum(R_peak<0.5,1)>3;

tot_score= int8(qrs_score)+int8(s_wave_score)+int8(ratio_score) +int8(s_start_score)+int8(R_wave_prog)+int8(r_score);

itemized_score=vertcat(int8(qrs_score),int8(s_wave_score),int8(ratio_score),int8(R_wave_prog),int8(s_start_score));

highest_scoring=find(tot_score>1);

long_APD= find(averaged_QRSs>200);
%plot only prolongued activations


long_params={};
set(0,'DefaultTextInterpreter','none')

leg_names={};
l=1;
figure()
for i =highest_scoring
    
    ecg=load([pop_path names(i+2).name]);
    
 %   [Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e]=ecgi_ECG_sim_plot_mat(ecg,V3e_f,l);
            full_ecg=ecgs{i};
        full_ecg={full_ecg(:,1),full_ecg(:,2),full_ecg(:,3),full_ecg(:,4),full_ecg(:,5),full_ecg(:,6),full_ecg(:,7),full_ecg(:,8),full_ecg(:,9),full_ecg(:,10),full_ecg(:,11),full_ecg(:,12)};
     plot_ecgs_2x6(full_ecg',names(3).name,0,l,1:size(full_ecg{1},1));
     hold on;
    long_params=vertcat(long_params,names(i+2).name);
    leg_names{l}=string(long_params{l});
    l=l+1;
end
leg=legend(leg_names,'Interpreter','none')
%set(leg,'Interpreter','none')





figure()
b=bar([horzcat(duration(1:4,:),duration_ecgi(1:4,9))]);
b(1).FaceColor='b';
b(2).FaceColor='c';
b(3).FaceColor='g';
b(4).FaceColor=[0.91 0.41 0.17];
b(5).FaceColor='r';

ylabel('TAD [ms]','FontSize',25)
set(gca,'xtickLabel',{'V1','V2','V3','V4'})
%leg=legend({string(long_params{1}(1:end-53)),string(long_params{2}(1:end-53)),string(long_params{3}(1:end-53)),string(long_params{4}(1:end-53)),string(long_params{5}(1:end-53)),'clinical recordings'})
set(gca,'FontSize',35)
%set(leg,'Interpreter','none')
legend(labs3,'clinical recordings')

QRS_ecgi(6,9)=89;

figure()
b=bar([horzcat(QRSs(1:6,:),QRS_ecgi(1:6,9))])
b(1).FaceColor='b';
b(2).FaceColor='c';
b(3).FaceColor='g';
b(4).FaceColor=[0.91 0.41 0.17];
b(5).FaceColor='r';
ylabel('QRS duration [ms]','FontSize',18)
set(gca,'xtickLabel',{'V1','V2','V3','V4','V5','V6'})
%leg=legend({string(long_params{1}(1:end-53)),string(long_params{2}(1:end-53)),string(long_params{3}(1:end-53)),string(long_params{4}(1:end-53)),string(long_params{5}(1:end-53)),'clinical recordings'})
set(gca,'FontSize',35)
%set(leg,'Interpreter','none')
labs3{5}='clinical QRS';
legend(labs3,'Location','northeastoutside')




figure()
bar([vertcat(ratio(1:4,1:4)',ratio_ecgi(1:4))]')
ylabel('R wave to S wave amplitude ratio','FontSize',30)
set(gca,'xtickLabel',{'V1','V2','V3','V4'})
%leg=legend({string(long_params{1}(1:end-53)),string(long_params{2}(1:end-53)),string(long_params{3}(1:end-53)),string(long_params{4}(1:end-53)),string(long_params{5}(1:end-53)),'clinical recordings'})
set(gca,'FontSize',30)
%set(leg,'Interpreter','none')
legend(labels,'clinical recordings')
ylim([0 3])








good_R_prog=find(R_wave_prog==1);
% figure()
% for i =long_APD
%     try
%     ecg=load([pop_path names(i+2).name]);
%     ecg=ecg.vars;
%     catch
%         continue
%     end
%     [Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e]=ecgi_ECG_sim_plot_mat(ecg,V3e_f,1);
%     long_params=vertcat(long_params,names(i+2).name);
% end

acts=dlmread(strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\results\',sim_name,'\','exp23\13.5x_3f_1.0fiber_roots_0_no_rv_pur_0\13.5x_3f_1.0fiber_roots_0_no_rv_pur_0_NEW_pred_ATMap','.csv'));
acts2=dlmread(strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\results\',sim_name,'\','exp17\11.0x_2f_septal_0.3fiber_roots_7_no_rv_pur_0_NEW_pred_ATMap','.csv'));
acts2=dlmread(strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\results\',sim_name,'\','exp20\11.0x_3f_0.3fiber_roots_7_no_rv_pur_0_NEW_pred_ATMap','.csv'));

sim_path='D:\ARVC meshing automatic\patients\patient06\results\Eikonal\figure_fib\';

folders_acts=dir(sim_path);
for i = [1:4]%3:nb

    names_acts=dir(strcat(sim_path,folders_acts(i+2).name));
    upstroke{i}=load(strcat(sim_path,folders_acts(i+2).name,'\',names_acts(3).name));
    
end

if plt_acts==1
 
H_mesh = read_CHASTE( strcat(patient_sim_folder,    'chaste',patient_nr,'\HEART'));
H_mesh2 = read_CHASTE( strcat(patient_sim_folder,    'chaste',patient_nr,'\HEART'));

 figure()
    headlight
    hold on
    patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',upstroke{2},'facecolor','interp')
    colorbar
    caxis([0 100])
    axis equal
end
pat_nr='06';
electrodes = load(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\mpp\ECG_ELECTRODES.mat'));

E = [ electrodes.ELECTRODES.LA ; electrodes.ELECTRODES.RA ; electrodes.ELECTRODES.LL ; electrodes.ELECTRODES.RL ; electrodes.ELECTRODES.V1 ; electrodes.ELECTRODES.V2 ; electrodes.ELECTRODES.V3 ; electrodes.ELECTRODES.V4 ; electrodes.ELECTRODES.V5 ; electrodes.ELECTRODES.V6 ];
E=E/10; %chaste is in cm
hplot3d( E(5:7,:) , '*m' ,'LineWidth',4,'MarkerSize',10 );




H_mesh.xyzacts=upstroke{1};
write_VTK(H_mesh,strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\results\',sim_name,'\','control_eikonal','.vtk'),'binary');
%plot all the candidate models from across pops
figure()
hold on
for l=1:size(keep_models,1)
    if l<3
        k=1;
    else
        k=1.5;
    end
    [Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e]=ecgi_ECG_sim_plot_mat(keep_models{l},V3e_f,k);
    legend( 'all roots','all roots','5 roots','5 roots')
end

names_cell=cellstr({names(3:end).name});
names_cell=names_cell;
model_nr=1:size(names_cell,2);
tbl=table(names_cell,averaged_QRSs,avg_s_wave,RS_ratio_mean,model_nr)


% \9
 figure()
 ecgi_plot('D:\Data_ARVC_raw_electrodes_vest\Data_ARVC\ARVC6\ARVC6_Baseline.mat','D:\ChrisReconstructions\ARVC6\ARVC6_Baseline.mat','D:\ARVC meshing automatic\patients\patient06\mpp\ECG_ELECTRODES.mat',6);
 hold on
 legend('recorded advanced ARVC phenotype')
 
 
 i=highest_scoring(1)+2;
    ecg=load([pop_path names(i).name]);
    ecg=ecg.vars;
  
    [Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e]=ecgi_ECG_sim_plot_mat(ecg,V3e_f,1);
    legend('clinical recording','best performing model')
    
figure()
plot_ecgi_rest('D:\ChrisReconstructions\ARVC6\ARVC6_Baseline.mat')
    caxis([40 46])



figure()
no_sodium=strfind(names_cell,'0f_1.0fiber');
idx=[];
ylabels={};
for i=1:size(no_sodium,2)
    if no_sodium{i}==7
        idx=[idx,i];
    end
end
h=heatmap(vertcat(averaged_QRSs(idx)/max(averaged_QRSs(idx)),avg_s_wave(idx)/max(avg_s_wave(idx)),RS_ratio_mean(idx)/max(RS_ratio_mean(idx))),names_cell(idx),[],[],'TickAngle',45,'ShowAllTicks',true)
set(gca,'ytick',[0,1,2,3])
set(gca,'yticklabel',{'','average QRS duration','average S - wave duration ','R-S wave ratio'})
colormap gray
colorbar


figure()
for i =idx(11)
    try
    ecg=load([pop_path names(i+2).name]);
    ecg=ecg.vars;
    catch
        continue
    end
    [Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e]=ecgi_ECG_sim_plot_mat(ecg,V3e_f,i/max(highest_scoring));

end

electrodes = load('D:\ARVC meshing automatic\patients\patient06\mpp\ECG_ELECTRODES.mat');


E = [ electrodes.ELECTRODES.LA ; electrodes.ELECTRODES.RA ; electrodes.ELECTRODES.LL ; electrodes.ELECTRODES.RL ; electrodes.ELECTRODES.V1 ; electrodes.ELECTRODES.V2 ; electrodes.ELECTRODES.V3 ; electrodes.ELECTRODES.V4 ; electrodes.ELECTRODES.V5 ; electrodes.ELECTRODES.V6 ];
E=E/10; %chaste is in cm





