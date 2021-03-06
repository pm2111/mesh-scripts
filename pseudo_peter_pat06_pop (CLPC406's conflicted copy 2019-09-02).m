addpath(genpath('C:/Users/petnov/Dropbox/'));
addpath(genpath('C:/Users/petnov/Dropbox/qrs/'));
addpath(genpath('C:/Users/petnov/Dropbox/sim 3d scripts/'));
addpath(genpath('C:/Users/petnov/Dropbox/shared - copy/mesh scripts/'));

addpath('C:\Users\peter\Dropbox\sim3d scripts\');
addpath('C:\Users\petnov\Dropbox\shared - copy\IO\');

clear all
close all
[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e,Ie_f,IIe_f,IIIe_f,aVRe_f,aVLe_f,aVFe_f,V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f,ts_ecgi]=ecgi_ECG06('D:\Data_ARVC_raw_electrodes_vest\Data_ARVC\ARVC6\ARVC6_Baseline.mat','D:\ChrisReconstructions\ARVC6\ARVC6_Baseline.mat','D:\ARVC meshing automatic\patients\patient06\mpp\ECG_ELECTRODES.mat',6);
ecg_ecgi=[V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f];
pop_path='D:\ARVC meshing automatic\patients\patient06\results\Eikonal\exp15_ECGs\';
plt_acts=0;
keep_models={};
patient_nr='06';
sim_name='Eikonal';
study='exp15';
patient_sim_folder =  strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\');
%H_mesh = read_CHASTE( strcat(patient_sim_folder,    'chaste',patient_nr,'\HEART'));
for z=[1]
    if z==1
        pop_path='D:\ARVC meshing automatic\patients\patient06\results\Eikonal\exp154_ECGs\';
    else
        pop_path='D:\ARVC meshing automatic\patients\patient06\results\Eikonal\exp9_ECGs\';
    end
names = dir(pop_path);
duration=[];
QRSs=[];
ratio=[];
R_peak=[];
nb=size(names,1);
figure()
for i = 3:nb
    try
    ecg=load([pop_path names(i).name]);
    ecg=ecg.vars;
    catch
        continue
    end
    [Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e]=ecgi_ECG_sim_plot_mat(ecg,V3e_f,1);
    precordial=[V1e,V2e,V3e,V4e,V5e,V6e];
    for j=1:6
        [QRSs(j,i-2)]=qrs_estimator(precordial(:,j));
        [~,ratio(j,i-2),~,~,~,~,~,R_peak(j,i-2)]=qrs_estimator(precordial(:,j));
        [ ~,~,duration(j,i-2),swave(j,i-2)]=s_wave_area(precordial(:,j));
    end
end


for j=1:6
      [QRSs_ecgi(j)]=qrs_estimator_bsp(ecg_ecgi(:,j))*(ts_ecgi(2)-ts_ecgi(1));
        [~,ratio_ecgi(j),~,~,~,~,~,~]=qrs_estimator_bsp(ecg_ecgi(:,j));
        [ ~,~,duration_ecgi(j)]=s_wave_area_bsp(ecg_ecgi(:,j),ts_ecgi); 
end

difference_QRS=[];
difference_S_wave_dur=[];
difference_RS_ratio=[];
for i=1:size(QRSs,2)
    difference_QRS(:,i)= QRSs(:,i)'-QRSs_ecgi;
    difference_S_wave_dur(:,i)=duration(:,i)-duration_ecgi';
    difference_RS_ratio(:,i)=ratio(:,i)-ratio_ecgi';
end

cum_RS_ratio= sum(abs(difference_RS_ratio),1);

averaged_QRSs=mean(QRSs,1);
averaged_QRSs_ecgi=mean(QRSs_ecgi);
averaged_S_wave=mean(duration(1:4,:),1);



diff_averaged_qrs=averaged_QRSs-averaged_QRSs_ecgi;

averaged_Swave_diff=mean(abs(difference_S_wave_dur),1);
averaged_QRS_diff=mean(abs(difference_QRS),1);
good_models=find(averaged_Swave_diff<25);
good_QRS_models=find(abs(diff_averaged_qrs)<15);
avg_s_wave=mean(swave(1:4,:),1);
%v_good_models=find(abs(diff_averaged_qrs)<25 & abs(averaged_Swave_diff)<15 & avg_s_wave>0.2);
deep_s_waves=find(avg_s_wave>0.5);
RS_ratio_mean=mean(ratio(1:4,:),1);
%score the models

qrs_score= abs(diff_averaged_qrs)<25;
s_wave_score= abs(averaged_Swave_diff)<25;
s_wave_ampl_score=avg_s_wave>0.5;
ratio_score=RS_ratio_mean<0.7;
R_prog_tot=zeros(5,size(averaged_QRS_diff,2));
for i=1:5
        Rprog{i}=R_peak(i+1,:)>R_peak(i,:);
          R_prog_tot(i,:)=Rprog{i};

 
end
for i=1:size(R_prog_tot,2)
    ids=strfind(R_prog_tot(:,i)',[1 1 1]);
    if isempty(ids)
        R_wave_prog(i)=0;
    
    else
        R_wave_prog(i)=1;
    end

end
tot_score= int8(qrs_score)+int8(s_wave_score)+int8(s_wave_ampl_score)+int8(ratio_score) +int8(R_wave_prog);

highest_scoring=find(tot_score==max(tot_score));

long_APD= find(averaged_QRSs>200);
%plot only prolongued activations
long_params={};
set(0,'DefaultTextInterpreter','none')

figure()
for i =highest_scoring
    try
    ecg=load([pop_path names(i+2).name]);
    ecg=ecg.vars;
    catch
        continue
    end
    [Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e]=ecgi_ECG_sim_plot_mat(ecg,V3e_f,i*1.5/max(highest_scoring));
    long_params=vertcat(long_params,names(i+2).name);
    keep_models=vertcat(keep_models,ecg);
end
leg=legend(string(long_params{1}(1:end-53)),string(long_params{2}(1:end-53)),string(long_params{3}(1:end-53)),string(long_params{4}(1:end-53)),string(long_params{5}(1:end-53)),'Interpreter','none')
set(leg,'Interpreter','none')

end


figure()
bar([vertcat(duration(1:4,highest_scoring)',duration_ecgi(1:4))]')
ylabel('ms','FontSize',18)
set(gca,'xtickLabel',{'V1','V2','V3','V4'})
leg=legend({string(long_params{1}(1:end-53)),string(long_params{2}(1:end-53)),string(long_params{3}(1:end-53)),string(long_params{4}(1:end-53)),string(long_params{5}(1:end-53)),'clinical recordings'})
set(gca,'FontSize',18)
set(leg,'Interpreter','none')

figure()
bar([vertcat(QRSs(1:6,highest_scoring)',QRSs_ecgi(1:6))]')
ylabel('ms','FontSize',18)
set(gca,'xtickLabel',{'V1','V2','V3','V4','V5','V6'})
leg=legend({string(long_params{1}(1:end-53)),string(long_params{2}(1:end-53)),string(long_params{3}(1:end-53)),string(long_params{4}(1:end-53)),string(long_params{5}(1:end-53)),'clinical recordings'})
set(gca,'FontSize',18)
set(leg,'Interpreter','none')

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

acts=dlmread(strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\results\',sim_name,'\','17.0x_0f_1.0fiber\17.0x_0f_1.0fiber_NEW_pred_ATMap','.csv'));

if plt_acts==1
H_mesh = read_CHASTE( strcat(patient_sim_folder,    'chaste',patient_nr,'\HEART'));
 figure()
    headlight
    hold on
    patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',acts,'facecolor','interp')
    colorbar
end


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
% figure()
% ecgi_plot('D:\Data_ARVC_raw_electrodes_vest\Data_ARVC\ARVC6\ARVC6_Baseline.mat','D:\ChrisReconstructions\ARVC6\ARVC6_Baseline.mat','D:\ARVC meshing automatic\patients\patient06\mpp\ECG_ELECTRODES.mat',6);
% 

plot_ecgi_rest('D:\ChrisReconstructions\ARVC6\ARVC6_Baseline.mat')

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





