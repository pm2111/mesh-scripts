addpath(genpath('C:/Users/petnov/Dropbox/'));
addpath(genpath('C:/Users/petnov/Dropbox/qrs/'));
addpath(genpath('C:/Users/petnov/Dropbox/shared - copy/mesh scripts/'));

addpath('C:\Users\peter\Dropbox\sim3d scripts\');
addpath('C:\Users\petnov\Dropbox\shared - copy\IO\');

close all
[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e,Ie_f,IIe_f,IIIe_f,aVRe_f,aVLe_f,aVFe_f,V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f,ts_ecgi]=ecgi_ECG_pat09('D:\Data_ARVC_raw_electrodes_vest\Data_ARVC\ARVC9\ARVC9_Baseline.mat','D:\ChrisReconstructions\ARVC9\ARVC9_Baseline.mat','D:\ARVC meshing automatic\patients\patient09\mpp\ECG_ELECTRODES.mat',9);
ecg_ecgi=[V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f];

sim_path='D:\ARVC meshing automatic\patients\patient09\results\Eikonal\exp7\';

plt_acts=0;
keep_models={};
patient_nr='09';
sim_name='Eikonal';
study='exp1';
patient_sim_folder =  strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\');
%H_mesh = read_CHASTE( strcat(patient_sim_folder,    'chaste',patient_nr,'\HEART'));
for z=[1]
    if z==1
        pop_path='D:\ARVC meshing automatic\patients\patient09\results\Eikonal\exp7_ECGs\';
    else
        pop_path='D:\ARVC meshing automatic\patients\patient06\results\Eikonal\exp9_ECGs\';
    end
names = dir(pop_path);
duration=[];
QRSs=[];
ratio=[];
R_peak=[];
nb=size(names,1);
ecgs={};
% ecg1=load(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',sim_name,'\',sim_type,'.mat'));
% ecg1=ecg1.ecg_full;
   figure('Position',[0 0 500 1000]); %left bottom width height
   amplitude_V1=[];
for i = 5:6
    
    ecg=load([pop_path names(i).name]);
    ecg=ecg.vars;
%     
%         ecg=load([pop_path names(i).name]);
%         ecg=ecg.ecg_full;
%     
    %[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e]=ecgi_ECG_sim_plot_mat(ecg,V3e_f,i-2);
        %[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e]=ecgi_ECG_sim_scale(ecg,V3e_f,i-2);
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
            hold off

        plot_ecgs_2x6(full_ecg',names(3).name,0,i-2,1:size(V1e,1));
        amplitude_V1(i-2)=max(V1e);
    precordial=[V1e,V2e,V3e,V4e,V5e,V6e];
    for j=1:6
       % [QRSs(j,i-2)]=qrs_estimator(precordial(:,j));
        [~,ratio(j,i-2),~,~,~,~,~,R_peak(j,i-2)]=qrs_estimator(precordial(:,j));
        [ ~,~,duration(j,i-2),swave(j,i-2),QRSs(j,i-2)]=s_wave_area(precordial(:,j));
    end
    
end
labs4={'no LGE','low density LGE ', 'intermediate density LGE','high density LGE'};

    legend(labs4(3:4),'FontSize',35)
    
    legend({'no LGE','low density LGE '},'FontSize',35)


long_params={};
set(0,'DefaultTextInterpreter','none')

tot_ecg=cell(1);
for i=1:12
    tot_ecg{i}=zeros(numel(full_ecg_ecgi{1}),1);
end
leg_names={};
l=1;
figure()
for m =3:nb
    
    ecg=load([pop_path names(m).name]);
    
 %   [Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e]=ecgi_ECG_sim_plot_mat(ecg,V3e_f,l);
        
            full_ecg=ecgs{m-2};
        full_ecg={full_ecg(:,1),full_ecg(:,2),full_ecg(:,3),full_ecg(:,4),full_ecg(:,5),full_ecg(:,6),full_ecg(:,7),full_ecg(:,8),full_ecg(:,9),full_ecg(:,10),full_ecg(:,11),full_ecg(:,12)};
    % plot_ecgs_2x6(full_ecg',names(3).name,0,l,1:size(full_ecg{1},1));
     hold on;
  %  long_params=vertcat(long_params,names(m).name);
    %leg_names{l}=string(long_params{l});
    l=l+1;

%leg=legend(leg_names,'Interpreter','none')
%set(leg,'Interpreter','none')

begin=15;
    finish=16;
    begin2=22;
    finish2=23;
filtered = cell(1);
i=1;
for x =[full_ecg{1},full_ecg{2},full_ecg{3},full_ecg{4},full_ecg{5},full_ecg{6},full_ecg{7},full_ecg{8},full_ecg{9},full_ecg{10},full_ecg{11},full_ecg{12}]
    F2=fft(x);
   if i==12
        F2(begin:finish)=0;
           F2(end-finish:end-begin)=0;
                   F2(begin2:finish2)=0;
    
    F2(end-finish2:end-begin2)=0;
   elseif i==5
       %figure()
       %plot(real(F2))
      F2(5:10)=0;
        
      F2(end-10:end-5)=0;
                F2(begin2:finish2)=0;
    
    F2(end-finish2:end-begin2)=0;
       
    else
        F2(begin:finish)=0;
        F2(end-finish:end-begin)=0;
                F2(begin2:finish2)=0;
    
    F2(end-finish2:end-begin2)=0;

    end

    filtered{i}=ifft(F2);
    filtered{i}(numel(full_ecg_ecgi{1}))=0;
  
        %tot_ecg{i}=tot_ecg{i}+filtered{i};
    
        i=i+1;

end
      plot_ecgs_2x6(filtered',names(3).name,0,5,1:size(full_ecg_ecgi{1},1));

end



    %labels={'no LGE','LGE max scaling factor 0.58 + septal CV slowing factor 0.5','LGE max scaling factor 0.5 + septal CV slowing factor 0.33','LGE max scaling factor 0.4 + septal CV slowing factor 0.25'};
    labels={'no LGE','LGE max scaling factor 0.58 ','LGE max scaling factor 0.5 ','LGE max scaling factor 0.4'};


    legend(labels{2},'FontSize',25)
    
    set(gca,'Fontsize',25)
    
    offset=230;
    offset2=size(Ie_f,1)-30;
    full_ecg_ecgi={real(Ie_f(offset:2:offset2)),real(IIe_f(offset:2:offset2)),real(IIIe_f(offset:2:offset2)),real(aVRe_f(offset:2:offset2)),real(aVLe_f(offset:2:offset2)),real(aVFe_f(offset:2:offset2)),real(V1e_f(offset:2:offset2)),real(V2e_f(offset:2:offset2)),real(V3e_f(offset:2:offset2)),real(V4e_f(offset:2:offset2)),real(V5e_f(offset:2:offset2)),real(V6e_f(offset:2:offset2))};
 %plot clinical ecgs
    figure('Position',[0 0 500 1000]); %left bottom width height
    
        plot_ecgs_2x6(full_ecg_ecgi,names(3).name,0,5,ts_ecgi(offset:2:offset2)-ts_ecgi(offset));

    
upstroke={};
folders_acts = dir(sim_path);

for i=1:4
    names_acts=dir(strcat(sim_path,folders_acts(i+2).name));
    upstroke{i}=load(strcat(sim_path,folders_acts(i+2).name,'\',names_acts(3).name));
end
    %labels={'no LGE','LGE max scaling factor 0.58 + septal CV slowing factor 0.5','LGE max scaling factor 0.5 + septal CV slowing factor 0.33','LGE max scaling factor 0.4 + septal CV slowing factor 0.25'};


E = [ electrodes.ELECTRODES.LA ; electrodes.ELECTRODES.RA ; electrodes.ELECTRODES.LL ; electrodes.ELECTRODES.RL ; electrodes.ELECTRODES.V1 ; electrodes.ELECTRODES.V2 ; electrodes.ELECTRODES.V3 ; electrodes.ELECTRODES.V4 ; electrodes.ELECTRODES.V5 ; electrodes.ELECTRODES.V6 ];
E=E/10; %chaste is in cm

    figure()
z=1;
for i=[1 2 3 4]
subplot_tight(2,2,z)
    title(strcat(labels{i}),'FontSize',25)

    hold on
    patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',[upstroke{i}],'facecolor','interp')
    [h,v] = IsoLine({ H_mesh.face, H_mesh.xyz},[upstroke{i}]',[10,20,30,40,50,60,70,80,90,100,110,120],'w');
    %catter3(E(5,1),E(5,2),E(5,3),105,'m','filled')
    axis off
    caxis([0 75])
    colormap(flipud(jet))
        headlight
  title(labs4{i},'FontSize',40,'Interpreter','none')
z=z+1;
end
colorbar('south')
set(gca,'Fontsize',25)
    
    
%ecg_diff=[ecgs{2}-ecgs{1}];
%figure()
%ecgi_ECG_sim_plot_mat(ecg_diff,V3e_f,i/4);


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
    difference_S_wave_dur(:,i)=duration(1:4,i)-duration_ecgi(1:4)';
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

itemized_score=vertcat(int8(qrs_score),int8(s_wave_score),int8(s_wave_ampl_score),int8(ratio_score),int8(R_wave_prog));

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
%leg=legend(string(long_params{1}(1:end-53)),string(long_params{2}(1:end-53)),string(long_params{3}(1:end-53)),string(long_params{4}(1:end-53)),'Interpreter','none')
%set(leg,'Interpreter','none')

end



    labels={'no LGE', 'LGE max CV scaling factor 0.58','LGE max CV scaling factor 0.5','LGE max CV scaling factor 0.4'};

figure()
bar([horzcat(duration(1:4,:),zeros(4,1))])
ylabel('terminal activation duration [ms]','FontSize',40)
set(gca,'xtickLabel',{'V1','V2','V3','V4'})
%leg=legend({string(long_params{1}(1:end-53)),string(long_params{2}(1:end-53)),string(long_params{3}(1:end-53)),string(long_params{4}(1:end-53)),string(long_params{5}(1:end-53)),'clinical recordings'})
set(gca,'FontSize',40)
%set(leg,'Interpreter','none')
legend(labels)


figure()
bar([vertcat(QRSs(1:6,1:4)',zeros(6,1)')]')
ylabel('QRS duration [ms]','FontSize',40)
set(gca,'xtickLabel',{'V1','V2','V3','V4','V5','V6'})
%leg=legend({string(long_params{1}(1:end-53)),string(long_params{2}(1:end-53)),string(long_params{3}(1:end-53)),string(long_params{4}(1:end-53)),string(long_params{5}(1:end-53)),'clinical recordings'})
set(gca,'FontSize',40)
%set(leg,'Interpreter','none')
legend(labels)




figure()
bar([vertcat(ratio(1:4,1:4)',ratio_ecgi(1:4))]')
ylabel('R wave to S wave amplitude ratio','FontSize',30)
set(gca,'xtickLabel',{'V1','V2','V3','V4'})
%leg=legend({string(long_params{1}(1:end-53)),string(long_params{2}(1:end-53)),string(long_params{3}(1:end-53)),string(long_params{4}(1:end-53)),string(long_params{5}(1:end-53)),'clinical recordings'})
set(gca,'FontSize',30)
%set(leg,'Interpreter','none')
legend(labels,'clinical recordings')
ylim([0 3])
electrodes = load('D:\ARVC meshing automatic\patients\patient09\mpp\ECG_ELECTRODES.mat');





good_R_prog=find(R_wave_prog==1);


acts=dlmread(strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\results\',sim_name,'\','exp16\11.0x_1f_1.0fiber_roots_7_no_rv_pur_0_NEW_pred_ATMap','.csv'));
acts2=dlmread(strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\results\',sim_name,'\','exp17\11.0x_2f_septal_0.3fiber_roots_7_no_rv_pur_0_NEW_pred_ATMap','.csv'));
acts2=dlmread(strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\results\',sim_name,'\','exp20\11.0x_3f_0.3fiber_roots_7_no_rv_pur_0_NEW_pred_ATMap','.csv'));

if plt_acts==1
    patient_sim_folder='D:\ARVC meshing automatic\patients\patient09\';
    patient_nr='09';
H_mesh = read_CHASTE( strcat(patient_sim_folder,    'chaste',patient_nr,'\HEART'));
 figure()
    headlight
    hold on
    patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',upstroke{2},'facecolor','interp')
    colorbar
    caxis([0 60])
    axis off
end
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
 ecgi_plot('D:\Data_ARVC_raw_electrodes_vest\Data_ARVC\ARVC9\ARVC9_Baseline.mat','D:\ChrisReconstructions\ARVC9\ARVC9_Baseline.mat','D:\ARVC meshing automatic\patients\patient09\mpp\ECG_ELECTRODES.mat',9);
 hold on
 
 
 
 i=highest_scoring(1)+2;
    ecg=load([pop_path names(i).name]);
    ecg=ecg.vars;
  
    [Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e]=ecgi_ECG_sim_plot_mat(ecg,V3e_f,1);
    legend('clinical recording','best performing model')
    
figure()
plot_ecgi_rest('D:\ChrisReconstructions\ARVC9\ARVC9_Baseline.mat')
caxis([22 70])



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



electrodes = load('D:\ARVC meshing automatic\patients\patient06\mpp\ECG_ELECTRODES.mat');


E = [ electrodes.ELECTRODES.LA ; electrodes.ELECTRODES.RA ; electrodes.ELECTRODES.LL ; electrodes.ELECTRODES.RL ; electrodes.ELECTRODES.V1 ; electrodes.ELECTRODES.V2 ; electrodes.ELECTRODES.V3 ; electrodes.ELECTRODES.V4 ; electrodes.ELECTRODES.V5 ; electrodes.ELECTRODES.V6 ];
E=E/10; %chaste is in cm


 figure()
subplot_tight(2,1,1)
    hold on

    patch('vertices',H_mesh.xyz*10,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',[upstroke{1}],'facecolor','interp')
    [h,v] = IsoLine({ H_mesh.face, H_mesh.xyz*10},[upstroke{1}]',5,'k');
    %catter3(E(5,1),E(5,2),E(5,3),105,'m','filled')
    axis equal
    axis off
    caxis([31 32])
    colormap(flipud(jet))
        headlight
        view([29, 32])
                zoom(1.5)

subplot_tight(2,1,2)

    hold on
      plot_ecgi_rest('D:\ChrisReconstructions\ARVC9\ARVC9_Baseline.mat')
    axis equal
    caxis([37 37.5])
        view([29, 32])
    colormap(flipud(jet))
        headlight


