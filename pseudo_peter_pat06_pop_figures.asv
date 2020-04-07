addpath(genpath('C:/Users/petnov/Dropbox/'));
addpath(genpath('C:/Users/petnov/Dropbox/qrs/'));
addpath(genpath('C:/Users/petnov/Dropbox/sim 3d scripts/'));
addpath(genpath('C:/Users/petnov/Dropbox/shared - copy/mesh scripts/'));

addpath('C:\Users\peter\Dropbox\sim3d scripts\');
addpath('C:\Users\petnov\Dropbox\shared - copy\IO\');

close all
[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e,Ie_f,IIe_f,IIIe_f,aVRe_f,aVLe_f,aVFe_f,V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f,ts_ecgi]=ecgi_ECG06('D:\Data_ARVC_raw_electrodes_vest\Data_ARVC\ARVC6\ARVC6_Baseline.mat','D:\ChrisReconstructions\ARVC6\ARVC6_Baseline.mat','D:\ARVC meshing automatic\patients\patient06\mpp\ECG_ELECTRODES.mat',6);
ecg_ecgi=[V1e_f,V2e_f,V3e_f,V4e_f,V5e_f,V6e_f];
sim_path='D:\ARVC meshing automatic\patients\patient06\results\Eikonal\exp30\';
plt_acts=0;
keep_models={};
patient_nr='06';
sim_name='Eikonal';
study='exp32';
patient_sim_folder =  strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\');
%H_mesh = read_CHASTE( strcat(patient_sim_folder,    'chaste',patient_nr,'\HEART'));
for z=[1]
    if z==1
        pop_path='D:\ARVC meshing automatic\patients\patient06\results\Eikonal\test_coarse26_ECGs\';
    else
        pop_path='D:\ARVC meshing automatic\patients\patient06\results\Eikonal\figure_fib_ecgs\';
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
sim_path='D:\ARVC meshing automatic\patients\patient06\results\Eikonal\exp_roots_light_pheno\';

folders_acts=dir(sim_path);
for i = [1]%3:nb

    names_acts=dir(strcat(sim_path,folders_acts(i+2).name));
    upstroke{i}=load(strcat(sim_path,folders_acts(i+2).name,'\',names_acts(3).name));
    
end



    offset=80;
    offset2=size(Ie_f,1)-270;
    full_ecg_ecgi={real(Ie_f(offset:2:offset2)),real(IIe_f(offset:2:offset2)),real(IIIe_f(offset:2:offset2)),real(aVRe_f(offset:2:offset2)),real(aVLe_f(offset:2:offset2)),real(aVFe_f(offset:2:offset2)),real(V1e_f(offset:2:offset2))-real(V1e_f(offset)),real(V2e_f(offset:2:offset2))-real(V2e_f(offset)),real(V3e_f(offset:2:offset2))-real(V3e_f(offset)),real(V4e_f(offset:2:offset2))-real(V4e_f(offset)),real(V5e_f(offset:2:offset2))-real(V5e_f(offset)),real(V6e_f(offset:2:offset2))-real(V6e_f(offset))};
 %plot clinical ecgs
    figure('Position',[0 0 500 1000]); %left bottom width height
    
        plot_ecgs_2x6(full_ecg_ecgi,names(3).name,0,5,ts_ecgi(offset:2:offset2)-ts_ecgi(offset));

    for j=1:6
      [QRSs_ecgi(j)]=abs(qrs_estimator_bsp(ecg_ecgi(:,j))*(ts_ecgi(2)-ts_ecgi(1)));
        [~,ratio_ecgi(j),~,~,~,~,~,~]=qrs_estimator_bsp(ecg_ecgi(:,j));
        [ ~,~,duration_ecgi(j),~,RS_grad(j)]=s_wave_area_bsp(ecg_ecgi(:,j),1:size(ecg_ecgi,1)); 
    end



labels= {  'no LGE and all roots present'   'LBBB control' 'no LGE global conduction scaling 0.7 + LBBB'    'LGE regional conduction scaling + LBBB'    };



k=1;
model_names ={};
%for o=2:nb/2
figure('Position',[0 0 500 1000]); %left bottom width height
z=1;
for i =[3 :6]%o+(1:2)%o+122*(0:1)%o+(1:4)%[-5+o+ 5*(1:5)]%[2+o+36*(0:3)]% [-3+o+ 6*(1:6)]%[-5+o+ 5*(1:5)] %[o+(0:3)*25]%[-2+ 5*(1:5)] %3:nb
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

        plot_ecgs_2x6(full_ecg',names(3).name,0,z,1:size(V1e,1));
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
labs={'control and 50% global INa reduction','control', 'slow purkinje and 50% INa reduction','slow Purkinje'};

legend(labs2,'Interpreter', 'none','Fontsize',25);

%     if isempty(idx)==0
%         i=idx;
%          ecg=load([pop_path names(i).name]);
%              ecg=ecg.vars;
%              model_names{k}=names(i).name;
%     %[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e]=ecgi_ECG_sim_plot_mat(ecg,V3e_f,i-2);
%         %[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e]=ecgi_ECG_sim_scale(ecg,V3e_f,i-2);
%           Ie=ecg(:,1)/10; 
%         IIe=ecg(:,2)/10;
%         IIIe=ecg(:,3)/10;
%         aVRe=ecg(:,4)/10;
%         aVLe=ecg(:,5)/10;
%         aVFe=ecg(:,6)/10;
%         V1e=ecg(:,7)/10;
%         V2e=ecg(:,8)/10;
%         V3e=ecg(:,9)/10;
%         V4e=ecg(:,10)/10;
%         V5e=ecg(:,11)/10;
%         V6e=ecg(:,12)/10;
%       % [Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e]=ecgi_ECG_sim_plot_mat(ecg1,V3e_f,0.1);
%     ecgs{i-2}=[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e];
%     full_ecg={};
%         full_ecg={Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e};
% 
%         plot_ecgs_2x6(full_ecg',names(3).name,i,z,1:size(V1e,1));
%             hold on
%     end
%     for j=1:6
%         [QRSs(j,i-2)]=qrs_estimator(precordial(:,j));
%         [~,ratio(j,i-2),~,~,~,~,~,R_peak(j,i-2)]=qrs_estimator(precordial(:,j));
%         [ ~,~,duration(j,i-2),swave(j,i-2)]=s_wave_area(precordial(:,j));
%     end
%     k=k+1;
%     z=2;
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
for i=3:nb
    idx=[];
    ecg=load([pop_path names(i).name]);
     model_names{i-2}=names(i).name;

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
            hold on
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

patient_sim_folder='D:\ARVC meshing automatic\patients\patient06\';
H_mesh = read_CHASTE( strcat(patient_sim_folder,    'chaste',patient_nr,'\HEART'));
H_mesh.xyzActs=upstroke{1};
H_mesh.xyzActsNew=upstroke{2};

write_VTK(H_mesh,strcat(sim_path,'control.vtk'),'binary');
sim_path='D:\ARVC meshing automatic\patients\patient06\results\Eikonal\figure_dense_acts\';

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
bar([vertcat(duration(1:4,:)',duration_ecgi(1:4))]')
ylabel('S wave duration [ms]','FontSize',25)
set(gca,'xtickLabel',{'V1','V2','V3','V4'})
%leg=legend({string(long_params{1}(1:end-53)),string(long_params{2}(1:end-53)),string(long_params{3}(1:end-53)),string(long_params{4}(1:end-53)),string(long_params{5}(1:end-53)),'clinical recordings'})
set(gca,'FontSize',25)
%set(leg,'Interpreter','none')
legend(labels,'clinical recordings')


figure()
bar([vertcat(QRSs(1:6,:)',QRSs_ecgi(1:6))]')
ylabel('QRS duration [ms]','FontSize',18)
set(gca,'xtickLabel',{'V1','V2','V3','V4','V5','V6'})
%leg=legend({string(long_params{1}(1:end-53)),string(long_params{2}(1:end-53)),string(long_params{3}(1:end-53)),string(long_params{4}(1:end-53)),string(long_params{5}(1:end-53)),'clinical recordings'})
set(gca,'FontSize',18)
%set(leg,'Interpreter','none')
legend(labels,'clinical recordings')




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





