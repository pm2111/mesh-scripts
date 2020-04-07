addpath(genpath('C:/Users/petnov/Dropbox/'));
addpath(genpath('C:/Users/petnov/Dropbox/qrs/'));
addpath(genpath('C:/Users/petnov/Dropbox/sim 3d scripts/'));
addpath(genpath('C:/Users/petnov/Dropbox/shared - copy/mesh scripts/'));

addpath('C:\Users\petnov\Dropbox\sim3d scripts\');
addpath(genpath('C:\Users\petnov\Dropbox\mesh scripts\'));
addpath(genpath('C:\Users\petnov\Dropbox\qrs\'));
addpath(genpath('C:\Users\petnov\Dropbox\tools\'));

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
        pop_path='D:\ARVC meshing automatic\patients\patient06\results\Eikonal_coarse\control_calibration_ECGs\';
    else
        pop_path='D:\ARVC meshing automatic\patients\patient06\results\Eikonal_coarse\rv_posterior_scar_long_baseline_qrs_AT\';
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
sim_path='D:\ARVC meshing automatic\patients\patient06\results\Eikonal_coarse\rv_mid_scar_short_qrs_dense_scar_AT\';

folders_acts=dir(sim_path);
o=1;
for i =10%3:7
    names_acts=dir(strcat(sim_path,folders_acts(i).name));
    upstroke{o}=load(strcat(sim_path,folders_acts(i).name,'\',names_acts(3).name));
    o=o+1;
    
end


    offset=80;
    offset2=size(Ie_f,1)-270;
    full_ecg_ecgi={real(Ie_f(offset:2:offset2)),real(IIe_f(offset:2:offset2)),real(IIIe_f(offset:2:offset2)),real(aVRe_f(offset:2:offset2)),real(aVLe_f(offset:2:offset2)),real(aVFe_f(offset:2:offset2)),real(V1e_f(offset:2:offset2))-real(V1e_f(offset)),real(V2e_f(offset:2:offset2))-real(V2e_f(offset)),real(V3e_f(offset:2:offset2))-real(V3e_f(offset)),real(V4e_f(offset:2:offset2))-real(V4e_f(offset)),real(V5e_f(offset:2:offset2))-real(V5e_f(offset)),real(V6e_f(offset:2:offset2))-real(V6e_f(offset))};
 %plot clinical ecgs
    figure('Position',[0 0 500 1000]); %left bottom width height
    
        plot_ecgs_2x6(full_ecg_ecgi',names(3).name,0,5,ts_ecgi(offset:2:offset2)-ts_ecgi(offset));
        %plot_ecgs_2x6(full_ecg_ecgi,names(3).name,0,2,1:size(ts_ecgi(offset:2:offset2)-ts_ecgi(offset),2));

    for j=1:6
      [QRSs_ecgi(j)]=abs(qrs_estimator_bsp(ecg_ecgi(:,j))*(ts_ecgi(2)-ts_ecgi(1)));
        [~,ratio_ecgi(j),s_wave_ecg(j),~,~,~,~,~]=qrs_estimator_bsp(ecg_ecgi(:,j));
        [ ~,~,duration_ecgi(j),~,RS_grad(j)]=s_wave_area_bsp(ecg_ecgi(:,j),ts_ecgi); 
    end



    labels={'slow Purkinje CV','intermediate Purkinje CV','control Purkinje CV'};



k=1;
model_names ={};
%for o=2:nb/2
figure('Position',[0 0 500 1000]); %left bottom width height
z=1;
start=3;
fig_end=5;
for i =start:fig_end%[15:16:256]%[3:nb]%[27,25,47,45]%o+(1:2)%o+122*(0:1)%o+(1:4)%[-5+o+ 5*(1:5)]%[2+o+36*(0:3)]% [-3+o+ 6*(1:6)]%[-5+o+ 5*(1:5)] %[o+(0:3)*25]%[-2+ 5*(1:5)] %3:nb
    idx=[];
    ecg=load([pop_path names(i).name]);
    ecg_control=load('D:\ARVC meshing automatic\patients\patient06\results\Eikonal_coarse\figure_ecgs_fib_lat\17.0x_0f_1.0fiber_roots_RV_ant5.0_RV_post_3.0_RV_septal0.0_NEW_pred_ATMap.csvpseudo_ECG.mat');
     model_names{z}=names(i).name(1:end-4);
scal=30;
    ecg=ecg.vars;
       Ie=ecg(:,1)/scal; 
        IIe=ecg(:,2)/scal;
        IIIe=ecg(:,3)/scal;
        aVRe=ecg(:,4)/scal;
        aVLe=ecg(:,5)/scal;
        aVFe=ecg(:,6)/scal;
        V1e=ecg(:,7)/scal;
        V2e=ecg(:,8)/scal;
        V3e=ecg(:,9)/scal;
        V4e=ecg(:,10)/scal;
        V5e=ecg(:,11)/scal;
        V6e=ecg(:,12)/scal;
      % [Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e]=ecgi_ECG_sim_plot_mat(ecg1,V3e_f,0.1);
    ecgs{i-2}=[Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e];
    full_ecg={};

        full_ecg={Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e};
        for n=1:12
            full_ecg_control{n} = ecg_control.vars(:,n)/scal;
        end
              if i==start
                       plot_ecgs_2x6(full_ecg_control',names(3).name,0,8,1:size(full_ecg_control{1},1));
              end
                        hold on

           plot_ecgs_2x6(full_ecg',names(3).name,0,z,1:size(full_ecg{1},1));


 precordial=[V1e,V2e,V3e,V4e,V5e,V6e];
    for j=1:6
       % [QRSs(j,i-2)]=qrs_estimator(precordial(:,j));
        %[~,ratio(j,i-2),~,~,~,~,~,R_peak(j,i-2)]=qrs_estimator(precordial(:,j));
      %  [ ~,~,duration(j,i-2),swave(j,i-2),QRSs(j,i-2),S_start(j,i-2)]=s_wave_area(precordial(:,j));
    end
   
%    
%     
%     else
%         dummy_name=model_names{k-1};
%         dummy_name(7)='2';
%         A=strfind({names.name},dummy_name);
%         idx=find(not(cellfun('isempty',A)));
z=z+1;
end
%labs={'control and 50% global INa reduction','control', 'slow purkinje and 50% INa reduction','slow Purkinje'};
%labs2={'control','control + LGE', 'slow Purkinje','slow Purkinje + LGE'};
%labs3={'Slow Purkinje Speed','Intermediate Purkinje Speed', 'Control Purkinje Speed'};
%labs4={'no LGE scaling','max LGE scaling ', 'Control Purkinje Speed'};
legend(model_names,'Interpreter', 'none','Fontsize',25);

for n=7:12
    [ ~,~,duration_con(n-6),swave_con(n-6),QRSs_con(n-6),S_start_con(n-6)]=s_wave_area(full_ecg_control{n});
end

labs5={};
labs5{1}='control';

for i=1:fig_end-start+1
    labs5{i+1} = [num2str(6+i), 'mm thick scarred fiber with',' 80% CV slowing'];
end
legend(labs5,'Interpreter', 'none','Fontsize',25);


%plot the QRS feature for a posterior scarred fiber 
xlabels={'V1','V2','V3','V4','V5','V6'};
cols_colorbar={[0 0 0], [0 0 1], [0.54 0.82 0.99], [0.91 0.41 0.17]};
figure()
b=bar(horzcat(QRSs_con', QRSs(1:6,start-2:fig_end-2)),1,'FaceColor','flat')
%xticks([1.5 ,2.5 ,3.5])
xticklabels(xlabels);
ylabel(' QRS duration [ms]')
set(gca,'FontSize',35)
legend({'control','7mm thick scar','8mm thick scar','9 mm thick scar'},'Fontsize',35)
for i=1:4
b(i).CData=cols_colorbar{i};
end


%plot the TAD feature for a posterior scarred fiber 
xlabels={'V1','V2','V3','V4'};
cols_colorbar={[0 0 0], [0 0 1], [0.54 0.82 0.99], [0.91 0.41 0.17]};
figure()
b=bar(horzcat(duration_con(1:4)', duration(1:4,start-2:fig_end-2)),1,'FaceColor','flat')
%xticks([1.5 ,2.5 ,3.5])
xticklabels(xlabels);
ylabel(' TAD [ms]')
set(gca,'FontSize',35)
legend({'control','7mm thick scar','8mm thick scar','9 mm thick scar'},'Fontsize',35,'Location','northeastoutside')
for i=1:4
b(i).CData=cols_colorbar{i};
end





o=1;
for i=[1 3 5 8]
    legs{o} = ['epi/endo scar density ratio ',num2str(round(1/scar_den(1,i,3),1))];
    o=o+1;
end
legend(legs,'Interpreter', 'none','Fontsize',25);


o=1;
for i=1:5:25
    legs{o} = ['endo/epi scar density ratio ',num2str(round(scar_den(i,1,3),1))];
    o=o+1;
end
legend(legs,'Interpreter', 'none','Fontsize',25);




legend({'intact endocardium','endocardial scar'},'Interpreter', 'none','Fontsize',25);



legend({'15ms delay in scar','40ms delay in scar','60ms delay in scar','120ms delay in scar'},'Interpreter', 'none','Fontsize',25);


legend({'60% reduction in RV Purkinje Tree CV','40% reduction in RV Purkinje Tree CV', '20% reduction in RV Purkinje Tree CV','control'},'Interpreter', 'none','Fontsize',25);

legend({'Slow Purkinje + 2 mm RV mean wall thickness','Slow Purkinje + 3 mm RV mean wall thickness'},'Interpreter', 'none','Fontsize',25);

legend({'no LGE','low density LGE','intermediate density LGE','high density LGE'},'Interpreter', 'none','Fontsize',25);

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


averaged_QRSs=mean(QRSs,1);

averaged_S_wave=mean(duration(1:4,:),1);





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



folders_acts=dir(sim_path);
upstroke={};
for i = [1:3]%3:nb

    names_acts=dir(strcat(sim_path,folders_acts(i+2).name));
    upstroke{i}=load(strcat(sim_path,folders_acts(i+2).name,'\',names_acts(3).name));
end

labels={'CV in LGE MRI core 50%','CV in LGE MRI core 25% ','CV in LGE MRI core 12.5% ','CV in LGE MRI core 6.25% ','CV in LGE MRI core 3.125% of outer fibrotic region'} 
labels={'no LGE','low density LGE','intermediate density LGE','high density LGE'};
labels={'slow Purkinje CV','intermediate Purkinje CV', 'control Purkinje CV'};

labels={'15ms delay in scar','40ms delay in scar','60ms delay in scar','120ms delay in scar'};


labels={'60% reduction in RV Purkinje Tree CV','40% reduction in RV Purkinje Tree CV', '20% reduction in RV Purkinje Tree CV','control'};
figure()
    %labels={'no LGE', 'LGE max conduction scaling factor 0.5','LGE max conduction scaling factor 0.4','LGE max conduction scaling factor 0.25'};
for i=[1:4]
subplot(2,2,i)

    patch('vertices',H_mesh2.xyz,'faces',H_mesh2.face,'edgecolor','none','FaceVertexCData',[upstroke{i}],'facecolor','interp')
    %patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',[upstroke_con]','facecolor','interp')
        hold on
            [h,v] = IsoLine({ H_mesh2.face, H_mesh2.xyz},[upstroke{i}]',[10,20,30,40,50,60,70,80,90,100,110,120],'w');
 
    axis equal
    axis off
   
    headlight
  title(model_names{i},'FontSize',10,'Interpreter','none')
    view ([ -6 58])
   colormap(flipud(jet))
   caxis([0 75])
end




colorbar('south')
set(gca,'Fontsize',25)




ecg_ecgi=[full_ecg_ecgi{6},full_ecg_ecgi{7},full_ecg_ecgi{8},full_ecg_ecgi{9},full_ecg_ecgi{10},full_ecg_ecgi{11},full_ecg_ecgi{12}];
for j=1:6
      [QRSs_ecgi(j)]=qrs_estimator_bsp(ecg_ecgi(:,j));
        [~,ratio_ecgi(j),~,~,~,~,~,~]=qrs_estimator_bsp(ecg_ecgi(:,j));
        [ ~,~,duration_ecgi(j)]=s_wave_area_bsp_new(ecg_ecgi(:,j),ts_ecgi); 
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

averaged_QRSs_ecgi=mean(QRSs_ecgi);






%diff_averaged_qrs=averaged_QRSs-averaged_QRSs_ecgi;

%averaged_Swave_diff=mean(abs(difference_S_wave_dur),1);
%averaged_QRS_diff=mean(abs(difference_QRS),1);
%good_models=find(averaged_Swave_diff<25);
%good_QRS_models=find(abs(diff_averaged_qrs)<15);
avg_s_wave=mean(swave(1:4,:),1);
%v_good_models=find(abs(diff_averaged_qrs)<25 & abs(averaged_Swave_diff)<15 & avg_s_wave>0.2);
deep_s_waves=find(avg_s_wave>5);
RS_ratio_mean=mean(ratio(1:4,:),1);
%score the models

qrs_score= abs(averaged_QRSs)>130;
s_wave_score= abs(averaged_S_wave)>55;
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

s_start_score=avg_s_wave_start>40;

r_score=sum(R_peak<0.5,1)>2;

avg_R_peak=mean(R_peak(1:4,:));

r_s_ratio_value=avg_R_peak./avg_s_wave;
r_s_ratio=avg_R_peak./avg_s_wave<0.9;

tot_score= int8(qrs_score)+int8(s_wave_score)+int8(ratio_score) +int8(s_start_score)+int8(R_wave_prog)+int8(r_score)+int8(r_s_ratio);

itemized_score=vertcat(int8(qrs_score),int8(s_wave_score),int8(ratio_score),int8(R_wave_prog),int8(s_start_score));

highest_scoring=find(tot_score>=max(tot_score)-1);

long_APD= find(averaged_QRSs>200);
%plot only prolongued activations


long_params={};
set(0,'DefaultTextInterpreter','none')

tot_ecg=cell(1);
for i=1:12
    tot_ecg{i}=zeros(numel(full_ecg_ecgi{1}),1);
end
leg_names={};
l=1;
figure()
for k =[15 79 143 207]%highest_scoring(5:10)
    
    ecg=load([pop_path names(k+2).name]);
    
 %   [Ie,IIe,IIIe,aVRe,aVLe,aVFe,V1e,V2e,V3e,V4e,V5e,V6e]=ecgi_ECG_sim_plot_mat(ecg,V3e_f,l);
        
          %  full_ecg=ecgs{k};
          full_ecg=ecg.vars;
        full_ecg={full_ecg(:,1),full_ecg(:,2),full_ecg(:,3),full_ecg(:,4),full_ecg(:,5),full_ecg(:,6),full_ecg(:,7),full_ecg(:,8),full_ecg(:,9),full_ecg(:,10),full_ecg(:,11),full_ecg(:,12)};
    % plot_ecgs_2x6(full_ecg',names(3).name,0,l,1:size(full_ecg{1},1));
     hold on;
    long_params=vertcat(long_params,names(k+2).name);
    leg_names{l}=string(long_params{l});
    l=l+1;

begin=5;
    finish=5.5;
    begin2=9.5;
    finish2=size(full_ecg{1},1)-20;
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

    filtered{i}=ifft(F2)/6.0;
    filtered{i}(numel(full_ecg_ecgi{1}))=0;
  
     %   tot_ecg{i}=tot_ecg{i}+filtered{i};
    
        i=i+1;

end
      plot_ecgs_2x6(filtered',names(3).name,0,l-1,1:size(filtered{1},1));

end


leg=legend({'mild ARVC control model','advanced ARVC control model'},'Interpreter','none','Fontsize',30)
leg=legend(long_params,'Fontsize',30);
set(leg,'Interpreter','none')



for k=1:12
    tot_ecg{k}=tot_ecg{k}./(size(highest_scoring,1)-1);
end
figure()
     plot_ecgs_2x6(tot_ecg',names(3).name,0,1,1:size(filtered{1},1));


 %plot biomarkers for entire pop
 order=[];
for i=1:size(names,1)
    idx=strfind(names(i).name,'6.8');
    if isempty(idx)==0
        order=[order,i-2];
    end
end
 order2=[];

for i=1:size(names,1)
    idx=strfind(names(i).name,'10.2');
    if isempty(idx)==0
        order2=[order2,i-2];
    end
end
order3=[];
for i=1:size(names,1)
    idx=strfind(names(i).name,'13.6');
    if isempty(idx)==0
        order3=[order3,i-2];
    end
end
order4=[];
for i=1:size(names,1)
    idx=strfind(names(i).name,'17.0');
    if isempty(idx)==0
        order4=[order4,i-2];
    end
end

scal=[0.4 0.4 0.4 0.4 0.6 0.6 0.6 0.6 0.8 0.8 0.8 0.8 1 1 1 1];
tot_scal=[];
for i=1:256/16
    tot_scal=horzcat(tot_scal,scal);
end

scal_rv=[1 0.8 0.6 0.4];
tot_scal_RV=[];
for i=1:256/4
    tot_scal_RV=horzcat(tot_scal_RV,scal_rv);
end

purk_models=[15 79 143 207];

%partial correlations below
parameters = ones(256,4);
parameters(1:64,1)=0.2;

parameters(65:128,1)=0.3;
parameters(129:192,1)=0.4;
parameters(192:end,1)=.1;
for i=0:3
    parameters((i*64)+17:17+(i*64)+15,2)=0.6;
    parameters((i*64)+33:33+(i*64)+15,2)=0.5;
    parameters((i*64)+49:49+(i*64)+15,2)=0.4;


end

parameters(:,3)=tot_scal';
parameters(:,4)=tot_scal_RV';


diground=@(x,d) round(x*10^d)/10^d;

% r_s_ratio_value(find(r_s_ratio_value>100))=10;
% 
% r_s_ratio_full=QRSs./duration;

rs_ratio_mean= mean(ratio(1:4,:));

rs_ratio_mean(rs_ratio_mean>10)=10;

y=[averaged_QRSs' ,averaged_S_wave',rs_ratio_mean'];

y_full=[QRSs' ,duration'];

parameters_full=[parameters];

selected = [13,77,141,205];

y=[averaged_QRSs(selected)' ,averaged_S_wave(selected)'];
y_full=[QRSs(selected)' ,duration(selected)'];

parameters=parameters(selected,:);

parameters_full=[parameters];


%use another corr fcn

rho=partialcorri(parameters,y);

rho_full=partialcorri(parameters_full,y_full,'Rows','complete','Type','Pearson');

clf
figure()
heatmap(rho',{'Purkinje CV','LGE regional CV scaling','Bulk tissue CV scaling','RV CV scaling'},{'QRS duration','TAD','R/S-wave amplitude ratio'},'%0.2f','TickAngle',10,'FontSize',30,'Colormap','money')
set(gca,'FontSize',30)



figure()
h=heatmap(rho_full',{'Purkinje CV','Fibrosis Density','Bulk tissue CV','RV CV'},{'QRS V1','QRS V2','QRS V3','QRS V4','QRS V5','QRS V6','dtS V1','dtS V2','dtS V3','dtS V4','dtS V5','dtS V6','R-S rat V1','R-S rat V2','R-S rat V3','R-S rat V4','R-S rat V5','R-S rat V6'},'%0.2f','FontSize',20,'Colormap','money','ShowAllTicks',true)
set(gca,'FontSize',25)


%change in QRS duration wrt baseline with increasing Purkinje CV modulation
% inds={};
% opts= [0.4,0.6,0.8,1];
% qrs_durations=[];
% 
% 
% inds = [11 ,5 ,9,13]+130;
% for i=1:4
%    
%     qrs_durations(i)=   mean(QRSs(:,inds(i)));
%     tads(i)=mean(duration(1:4,inds(i)));
%     rs(i)=rs_ratio_mean(inds(i));
%     
% end
% qrs_durations_change = qrs_durations(4)-qrs_durations(1:3);
% tads_change = tads(4)-tads(1:3);
% rs_change = (rs(4)-rs(1:3))/rs(4) *100;
% 
% xlabels={'40% global CV slowing','50% global CV slowing','60% global CV slowing'};
% figure()
% bar(vertcat(qrs_durations_change,tads_change)',1)
% xticks([1.5 ,2.5 ,3.5])
% xticklabels(xlabels);
% xtickangle(15)
% ylabel(' change in feature value [ms]')
% set(gca,'FontSize',35)
% legend({'QRS duration','TAD'},'Fontsize',35)
% 
% %rs change % 
% figure()
% bar(vertcat(rs_change)',0.4,'r')
% xticks([1.5 ,2.5 ,3.5])
% xticklabels(xlabels);
% ylabel(' % change in R-S ratio')
% xtickangle(15)
% set(gca,'FontSize',35)
% legend({'R-S wave amplitude ratio'},'Fontsize',35)
% xlim([0.5,3.5])





%change in QRS duration wrt baseline with global bulk speed
inds={};
opts= [0.4,0.6,0.8,1];
qrs_durations=[];


inds = [1 ,5 ,9,13]+130; %increasing CV ?, yes 4 is control
  

j=4;

for i=1:4
    qrs_durations(j)=   mean(QRSs(:,inds(i)));
    tads(j)=mean(duration(1:4,inds(i)));
    rs(j)=rs_ratio_mean(inds(i));
    j=j-1;
    
end
qrs_durations_change = qrs_durations(2:4)-qrs_durations(1);
tads_change = tads(2:4)-tads(1);
rs_change = (rs(2:4)-rs(1))/rs(1) *100;


xlabels={'40% global CV slowing','50% global CV slowing','60% global CV slowing'};
figure()
bar(vertcat(qrs_durations_change,tads_change)',1)
xticks([1.5 ,2.5 ,3.5])
xticklabels(xlabels);
xtickangle(15)
ylabel(' change in feature value [ms]')
set(gca,'FontSize',35)
legend({'QRS duration','TAD'},'Fontsize',35)

%rs change % 
figure()
bar(vertcat(rs_change)',0.5,'r')
xticks([1.5 ,2.5 ,3.5])
xticklabels(xlabels);
ylabel(' % change in R-S ratio')
xtickangle(15)
set(gca,'FontSize',35)
legend({'R-S wave amplitude ratio'},'Fontsize',35)
xlim([0.5,3.5])




figure()
for i=1:4
hold on
plot(ecgs{inds(i)}(:,9))
end
legend




%plot s wave duration vs qrs duration on scatter line
%params: purk, fib,global scal, RV scal
figure()
hold on
scatter(QRSs(3,order),duration(3,order),155,'r','filled')
scatter(QRSs(3,order2),duration(3,order2),155,'g','filled')
scatter(QRSs(3,order3),duration(3,order3),155,'b','filled')
scatter(QRSs(3,order4),duration(3,order4),155,'m','filled')

set(gca,'FontSize',45)
xlabel('QRS duration [ms]','FontSize',45)
ylabel('S wave upstroke duration [ms]','FontSize',45)
legend({'40% Purkinje Speed','60% Purkinje Speed','80% Purkinje Speed','Normal Purkinje Speed'},'Fontsize',40)
xlim([20 150])
ylim([20 80])

figure()
bar([duration(1:4,:)']')
ylabel('Terminal activation duration [ms]','FontSize',30)
set(gca,'xtickLabel',{'V1','V2','V3','V4'})
%leg=legend({string(long_params{1}(1:end-53)),string(long_params{2}(1:end-53)),string(long_params{3}(1:end-53)),string(long_params{4}(1:end-53)),string(long_params{5}(1:end-53)),'clinical recordings'})
set(gca,'FontSize',40)
%set(leg,'Interpreter','none')
%legend(labs3)
     
     
figure()
bar([vertcat(QRSs(1:6,order-2)')]')
ylabel('QRS duration [ms]','FontSize',45)
set(gca,'xtickLabel',{'V1','V2','V3','V4','V5','V6'})
%leg=legend({string(long_params{1}(1:end-53)),string(long_params{2}(1:end-53)),string(long_params{3}(1:end-53)),string(long_params{4}(1:end-53)),string(long_params{5}(1:end-53)),'clinical recordings'})
set(gca,'FontSize',45)
%set(leg,'Interpreter','none')
%legend(labs3)
     



figure()
bar([duration(1:4,:)']')
ylabel('S wave upstroke duration [ms]','FontSize',45)
set(gca,'xtickLabel',{'V1','V2','V3','V4'})
%leg=legend({string(long_params{1}(1:end-53)),string(long_params{2}(1:end-53)),string(long_params{3}(1:end-53)),string(long_params{4}(1:end-53)),string(long_params{5}(1:end-53)),'clinical recordings'})
set(gca,'FontSize',45)
%set(leg,'Interpreter','none')
%legend({'Slow Purkinje Speed','Intermediate Purkinje Speed','Fast Purkinje Speed'},'Fontsize',25,'Location','northeastoutside')


figure()
bar([vertcat(QRSs(1:6,:)')]')
ylabel('QRS duration [ms]','FontSize',45)
set(gca,'xtickLabel',{'V1','V2','V3','V4','V5','V6'})
%leg=legend({string(long_params{1}(1:end-53)),string(long_params{2}(1:end-53)),string(long_params{3}(1:end-53)),string(long_params{4}(1:end-53)),string(long_params{5}(1:end-53)),'clinical recordings'})
set(gca,'FontSize',45)
%set(leg,'Interpreter','none')
legend(labs3)




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

sim_path='D:\ARVC meshing automatic\patients\patient06\results\Eikonal_coarse\figure_acts2\';

folders_acts=dir(sim_path);
for i = [1]%3:nb

    names_acts=dir(strcat(sim_path,folders_acts(i+2).name));
    upstroke{i}=load(strcat(sim_path,folders_acts(i+2).name,'\',names_acts(3).name));
    
end

if plt_acts==1
 
H_mesh = read_CHASTE( strcat(patient_sim_folder,    'chaste',patient_nr,'\HEART'));
H_mesh2 = read_CHASTE( strcat(patient_sim_folder,    'chaste06_coarse','\HEART'));
H_mesh3 = read_CHASTE( strcat(patient_sim_folder,    'chaste06_edge_0.15_RV_0.2','\HEART'));



folders_acts=dir(sim_path);
o=1;
for i =3%3:7
    names_acts=dir(strcat(sim_path,folders_acts(i).name));
    upstroke{o}=load(strcat(sim_path,folders_acts(i).name,'\',names_acts(3).name));
    o=o+1;
    
end


 figure()
    headlight
    hold on
    patch('vertices',H_mesh2.xyz,'faces',H_mesh2.face,'edgecolor','none','FaceVertexCData',upstroke{1},'facecolor','interp')
    colorbar
    caxis([0 200])
    axis off
    axis equal
       colormap(flipud(jet))

end
pat_nr='06';
electrodes = load(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\mpp\ECG_ELECTRODES.mat'));

E = [ electrodes.ELECTRODES.LA ; electrodes.ELECTRODES.RA ; electrodes.ELECTRODES.LL ; electrodes.ELECTRODES.RL ; electrodes.ELECTRODES.V1 ; electrodes.ELECTRODES.V2 ; electrodes.ELECTRODES.V3 ; electrodes.ELECTRODES.V4 ; electrodes.ELECTRODES.V5 ; electrodes.ELECTRODES.V6 ];
E=E/10; %chaste is in cm
hplot3d( E(5:7,:) , '*m' ,'LineWidth',4,'MarkerSize',10 );




H_mesh2.xyzacts=upstroke{1};
write_VTK(H_mesh2,strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\results\',sim_name,'\','control_eikonal','.vtk'),'binary');
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
    caxis([43 46])



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



 figure()
subplot_tight(2,1,1)
    hold on
    patch('vertices',H_mesh3.xyz*10,'faces',H_mesh3.face,'edgecolor','none','FaceVertexCData',[upstroke{i}],'facecolor','interp')
    [h,v] = IsoLine({ H_mesh3.face, H_mesh3.xyz*10},[upstroke{i}]',5,'k');
    %catter3(E(5,1),E(5,2),E(5,3),105,'m','filled')
    axis equal
    axis off
    caxis([5 75])
    colormap(flipud(jet))
        headlight
        view([29, 32])
subplot_tight(2,1,2)

    hold on
      plot_ecgi_rest('D:\ChrisReconstructions\ARVC6\ARVC6_Baseline.mat')
    axis equal
    caxis([38 65])
            view([29, 32])
    colormap(flipud(jet))
        headlight

 figure()
subplot_tight(2,1,1)
    hold on
    patch('vertices',H_mesh3.xyz*10,'faces',H_mesh3.face,'edgecolor','none','FaceVertexCData',[upstroke{i}],'facecolor','interp')
    [h,v] = IsoLine({ H_mesh3.face, H_mesh3.xyz*10},[upstroke{i}]',5,'k');
    %catter3(E(5,1),E(5,2),E(5,3),105,'m','filled')
    axis equal
    axis off
    caxis([20 30])
    colormap(flipud(jet))
        headlight
        view([29, 32])
subplot_tight(2,1,2)

    hold on
      plot_ecgi_rest('D:\ChrisReconstructions\ARVC6\ARVC6_Baseline.mat')
    axis equal
    caxis([43 44.5])
            view([29, 32])
    colormap(flipud(jet))
        headlight
        
        
        
 figure()
subplot_tight(2,1,1)
    hold on
    patch('vertices',H_mesh3.xyz*10,'faces',H_mesh3.face,'edgecolor','none','FaceVertexCData',[upstroke{i}],'facecolor','interp')
    [h,v] = IsoLine({ H_mesh3.face, H_mesh3.xyz*10},[upstroke{i}]',5,'k');
    %catter3(E(5,1),E(5,2),E(5,3),105,'m','filled')
    axis equal
    axis off
    caxis([20 30])
    colormap(flipud(jet))
        headlight
            view([-109, -60])
subplot_tight(2,1,2)

    hold on
      plot_ecgi_rest('D:\ChrisReconstructions\ARVC6\ARVC6_Baseline.mat')
    axis equal
    caxis([42 42.5])
            view([-109, -60])
    colormap(flipud(jet))
        headlight



            fid=fopen(strcat('D:\ARVC meshing automatic\patients\patient06','\results\TestHeart06NewRootsOldFibers\','results_maps.h5'),'rb');
           acts_bid = fread(fid,'double');
            permutation=load(strcat('D:\ARVC meshing automatic\patients\patient06','\results\TestHeart06NewRootsOldFibers\','permutation.txt'));
            permutation = permutation+1;

            acts_bid = acts_bid(permutation); %select only the heart activation times!! 
            acts_bid(find(acts_bid>200))=200;
            
 figure()
    hold on
    patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',acts_bid,'facecolor','interp')
    [h,v] = IsoLine({ H_mesh.face, H_mesh.xyz},[acts_bid]',5,'k');
    headlight
    
    fast_endo_activation=[47,37];
    slow_endo_activation=[111,89];