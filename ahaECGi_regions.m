close all
%aha style plots
addpath('/dataDropPeter/Dropbox/mesh scripts/AHABullseye/');
addpath(genpath('/data/DropPeter/Dropbox/mesh scripts/'));
%%%%%%%%%%%%%
filename = '/data/ChrisReconstructions/baselines/ARVC09_Baseline.mat';
regions_filename = '/home/scratch/meshes ready/patient02/regpatient aha ECGi.mat';
%%%%%%%%%%%%%
ecgi = load(filename);
pat_name = fields(ecgi);
regions = load(regions_filename);
regions = regions.regions_complete;
region_idx = [];
for i=1:26
    region_idx = vertcat(region_idx,regions{1,i}(:,1));
end

[GadLV,GadRV] =   read_csv('/data/sim3d/late_gad.csv',1);

GadRV=str2double(GadRV);
GadLV=str2double(GadLV);

for j=[7 18]%[3,4,11,12,15,19]%[1,5,6,7,8,17]%[2,5,9,10,14,18]
    
patient_nr=j;
for i=1:20
    pat_LGE_LV(i)=GadLV(i+20*(patient_nr-1));
        pat_LGE_RV(i)=GadRV(i+20*(patient_nr-1));

end
%  Sample Bullseye Creation
pat_LGE_RV=pat_LGE_RV(11:end)*100/3;
% Create the AHA 17-segment bullseye
figure;
c = createBullseye([0 0.5 2 0; 0.5 1 6 0; 1 1.5 6 0]);
set(c,'Color','k','LineWidth',4)
set(gca,'Fontsize',22)
title(strcat('LV LGE pat nr ' ,num2str(j)));

% Example 1 of filling the bullseye, vector by vector
fillBullseye(pat_LGE_LV(9),0.01,0.49,0.5,179.5); %apex
fillBullseye(pat_LGE_LV(10),0.01,0.49,180.5,359.5); %apex


fillBullseye([pat_LGE_LV(4)],0.51,0.99,0.5,59.5); %lv mid posterolateral wall
fillBullseye([pat_LGE_LV(2)],0.51,0.99,60.5,119.5); %lv mid posterior wall
fillBullseye(nan,0.51,0.99,120.5,179.5); %septal posterior
fillBullseye(nan,0.51,0.99,180.5,239.5); %septal anterior
fillBullseye(pat_LGE_LV(8),0.51,0.99,240.5,299.5); %lv mid  anterior
fillBullseye(pat_LGE_LV(6),0.51,0.99,300.5,359.5); %lv mid  anterolateral

fillBullseye([100],1.01,1.49,0.5,59.5); %lv basal posterolateral wall
fillBullseye([pat_LGE_LV(3)],1.01,1.49,0.5,59.5); %lv basal posterolateral wall
fillBullseye([pat_LGE_LV(1)],1.01,1.49,60.5,119.5); %lv basal posterior wall
fillBullseye(nan,1.01,1.49,120.5,179.5); %septal basal posterior
fillBullseye(nan,1.01,1.49,180.5,239.5); %septal basal anterior
fillBullseye(pat_LGE_LV(7),1.01,1.49,240.5,299.5); %lv basal  anterior
fillBullseye(pat_LGE_LV(5),1.01,1.49,300.5,359.5); %lv basal  anterolateral
uistack(c,'top');
c=colorbar()
c.Label.String='LGE MRI';
c.Limits=[0,100];
axis off

% Create the AHA 17-segment bullseye
figure;
title(strcat('RV LGE pat nr ' ,num2str(j)));

c = createBullseye([0 0.5 2 0; 0.5 1 6 0; 1 1.5 6 0]);
set(c,'Color','k','LineWidth',4)
set(gca,'Fontsize',22)
%pat_LGE_RV=zeros(1,10)+0.00000001;

% Example 1 of filling the bullseye, vector by vector
fillBullseye(pat_LGE_RV(9),0.01,0.49,0.5,179.5); %apex
fillBullseye(pat_LGE_RV(10),0.01,0.49,180.5,359.5); %apex


fillBullseye(nan,0.51,0.99,0.5,59.5); %lv mid posterolateral wall
fillBullseye([pat_LGE_RV(2)],0.51,0.99,60.5,119.5); %lv mid posterior wall
fillBullseye(pat_LGE_RV(4),0.51,0.99,120.5,179.5); %septal posterior
fillBullseye(pat_LGE_RV(6),0.51,0.99,180.5,239.5); %septal anterior
fillBullseye(pat_LGE_RV(8),0.51,0.99,240.5,299.5); %lv mid  anterior
fillBullseye(nan,0.51,0.99,300.5,359.5); %lv mid  anterolateral


fillBullseye(nan,1.01,1.49,0.5,59.5); %lv basal posterolateral wall
fillBullseye([pat_LGE_RV(1)],1.01,1.49,60.5,119.5); %lv basal posterior wall
fillBullseye(pat_LGE_RV(3),1.01,1.49,120.5,179.5); %septal basal posterior
fillBullseye(pat_LGE_RV(5),1.01,1.49,180.5,239.5); %septal basal anterior
fillBullseye(pat_LGE_RV(7),1.01,1.49,240.5,299.5); %lv basal  anterior
fillBullseye(nan,1.01,1.49,300.5,359.5); %lv basal  anterolateral
uistack(c,'top');
col=colorbar()
col.Label.String='LGE MRI';
col.Limits=[0,100]
axis off
axis equal


end





figure;
title('Time of activation of last site in region')
c = createBullseye([0 0.5 1 0; 0.5 1 4 45; 1 1.5 6 0; 1.5 2 6 0]);
set(c,'Color','k','LineWidth',2)
set(gca,'Fontsize',22)

% Example 1 of filling the bullseye, vector by vector
fillBullseye(lv_summary{12,3},0,0.5,-45,315); %apex

fillBullseye([lv_summary{11,3}],0.5,1,45,135); %lv apical anterior
fillBullseye([lv_summary{10,3}],0.5,1,135,225); %lv apical free wall
fillBullseye([lv_summary{9,3}],0.5,1,-135,-45);%lv apical inferior

fillBullseye([lv_summary{7,3}],1,1.5,120,180); %lv mid anterolateral wall
fillBullseye([lv_summary{8,3}],1,1.5,60,120); %lv mid anterior
fillBullseye([lv_summary{5,3}],1,1.5,-180,-120); %lv mid inferior
fillBullseye([lv_summary{6,3}],1,1.5,-120,-60); %lv mid inferolateral wall 

fillBullseye([lv_summary{3,3}],1.5,2,120,180); %lv basal anterolateral wall
fillBullseye([lv_summary{4,3}],1.5,2,60,120); %lv basal anterior
fillBullseye([lv_summary{1,3}],1.5,2,-180,-120); %lv basal inferior
fillBullseye([lv_summary{2,3}],1.5,2,-120,-60); %lv basal inferolateral wall 


uistack(c,'top');
col=colorbar()
col.Label.String='activation time [ms]';



figure;
title('Maximum regional activation gradient')
c = createBullseye([0 0.5 1 0; 0.5 1 4 45; 1 1.5 6 0; 1.5 2 6 0]);
set(c,'Color','k','LineWidth',2)
set(gca,'Fontsize',22)

% Example 1 of filling the bullseye, vector by vector
fillBullseye(lv_summary{12,4},0,0.5,-45,315); %apex

fillBullseye([lv_summary{11,4}],0.5,1,45,135); %lv apical anterior
fillBullseye([lv_summary{10,4}],0.5,1,135,225); %lv apical free wall
fillBullseye([lv_summary{9,4}],0.5,1,-135,-45);%lv apical inferior

fillBullseye([lv_summary{7,4}],1,1.5,120,180); %lv mid anterolateral wall
fillBullseye([lv_summary{8,4}],1,1.5,60,120); %lv mid anterior
fillBullseye([lv_summary{5,4}],1,1.5,-180,-120); %lv mid inferior
fillBullseye([lv_summary{6,4}],1,1.5,-120,-60); %lv mid inferolateral wall 

fillBullseye([lv_summary{3,4}],1.5,2,120,180); %lv basal anterolateral wall
fillBullseye([lv_summary{4,4}],1.5,2,60,120); %lv basal anterior
fillBullseye([lv_summary{1,4}],1.5,2,-180,-120); %lv basal inferior
fillBullseye([lv_summary{2,4}],1.5,2,-120,-60); %lv basal inferolateral wall 


uistack(c,'top');
c=colorbar()
c.Label.String='activation gradient [ms/mm]';



figure;
title('Mean regional activation gradient')
c = createBullseye([0 0.5 1 0; 0.5 1 4 45; 1 1.5 6 0; 1.5 2 6 0]);
set(c,'Color','k','LineWidth',2)
set(gca,'Fontsize',22)

% Example 1 of filling the bullseye, vector by vector
fillBullseye(lv_summary{12,6},0,0.5,-45,315); %apex

fillBullseye([lv_summary{11,6}],0.5,1,45,135); %lv apical anterior
fillBullseye([lv_summary{10,6}],0.5,1,135,225); %lv apical free wall
fillBullseye([lv_summary{9,6}],0.5,1,-135,-45);%lv apical inferior

fillBullseye([lv_summary{7,6}],1,1.5,120,180); %lv mid anterolateral wall
fillBullseye([lv_summary{8,6}],1,1.5,60,120); %lv mid anterior
fillBullseye([lv_summary{5,6}],1,1.5,-180,-120); %lv mid inferior
fillBullseye([lv_summary{6,6}],1,1.5,-120,-60); %lv mid inferolateral wall 

fillBullseye([lv_summary{3,6}],1.5,2,120,180); %lv basal anterolateral wall
fillBullseye([lv_summary{4,6}],1.5,2,60,120); %lv basal anterior
fillBullseye([lv_summary{1,6}],1.5,2,-180,-120); %lv basal inferior
fillBullseye([lv_summary{2,6}],1.5,2,-120,-60); %lv basal inferolateral wall 

uistack(c,'top');
c=colorbar()
c.Label.String='activation gradient [ms/mm]';




figure;
title('Mean regional activation time')
c = createBullseye([0 0.5 1 0; 0.5 1 4 45; 1 1.5 6 0; 1.5 2 6 0]);
set(c,'Color','k','LineWidth',2)
set(gca,'Fontsize',22)

% Example 1 of filling the bullseye, vector by vector
fillBullseye(lv_summary{12,2},0,0.5,-45,315); %apex

fillBullseye([lv_summary{11,5}],0.5,1,45,135); %lv apical anterior
fillBullseye([lv_summary{10,5}],0.5,1,135,225); %lv apical free wall
fillBullseye([lv_summary{9,5}],0.5,1,-135,-45);%lv apical inferior

fillBullseye([lv_summary{7,5}],1,1.5,120,180); %lv mid anterolateral wall
fillBullseye([lv_summary{8,5}],1,1.5,60,120); %lv mid anterior
fillBullseye([lv_summary{5,5}],1,1.5,-180,-120); %lv mid inferior
fillBullseye([lv_summary{6,5}],1,1.5,-120,-60); %lv mid inferolateral wall 

fillBullseye([lv_summary{3,5}],1.5,2,120,180); %lv basal anterolateral wall
fillBullseye([lv_summary{4,5}],1.5,2,60,120); %lv basal anterior
fillBullseye([lv_summary{1,5}],1.5,2,-180,-120); %lv basal inferior
fillBullseye([lv_summary{2,5}],1.5,2,-120,-60); %lv basal inferolateral wall 


uistack(c,'top');
c=colorbar()
c.Label.String='activation time [ms]';
