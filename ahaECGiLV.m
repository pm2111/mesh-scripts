close all
clear all
%aha style plots
addpath('/users/petnov/Dropbox/mesh scripts/AHABullseye/');
addpath('/users/petnov/Dropbox/mesh scripts/');
%%%%%%%%%%%%%
filename = '/data/ChrisReconstructions/baselines/ARVC02_Baseline.mat';
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

 lv_summary =   lv_stats(filename,regions_filename,2);

%  Sample Bullseye Creation

% Create the AHA 17-segment bullseye
figure;
title('Time of activation of first site in region')
c = createBullseye([0 0.5 1 0; 0.5 1 4 45; 1 1.5 6 0; 1.5 2 6 0]);
set(c,'Color','k','LineWidth',2)
set(gca,'Fontsize',22)

% Example 1 of filling the bullseye, vector by vector
fillBullseye(lv_summary{12,2},0,0.5,-45,315); %apex

fillBullseye([lv_summary{11,2}],0.5,1,45,135); %lv apical anterior
fillBullseye([lv_summary{10,2}],0.5,1,-45,45); %lv apical free wall
fillBullseye([lv_summary{9,2}],0.5,1,-135,-45);%lv apical inferior

fillBullseye([lv_summary{7,2}],1,1.5,0,60); %lv mid anterolateral wall
fillBullseye([lv_summary{8,2}],1,1.5,60,120); %lv mid anterior
fillBullseye([lv_summary{5,2}],1,1.5,-60,0); %lv mid inferior
fillBullseye([lv_summary{6,2}],1,1.5,-120,-60); %lv mid inferolateral wall 

fillBullseye([lv_summary{3,2}],1.5,2,0,60); %lv basal anterolateral wall
fillBullseye([lv_summary{4,2}],1.5,2,60,120); %lv basal anterior
fillBullseye([lv_summary{1,2}],1.5,2,-60,0); %lv basal inferior
fillBullseye([lv_summary{2,2}],1.5,2,-120,-60); %lv basal inferolateral wall 

uistack(c,'top');
c=colorbar()
c.Label.String='activation time [ms]';


figure;
title('Time of activation of last site in region')
c = createBullseye([0 0.5 1 0; 0.5 1 4 45; 1 1.5 6 0; 1.5 2 6 0]);
set(c,'Color','k','LineWidth',2)
set(gca,'Fontsize',22)

% Example 1 of filling the bullseye, vector by vector
fillBullseye(lv_summary{12,3},0,0.5,-45,315); %apex

fillBullseye([lv_summary{11,3}],0.5,1,45,135); %lv apical anterior
fillBullseye([lv_summary{10,3}],0.5,1,-45,45); %lv apical free wall
fillBullseye([lv_summary{9,3}],0.5,1,-135,-45);%lv apical inferior

fillBullseye([lv_summary{7,3}],1,1.5,0,60); %lv mid anterolateral wall
fillBullseye([lv_summary{8,3}],1,1.5,60,120); %lv mid anterior
fillBullseye([lv_summary{5,3}],1,1.5,-60,0); %lv mid inferior
fillBullseye([lv_summary{6,3}],1,1.5,-120,-60); %lv mid inferolateral wall 

fillBullseye([lv_summary{3,3}],1.5,2,0,60); %lv basal anterolateral wall
fillBullseye([lv_summary{4,3}],1.5,2,60,120); %lv basal anterior
fillBullseye([lv_summary{1,3}],1.5,2,-60,0); %lv basal inferior
fillBullseye([lv_summary{2,3}],1.5,2,-120,-60); %lv basal inferolateral wall 

uistack(c,'top');
c=colorbar()
c.Label.String='activation time [ms]';



figure;
title('Maximum regional activation gradient')
c = createBullseye([0 0.5 1 0; 0.5 1 4 45; 1 1.5 6 0; 1.5 2 6 0]);
set(c,'Color','k','LineWidth',2)
set(gca,'Fontsize',22)

% Example 1 of filling the bullseye, vector by vector
fillBullseye(lv_summary{12,4},0,0.5,-45,315); %apex

fillBullseye([lv_summary{11,4}],0.5,1,45,135); %lv apical anterior
fillBullseye([lv_summary{10,4}],0.5,1,-45,45); %lv apical free wall
fillBullseye([lv_summary{9,4}],0.5,1,-135,-45);%lv apical inferior

fillBullseye([lv_summary{7,4}],1,1.5,0,60); %lv mid anterolateral wall
fillBullseye([lv_summary{8,4}],1,1.5,60,120); %lv mid anterior
fillBullseye([lv_summary{5,4}],1,1.5,-60,0); %lv mid inferior
fillBullseye([lv_summary{6,4}],1,1.5,-120,-60); %lv mid inferolateral wall 

fillBullseye([lv_summary{3,4}],1.5,2,0,60); %lv basal anterolateral wall
fillBullseye([lv_summary{4,4}],1.5,2,60,120); %lv basal anterior
fillBullseye([lv_summary{1,4}],1.5,2,-60,0); %lv basal inferior
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
fillBullseye([lv_summary{10,6}],0.5,1,-45,45); %lv apical free wall
fillBullseye([lv_summary{9,6}],0.5,1,-135,-45);%lv apical inferior

fillBullseye([lv_summary{7,6}],1,1.5,0,60); %lv mid anterolateral wall
fillBullseye([lv_summary{8,6}],1,1.5,60,120); %lv mid anterior
fillBullseye([lv_summary{5,6}],1,1.5,-60,0); %lv mid inferior
fillBullseye([lv_summary{6,6}],1,1.5,-120,-60); %lv mid inferolateral wall 

fillBullseye([lv_summary{3,6}],1.5,2,0,60); %lv basal anterolateral wall
fillBullseye([lv_summary{4,6}],1.5,2,60,120); %lv basal anterior
fillBullseye([lv_summary{1,6}],1.5,2,-60,0); %lv basal inferior
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
fillBullseye(lv_summary{12,5},0,0.5,-45,315); %apex

fillBullseye([lv_summary{11,5}],0.5,1,45,135); %lv apical anterior
fillBullseye([lv_summary{10,5}],0.5,1,-45,45); %lv apical free wall
fillBullseye([lv_summary{9,5}],0.5,1,-135,-45);%lv apical inferior

fillBullseye([lv_summary{7,5}],1,1.5,0,60); %lv mid anterolateral wall
fillBullseye([lv_summary{8,5}],1,1.5,60,120); %lv mid anterior
fillBullseye([lv_summary{5,5}],1,1.5,-60,0); %lv mid inferior
fillBullseye([lv_summary{6,5}],1,1.5,-120,-60); %lv mid inferolateral wall 

fillBullseye([lv_summary{3,5}],1.5,2,0,60); %lv basal anterolateral wall
fillBullseye([lv_summary{4,5}],1.5,2,60,120); %lv basal anterior
fillBullseye([lv_summary{1,5}],1.5,2,-60,0); %lv basal inferior
fillBullseye([lv_summary{2,5}],1.5,2,-120,-60); %lv basal inferolateral wall 

uistack(c,'top');
c=colorbar()
c.Label.String='activation time [ms]';
