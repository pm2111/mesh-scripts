function [LV_summary]=rv_stats(mesh_location,aha_regions,patient_nr)
addpath('/users/petnov/Dropbox/mesh scripts/');
filenames = mesh_location;
data = load(mesh_location);


regions = load(aha_regions);
regions = regions.regions_complete;

nr_lv_regions = 13;
nr_rv_regions = size(regions,2)-nr_lv_regions;
LV_summary = cell(1);
name_dummy = [];
for i=1
    fields = fieldnames(data);
    activation = data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart;
    
    %find 1st 4 breakthroughs
    %do a regional based analysis
    region_minima_ind = [];
    region_minima_value = [];
    for j=nr_lv_regions+1:size(regions,2)
       regions{1,j}(:,1) = regions{1,j}(:,1);
       activation_region = zeros(1002,1)+10000;
       activation_region(regions{1,j}(:,1))= activation( regions{1,j}(:,1));
       
       LV_summary{j-nr_lv_regions,5} = mean(activation_region);
       [value,id] = min(activation_region);
       region_minima_ind(j-nr_lv_regions)=id;
       region_minima_value(j-nr_lv_regions) = value;
       LV_summary{j-nr_lv_regions,1} = j;
       LV_summary{j-nr_lv_regions,2} = value; %minima of activation
    end
    


    %find top 3 earliest activations
    region_dummy = region_minima_value;
    earliest_activation_idx = [];
    earliest_activation_value = [];
    for j=1:4
        [value,idx] = min(region_dummy);
        region_idx(j) = idx;
        earliest_activation_idx(j)= region_minima_ind(idx);
        earliest_activation_value(j) = value;
        region_dummy(idx) = 1000;
        name_dummy = [name_dummy,idx]; 
    end
end



%now find the maxima of activation, namely the regions of slow conduction
for i=1
    fields = fieldnames(data);
    activation = data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart;
    
    %find 1st 4 slow regions
    %do a regional based analysis
    region_minima_ind = [];
    region_minima_value = [];
    for j=nr_lv_regions+1:size(regions,2)
       activation_region = zeros(1002,1);
       activation_region(regions{1,j}(:,1))= activation( regions{1,j}(:,1));
       [value,id] = max(activation_region);
       region_minima_ind(j-nr_lv_regions)=id;
       region_minima_value(j-nr_lv_regions) = value;
       LV_summary{j-nr_lv_regions,3} = value; %maxima of activation
    end

    %find top 3 earliest activations
    region_dummy = region_minima_value;
    earliest_activation_idx = [];
    for j=1:4
        [value,idx] = max(region_dummy);
        region_idx_maxima(j) = idx;
        latest_activation_value(j) =value;
        earliest_activation_idx(j)= region_minima_ind(idx);
        region_dummy(idx) = 0;
        name_dummy = [name_dummy,idx]; 
    end
end
   

%gradients

region_grad_max = [];
filenames_grad = dir('/data/ecgi/gradients/');
for i=1
    
    grad = csvread(strcat('/data/ecgi/gradients/', filenames_grad(patient_nr+2).name),1);
    grad_mag = sqrt( power(grad(:,2),2) + power(grad(:,3),2) + power(grad(:,4),2)); %find the magnitude of the powers

    for j=nr_lv_regions+1:size(regions,2)
       %grad_region = zeros(1002,1)+10000;
       grad_region= grad_mag( regions{1,j}(:,1));
       [value,id] = max(grad_region);   
       region_grad_max(i,j) = value;
       LV_summary{j-nr_lv_regions,4} = value;
       LV_summary{j-nr_lv_regions,6} = mean(grad_region);

    end
end
% for i=1
%     [mx,idx] = max(region_grad_max(i,1:nr_lv_regions));
%     
%     LV_summary{i,5} = mx;
%     LV_summary{i,6} = idx+1;
% end
% 
% for i=1
%     LV_summary{i,2} = region_names([LV_summary{i,2}]);
%     LV_summary{i,4} = region_names([LV_summary{i,4}]);
%     LV_summary{i,6} = region_names([LV_summary{i,6}]);
% end


%summary table, LV
