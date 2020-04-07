cd /data/ecgi/
addpath('/users/petnov/Dropbox/mesh scripts/');
filenames = dir('/data/ChrisReconstructions/baselines/');

LV_summary = cell(1);
name_dummy = [];
for i=1:size(filenames,1)-2
    data = load(strcat(filenames(i+2).folder,'/',filenames(i+2).name));
    fields = fieldnames(data);
    region_names = fieldnames(data.(fields{1}).labels);
    activation = data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart;
    
    %find 1st 4 breakthroughs
    %do a regional based analysis
    region_names=fieldnames(data.(fields{1}).labels);
    region_minima_ind = [];
    region_minima_value = [];
    for j=2:11
       activation_region = zeros(1002,1)+10000;
       activation_region(data.(fields{1}).labels.(region_names{j}))= activation( data.(fields{1}).labels.(region_names{j}));
       [value,id] = min(activation_region);
       region_minima_ind(j)=id;
       region_minima_value(j) = value;
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
    LV_summary{i,1} =  earliest_activation_value(2); %1 gives zero because of j going from 2:11 (1st element is empty)
    LV_summary{i,2} = region_idx(2);
end



%now find the maxima of activation, namely the regions of slow conduction
for i=1:size(filenames,1)-2
    data = load(strcat(filenames(i+2).folder,'/',filenames(i+2).name));
    fields = fieldnames(data);
    region_names = fieldnames(data.(fields{1}).labels);
    activation = data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart;
    
    %find 1st 4 slow regions
    %do a regional based analysis
    region_names=fieldnames(data.(fields{1}).labels);
    region_minima_ind = [];
    region_minima_value = [];
    for j=2:11
       activation_region = zeros(1002,1);
       activation_region(data.(fields{1}).labels.(region_names{j}))= activation( data.(fields{1}).labels.(region_names{j}));
       [value,id] = max(activation_region);
       region_minima_ind(j)=id;
       region_minima_value(j) = value;
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
    LV_summary{i,3} =  latest_activation_value(2); %1 gives zero because of j going from 2:11 (1st element is empty)
    LV_summary{i,4} = region_idx_maxima(2);
end
   

%gradients

region_grad_max = [];
filenames_grad = dir('/data/ecgi/gradients/');
for i=1:size(filenames_grad,1)-2
    
    grad = csvread(strcat('/data/ecgi/gradients/', filenames_grad(i+2).name),1);
    grad_mag = sqrt( power(grad(:,2),2) + power(grad(:,3),2) + power(grad(:,4),2)); %find the magnitude of the powers

    for j=1:21
       %grad_region = zeros(1002,1)+10000;
       grad_region= grad_mag( data.(fields{1}).labels.(region_names{j}));
       [value,id] = max(grad_region);   
       region_grad_max(i,j) = value;
    end
end
for i=1:20
    [mx,idx] = max(region_grad_max(i,2:11));
    
    LV_summary{i,5} = mx;
    LV_summary{i,6} = idx+1;
end

for i=1:20
    LV_summary{i,2} = region_names([LV_summary{i,2}]);
    LV_summary{i,4} = region_names([LV_summary{i,4}]);
    LV_summary{i,6} = region_names([LV_summary{i,6}]);
end


%summary table, LV
