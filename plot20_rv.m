cd /data/ecgi/
addpath('/users/petnov/Dropbox/mesh scripts/');
filenames = dir('/data/ChrisReconstructions/baselines/');

RV_summary = cell(1);
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
    for j=12:21
       activation_region = zeros(1002,1)+10000;
       activation_region(data.(fields{1}).labels.(region_names{j}))= activation( data.(fields{1}).labels.(region_names{j}));
       [value,id] = min(activation_region);
       region_minima_ind(j)=id;
       region_minima_value(j) = value;
    end
    
    region_minima_value(1:11) = 10000;

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
    RV_summary{i,1} =  earliest_activation_value(2); %1 gives zero because of j going from 2:11 (1st element is empty)
    RV_summary{i,2} = region_idx(2);
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
    for j=12:21
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
    RV_summary{i,3} =  latest_activation_value(2); %1 gives zero because of j going from 2:11 (1st element is empty)
    RV_summary{i,4} = region_idx_maxima(2);
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
    [mx,idx] = max(region_grad_max(i,12:21));
    
    RV_summary{i,5} = mx;
    RV_summary{i,6} = idx+11;
end

for i=1:20
    RV_summary{i,2} = region_names([RV_summary{i,2}]);
    RV_summary{i,4} = region_names([RV_summary{i,4}]);
    RV_summary{i,6} = region_names([RV_summary{i,6}]);
end

tot =0;
time_lv = 0;
time_rv = 0;
for i=1:20
    a = RV_summary{i,1} < LV_summary{i,1};
    time_lv = time_lv + LV_summary{i,3}- LV_summary{i,1};
    time_rv = time_rv + RV_summary{i,3} - RV_summary{i,1};
    tot = tot+a;
end

for i=1:20
    diff(i)=[RV_summary{i,1}-LV_summary{i,1}];
end
optimal_activation = [];
optimal_list=[];
for i=1:20
    optimal_list =vertcat(optimal_list,[RV_summary{i,1} - LV_summary{i,1} < 0 , RV_summary{i,5} <3 , LV_summary{i,5} <3 , RV_summary{i,3} < 60 , LV_summary{i,3} <60]);
    if RV_summary{i,1} - LV_summary{i,1} < 0 && RV_summary{i,5} <3 && LV_summary{i,5} <3 && RV_summary{i,3} < 60 && LV_summary{i,3} <60

        optimal_activation = [optimal_activation,i];
    end
end


sub_optimal_activation = [];
sub_optimal_list=[];
for i=1:20
     sub_optimal_list =vertcat(sub_optimal_list,[RV_summary{i,1} - LV_summary{i,1} < 10 , RV_summary{i,5} <4 , LV_summary{i,5} <4 , RV_summary{i,3} < 70 , LV_summary{i,3} <70 , ismember(i,optimal_activation)==0])
    if  RV_summary{i,1} - LV_summary{i,1} < 10 && RV_summary{i,5} <5 && LV_summary{i,5} <5 && RV_summary{i,3} < 70 && LV_summary{i,3} <70 && ismember(i,optimal_activation)==0
        
        sub_optimal_activation = [sub_optimal_activation,i];
    end
end

sub_sub_optimal_activation = [];
for i=1:20
    if  RV_summary{i,1} - LV_summary{i,1} < 20 && RV_summary{i,5} <7 && LV_summary{i,5} <7 && RV_summary{i,3} < 75 && LV_summary{i,3} <75 && ismember(i,optimal_activation)==0

        sub_sub_optimal_activation = [sub_sub_optimal_activation,i];
    end
end


sub_sub_sub_optimal_activation = [];

for i=1:20
    if  ismember(i,optimal_activation)==0 && ismember(i,sub_optimal_activation)==0 && ismember(i,sub_sub_optimal_activation)==0

        sub_sub_sub_optimal_activation = [sub_sub_sub_optimal_activation,i];
    end
end

%try some clustering on this dataa

X = [LV_summary{:,1};RV_summary{:,1} ;LV_summary{:,5}; RV_summary{:,5} ;LV_summary{:,3} ;RV_summary{:,3}]';
idx = kmeans(X,2);

%below try hierarchical clustering as k means doesnt produce good results
Y=pdist(X(:,1:2));
Z=linkage(Y);
%figures of category making

figure()
hold on
plot([RV_summary{optimal_activation,1}] -[ LV_summary{optimal_activation,1}],[RV_summary{optimal_activation,5}],'r*')
plot([RV_summary{sub_optimal_activation,1}] -[ LV_summary{sub_optimal_activation,1}],[RV_summary{sub_optimal_activation,5}],'m*')


figure()
hist([RV_summary{:,1}] - [LV_summary{:,1}])
xlabel('Difference between RV and LV Breakthrough times [ms]')
ylabel('Count')
title('Criteria 1: Time delay of RV epicardial breakthrough relative to LV breakthrough')
set(gca, 'FontSize', 20)

figure()
hist([RV_summary{:,3}] )
xlabel('RV complete depolarisation times [ms]')
ylabel('Count')
title('Criteria 2a: Time to completely depolarise the RV')
set(gca, 'FontSize', 20)

figure()
hist([LV_summary{:,3}] )
xlabel('LV complete depolarisation times [ms]')
ylabel('Count')
title('Criteria 2b: Time to completely depolarise the LV')
set(gca, 'FontSize', 20)


figure()
hist([RV_summary{:,5}] )
xlabel('RV maximum activation gradients [ms/cm]')
ylabel('Count')
title('Criteria 3a: Maximum activation gradients on the RV epicardial surface')
set(gca, 'FontSize', 20)

figure()
hist([LV_summary{:,5}] )
xlabel('LV maximum activation gradients [ms/cm]')
ylabel('Count')
title('Criteria 3b: Maximum activation gradients on the LV epicardial surface')
set(gca, 'FontSize', 20)

figure()
hold on
scatter(1:sum(idx==1),X(idx==1,1),'g','filled')
scatter(1:sum(idx==2),X(idx==2,1),'b','filled')
scatter(1:sum(idx==3),X(idx==3,1),'r','filled')


figure()
hold on
scatter(1:sum(idx==1),X(idx==1,2),'g','filled')
scatter(1:sum(idx==2),X(idx==2,2),'b','filled')
scatter(1:sum(idx==3),X(idx==3,2),'r','filled')


figure()
hold on
scatter(1:sum(idx==1),X(idx==1,3),'g','filled')
scatter(1:sum(idx==2),X(idx==2,3),'b','filled')
scatter(1:sum(idx==3),X(idx==3,3),'r','filled')

figure()
hold on
scatter(1:sum(idx==1),X(idx==1,4),'g','filled')
scatter(1:sum(idx==2),X(idx==2,4),'b','filled')
scatter(1:sum(idx==3),X(idx==3,4),'r','filled')


%summary table, LV
