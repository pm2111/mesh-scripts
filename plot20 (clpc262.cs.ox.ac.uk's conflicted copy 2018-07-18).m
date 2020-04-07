cd /data/ecgi/
addpath('/users/petnov/Dropbox/mesh scripts/');
filenames = dir('/data/ChrisReconstructions/baselines/');

LV_summary = cell(1);

figure()
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
    for j=1:21
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
    for j=1:5
        [value,idx] = min(region_dummy);
        earliest_activation_idx(j)= region_minima_ind(idx);
        earliest_activation_value(j) = value;
        region_dummy(idx) = 1000;
        name_dummy = [name_dummy,idx]; 
    end
    
    subplot('Position',[(mod(i-1,5))/5 0.9-(ceil(i/5))/5 1/5 1/5])
    
    %patch('vertices', data.arvc2_baseline.nodes, 'faces', data.arvc2_baseline.mesh, 'facecolor', 'interp', 'cdata', labels(:,1), 'edgecolor', 'none') ;
    p = patch('vertices', data.(fields{1}).nodes, 'faces', data.(fields{1}).mesh,'facecolor', 'interp', 'cdata', data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart,'edgecolor', 'none') ;
    %p.Facecolor = 'Flat';
    cb =colorbar;
    caxis([0 70])
    axis equal
    %view([50,75]) %front view
    view([225,-90]) %back view
    %%view([125,20]) % lv side view
    %view([-35,-20]) % rv side view
    axis off
    set(cb,'position',[0.95 .35 .03 .3])
    %title(strcat(filepath(end-18:end-13),' activation times [ms] lv view'))
    title(filenames(i+2).name(end-18:end-13))
    [h,v] = IsoLine({data.(fields{1}).mesh, data.(fields{1}).nodes},data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart,10,'k');
    plot3(data.(fields{1}).nodes(earliest_activation_idx,1),data.(fields{1}).nodes(earliest_activation_idx,2),data.(fields{1}).nodes(earliest_activation_idx,3),'r*','markers',16);

end

figure()

hist(name_dummy,21);
xticks( linspace(1,21,21)-0.5);
xticklabels( (region_names));
xt = get(gca, 'XTick');
set(gca, 'FontSize', 10)


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
    for j=1:21
       activation_region = zeros(1002,1);
       activation_region(data.(fields{1}).labels.(region_names{j}))= activation( data.(fields{1}).labels.(region_names{j}));
       [value,id] = max(activation_region);
       region_minima_ind(j)=id;
       region_minima_value(j) = value;
    end

    %find top 5 earliest activations
    region_dummy = region_minima_value;
    earliest_activation_idx = [];
    for j=1:5
        [value,idx] = max(region_dummy);
        earliest_activation_idx(j)= region_minima_ind(idx);
        region_dummy(idx) = 0;
        name_dummy = [name_dummy,idx]; 
    end

    
    subplot('Position',[(mod(i-1,5))/5 0.9-(ceil(i/5))/5 1/5 1/5])
    
    %patch('vertices', data.arvc2_baseline.nodes, 'faces', data.arvc2_baseline.mesh, 'facecolor', 'interp', 'cdata', labels(:,1), 'edgecolor', 'none') ;
    p = patch('vertices', data.(fields{1}).nodes, 'faces', data.(fields{1}).mesh,'facecolor', 'interp', 'cdata', data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart,'edgecolor', 'none') ;
    %p.Facecolor = 'Flat';
    cb =colorbar;
    caxis([0 70])
    axis equal
    %view([50,75]) %front view
    view([225,-90]) %back view
    %%view([125,20]) % lv side view
    %view([-35,-20]) % rv side view
    axis off
    set(cb,'position',[0.95 .35 .03 .3])
    %title(strcat(filepath(end-18:end-13),' activation times [ms] lv view'))
    title(filenames(i+2).name(end-18:end-13))
    [h,v] = IsoLine({data.(fields{1}).mesh, data.(fields{1}).nodes},data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart,10,'k');
    plot3(data.(fields{1}).nodes(earliest_activation_idx,1),data.(fields{1}).nodes(earliest_activation_idx,2),data.(fields{1}).nodes(earliest_activation_idx,3),'r*','markers',16);

end
%find the mean activation per region
patients_mean_activation = [];
patients_median_activation = [];
patients_min_activation = [];
patients_max_activation = [];

for i=1:size(filenames,1)-2
    data = load(strcat(filenames(i+2).folder,'/',filenames(i+2).name));
    fields = fieldnames(data);
    region_names = fieldnames(data.(fields{1}).labels);
    activation = data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart;
    
    %do a regional based analysis
    region_names=fieldnames(data.(fields{1}).labels);

    for j=1:21
       activation_region = activation(data.(fields{1}).labels.(region_names{j}));
       mean_region = mean(activation_region);
       median_region = median(activation_region);
       min_region = min(activation_region);
       max_region = max(activation_region);
       patients_mean_activation(i,j) = mean_region;
       patients_median_activation(i,j) = median_region;
       patients_min_activation(i,j) = min_region;
       patients_max_activation(i,j) = max_region;

       
    end

end

stds_mean = std(patients_mean_activation);
mean_means = mean(patients_mean_activation);
mean_maxs = mean(patients_max_activation);

mean_mins = mean(patients_min_activation);
figure()
hold on
for i =1:20
    scatter(linspace(1,50,21),patients_mean_activation(i,:),50);
end
ylabel('region activation times [ms]')
title('mean activation')
xticklabels( (region_names));
xticks( linspace(1,50,21));
xtickangle(45)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)


figure()
hold on
for i =1:20
    scatter(linspace(1,50,21),patients_min_activation(i,:),50,'b');
end
for i =1:20
    scatter(linspace(1,50,21),mean_mins,50,'r','filled');
end

ylabel('region activation times [ms]')
title('1st activation site of region')
xticklabels( (region_names));
xticks( linspace(1,50,21));
xtickangle(45)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)


%last activation point
figure()
hold on
for i =1:20
    scatter(linspace(1,50,21),patients_max_activation(i,:),50,'b');
end
for i =1:20
    scatter(linspace(1,50,21),mean_maxs,50,'r','filled');
end

ylabel('region activation times [ms]')
title('last activation cite of region')
xticklabels( (region_names));
xticks( linspace(1,50,21));
xtickangle(45)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)





figure()
hold on
for i =1:20
    scatter(linspace(1,21,21),patients_median_activation(i,:));
end



figure()

hist(name_dummy,22,'BarWidth',1);
xticks( linspace(1,21,21)-0.9);
xticklabels( (region_names));
xt = get(gca, 'XTick');
set(gca, 'FontSize', 9)
title('regions of late conduction')



%look at gradients of activation
gradient_activation=diff(data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart);
gradient_activation(1002) = 0;

%need to interpolate the activation on a uniform grid in order to take the
%gradient

fInterpolated = scatteredInterpolant(data.(fields{1}).nodes(:,1), data.(fields{1}).nodes(:,2),data.(fields{1}).nodes(:,3),activation); %interpolate function
for i=1:size(activation,1)
   %x_edge_lengths = 
end

% 
% for i=1:size(filenames,1)-2
%     data = load(strcat(filenames(i+2).folder,'/',filenames(i+2).name));
%     fields = fieldnames(data);
%     
% 
%     M = struct;
% 
%     M.xyz = data.(fields{1}).nodes;
%     M.tri = data.(fields{1}).mesh;
%     M.xyzATTRIBUTE =  data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart;
%     write_VTK_POLYDATA(M,strcat(fields{1},'.vtk'))
% 
% end

grad_names = dir('/data/ecgi/gradients/');

figure()

for i=1:size(filenames,1)-2
    data = load(strcat(filenames(i+2).folder,'/',filenames(i+2).name));
    fields = fieldnames(data);
    grad = csvread(strcat(grad_names(1).folder,'/',grad_names(i+2).name),1);
    grad_mag = sqrt( power(grad(:,2),2) + power(grad(:,3),2) + power(grad(:,4),2)); %find the magnitude of the powers

    subplot('Position',[(mod(i-1,5))/5 0.9-(ceil(i/5))/5 1/5 1/5])
    
    %patch('vertices', data.arvc2_baseline.nodes, 'faces', data.arvc2_baseline.mesh, 'facecolor', 'interp', 'cdata', labels(:,1), 'edgecolor', 'none') ;
    p = patch('vertices', data.(fields{1}).nodes, 'faces', data.(fields{1}).mesh,'FaceVertexCData',grad_mag,'facecolor', 'interp','edgecolor', 'none') ;
    %p.Facecolor = 'Flat';
    cb =colorbar;
    caxis([0 5])
    axis equal
    %view([50,75]) %front view
    view([225,-90]) %back view
    %%view([125,20]) % lv side view
    %view([-35,-20]) % rv side view
    axis off
    set(cb,'position',[0.95 .35 .03 .3])
    %title(strcat(filepath(end-18:end-13),' activation times [ms] lv view'))
    title(filenames(i+2).name(end-18:end-13))
    %[h,v] = IsoLine({data.(fields{1}).mesh, data.(fields{1}).nodes},data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart,10,'k');

end


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


mean_max_grad = mean(region_grad_max);

figure()
hold on
for i =1:20
    scatter(linspace(1,50,21),region_grad_max(i,:),50,'b');
end
for i =1:20
    scatter(linspace(1,50,21),mean_max_grad,50,'r','filled');
end

ylabel('activation time gradient magnitude [ms/mm]')
title('max activation gradient magnitude [ms/mm]')
xticklabels( (region_names));
xticks( linspace(1,50,21));
xtickangle(45)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)



%summary table, LV
