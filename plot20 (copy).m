cd /data/ecgi/
addpath('/users/petnov/Dropbox/mesh scripts/');
filenames = dir('/data/ChrisReconstructions/baselines/');

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
    for j=1:4
        [value,idx] = min(region_dummy);
        earliest_activation_idx(j)= region_minima_ind(idx);
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


gradient_activation=diff(data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart);
gradient_activation(1002) = 0;
% 
% figure()
% 
% for i=1:size(filenames,1)-2
%     data = load(strcat(filenames(i+2).folder,'/',filenames(i+2).name));
%     fields = fieldnames(data);
% 
%     subplot('Position',[(mod(i-1,5))/5 0.9-(ceil(i/5))/5 1/5 1/5])
%     
%     %patch('vertices', data.arvc2_baseline.nodes, 'faces', data.arvc2_baseline.mesh, 'facecolor', 'interp', 'cdata', labels(:,1), 'edgecolor', 'none') ;
%     p = patch('vertices', data.(fields{1}).nodes, 'faces', data.(fields{1}).mesh,'facecolor', 'interp', 'cdata', gradient_activation,'edgecolor', 'none') ;
%     %p.Facecolor = 'Flat';
%     cb =colorbar;
%     %caxis([-40 10])
%     axis equal
%     %view([50,75]) %front view
%     view([225,-90]) %back view
%     %%view([125,20]) % lv side view
%     %view([-35,-20]) % rv side view
%     axis off
%     set(cb,'position',[0.95 .35 .03 .3])
%     %title(strcat(filepath(end-18:end-13),' activation times [ms] lv view'))
%     title(filenames(i+2).name(end-18:end-13))
%     %[h,v] = IsoLine({data.(fields{1}).mesh, data.(fields{1}).nodes},data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart,10,'k');
% 
% end