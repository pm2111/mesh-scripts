cd /data/ecgi/
addpath(genpath('/home/scratch/petnov/Dropbox/'));
filenames = dir('/data/ChrisReconstructions/baselines/');

figure()

for i=1:size(filenames,1)-2
    data = load(strcat(filenames(i+2).folder,'/',filenames(i+2).name));
    fields = fieldnames(data);

    subplot('Position',[(mod(i-1,5))/5 0.9-(ceil(i/5))/5 1/5 1/5])
    
        %patch('vertices', data.arvc2_baseline.nodes, 'faces', data.arvc2_baseline.mesh, 'facecolor', 'interp', 'cdata', labels(:,1), 'edgecolor', 'none') ;
        p = patch('vertices', data.(fields{1}).nodes, 'faces', data.(fields{1}).mesh,'facecolor', 'interp', 'cdata', data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart,'edgecolor', 'none') ;
        %p.Facecolor = 'Flat';
        cb =colorbar;
        caxis([0 70])
        axis equal
        view([50,75]) %front view
        %view([225,-90]) %back view
        %%view([125,20]) % lv side view
        %view([-35,-20]) % rv side view
        axis off
        set(cb,'position',[0.75 .3 .03 .3])
        set(gca,'FontSize',22)
        %title(strcat(filepath(end-18:end-13),' activation times [ms] lv view'))
        title(filenames(i+2).name(end-18:end-13))
        [h,v] = IsoLine({data.(fields{1}).mesh, data.(fields{1}).nodes},data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart,10,'k');
        cb.Label.String = 'activation times [ms] ' ;
        
end

figure()
i=1
data = load(strcat(filenames(i+2).folder,'/',filenames(i+2).name));
    fields = fieldnames(data);

    
        %patch('vertices', data.arvc2_baseline.nodes, 'faces', data.arvc2_baseline.mesh, 'facecolor', 'interp', 'cdata', labels(:,1), 'edgecolor', 'none') ;
        p = patch('vertices', data.(fields{1}).nodes, 'faces', data.(fields{1}).mesh,'facecolor', 'interp', 'cdata', data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart,'edgecolor', 'none') ;
        %p.Facecolor = 'Flat';
        cb =colorbar;
        caxis([0 70])
        axis equal
        view([50,75]) %front view
        %view([225,-90]) %back view
        %%view([125,20]) % lv side view
        %view([-35,-20]) % rv side view
        axis off
        set(cb,'position',[0.75 .3 .03 .3])
        set(gca,'FontSize',22)
        %title(strcat(filepath(end-18:end-13),' activation times [ms] lv view'))
        title(filenames(i+2).name(end-18:end-13))
        [h,v] = IsoLine({data.(fields{1}).mesh, data.(fields{1}).nodes},data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart,10,'k');
        cb.Label.String = 'activation times [ms] ' ;


gradient_activation=diff(data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart);
gradient_activation(1002) = 0;

figure()

    data = load(strcat(filenames(i+2).folder,'/',filenames(i+2).name));
    fields = fieldnames(data);

  %  subplot('Position',[(mod(i-1,5))/5 0.9-(ceil(i/5))/5 1/5 1/5])
    
    %patch('vertices', data.arvc2_baseline.nodes, 'faces', data.arvc2_baseline.mesh, 'facecolor', 'interp', 'cdata', labels(:,1), 'edgecolor', 'none') ;
    p = patch('vertices', data.(fields{1}).nodes, 'faces', data.(fields{1}).mesh,'facecolor', 'interp', 'cdata', gradient_activation,'edgecolor', 'none') ;
    %p.Facecolor = 'Flat';
    cb =colorbar;
    %caxis([-40 10])
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

