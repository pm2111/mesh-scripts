function plot_ecgi_rest(filepath)
    data = load(filepath);
    fields = fieldnames(data);

    
    
    %figure()

    %patch('vertices', data.arvc2_baseline.nodes, 'faces', data.arvc2_baseline.mesh, 'facecolor', 'interp', 'cdata', labels(:,1), 'edgecolor', 'none') ;
    p = patch('vertices', data.(fields{1}).nodes, 'faces', data.(fields{1}).mesh,'facecolor', 'interp', 'cdata', data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart,'edgecolor', 'none') ;
    [h,v] = IsoLine({data.(fields{1}).mesh, data.(fields{1}).nodes},data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart,30,'k');

    %p.Facecolor = 'Flat';
   % cb =colorbar;
    %cb.Label.String='activation time [ms]';
   % cb.Location= 'east';
    %cb.Limits=[0.25,0.75];
   % caxis([0 70])
    axis equal
    view([50,75]) %front view
    %view([225,-90]) %back view
    %view([125,20]) % lv side view
    %view([-35,-20]) % rv side view
    axis off
    colormap(flipud(jet));
    %set(cb,'position',[0.75 .35 .03 .3])
    %title(strcat(filepath(end-18:end-13),' activation times [ms] lv view'))
    set(gca,'Fontsize',25)
    
end