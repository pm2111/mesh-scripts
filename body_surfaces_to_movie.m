%body surface to movie

electrode_locations = load('/data/body surface potentials/ARVC_electrode_locations.mat');

torso = load('/home//scratch/meshes ready/patient02/fittedVEST.mat');

bodysp = load('/data/Data_ARVC_raw_electrodes_vest/Data_ARVC/ARVC2/ARVC2_Baseline.mat');


%one potential lead is way out of range
bodysp_filtered = bodysp.signals_raw;
for i=1:size(bodysp.signals_raw,2)
    
    maxima(i) = max(bodysp.signals_raw(:,i));
    
    if maxima(i) < 0.2 * power(10,7)
        bodysp_filtered(:,i) = bodysp.signals_raw(:,i);
        
    else
        bodysp_filtered(:,i) = [];%zeros(size(bodysp.signals_raw,1));
    end
        
end


for i=1:300
    h = figure();
    view([0 0]);

    set(h, 'Visible', 'off');
    patch('vertices', torso.BODY.xyz, 'faces', torso.BODY.tri,'edgecolor','none')
    camlight
    alpha(0.3)                % set all patches transparency to 0.3
    hold all
    scatter3(electrode_locations.arvc2_electrodes(:,1),electrode_locations.arvc2_electrodes(:,2),electrode_locations.arvc2_electrodes(:,3),150,bodysp_filtered(43000+i*10,1:256),'filled')

    
    F(i) = getframe(h);

end


fig = figure;
movie(fig,F,1)

v = VideoWriter('arvc2_baseline_body_surfaces1','Motion JPEG AVI');

open(v)
writeVideo(v,F);

close(v)