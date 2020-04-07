%%%% ECGi to movie%%%%%%%%%%%

% need to load the arvc pat data(here pat 14 & baseline)
close all
clear all
data = load('/data/ChrisReconstructions/ARVC18/ARVC18_Baseline.mat');
nr_frames = size(data.arvc18_baseline.epipots.potential,2);
for i=1:nr_frames
    h = figure();
    set(h, 'Visible', 'off');
    h = patch('vertices', data.arvc18_baseline.nodes, 'faces', data.arvc18_baseline.mesh, 'facecolor', 'interp', 'cdata', data.arvc18_baseline.epipots.potential(:,i), 'edgecolor', 'none') ;
    
    F(i) = getframe(h);

end


fig = figure;
movie(fig,F,1)

v = VideoWriter('/data/ecgi videos/arvc18Baseline_full','Motion JPEG AVI');

open(v)
writeVideo(v,F);

close(v)

