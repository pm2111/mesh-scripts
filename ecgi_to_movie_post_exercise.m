%%%% ECGi to movie%%%%%%%%%%%

% need to load the arvc pat data(here pat 14 & baseline)

nr_frames = size(arvc14_post_exercise.epipots.potential);

for i=1:200
    h = figure();
    set(h, 'Visible', 'off');
    
    patch('vertices', arvc14_post_exercise.nodes, 'faces', arvc14_post_exercise.mesh, 'facecolor', 'interp', 'cdata', arvc14_post_exercise.epipots.potential(:,i), 'edgecolor', 'none') ;
    
    F(i) = getframe(h);

end


fig = figure;
movie(fig,F,1)

v = VideoWriter('arvc14PostExercise','Motion JPEG AVI');

open(v)
writeVideo(v,F);

close(v)