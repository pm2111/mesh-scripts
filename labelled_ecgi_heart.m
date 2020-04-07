%atlas of the heart regions
% data = load('/data/ChrisReconstructions/ARVC2/ARVC2_Baseline.mat');
% plot_ventricle_regions(data)
% 

function labelled_ecgi_heart(filepath)
    data = load(filepath);
    fields = fieldnames(data);

    labels = zeros(size(data.(fields{1}).nodes));
    
    labels(data.(fields{1}).labels.valve_plane) = 0;
    labels(data.(fields{1}).labels.basal_anterior_lv) = 1;
    labels(data.(fields{1}).labels.mid_anterior_lv) = 2;
    labels(data.(fields{1}).labels.basal_anterolateral_lv) = 3;
    labels(data.(fields{1}).labels.mid_anterolateral_lv) = 4;
    labels(data.(fields{1}).labels.basal_inferolateral_lv) = 5;
    labels(data.(fields{1}).labels.mid_inferolateral_lv) = 6;
    labels(data.(fields{1}).labels.basal_inferior_lv) = 7;
    labels(data.(fields{1}).labels.mid_inferior_lv) = 8;
    labels(data.(fields{1}).labels.apical_anterior_lv) = 9;
    labels(data.(fields{1}).labels.apical_inferior_lv) = 10;
    labels(data.(fields{1}).labels.basal_anterior_rv) = 11;
    labels(data.(fields{1}).labels.mid_anterior_rv) = 12;
    labels(data.(fields{1}).labels.basal_anterolateral_rv) = 13;
    labels(data.(fields{1}).labels.mid_anterolateral_rv) = 14;
    labels(data.(fields{1}).labels.basal_inferolateral_rv) = 15;
    labels(data.(fields{1}).labels.mid_inferolateral_rv) = 16;
    labels(data.(fields{1}).labels.basal_inferior_rv) = 17;
    labels(data.(fields{1}).labels.mid_inferior_rv) = 18;
    labels(data.(fields{1}).labels.apical_anterior_rv) = 19;
    labels(data.(fields{1}).labels.apical_inferior_rv) = 20;


    figure()

    %patch('vertices', data.arvc2_baseline.nodes, 'faces', data.arvc2_baseline.mesh, 'facecolor', 'interp', 'cdata', labels(:,1), 'edgecolor', 'none') ;
    p = patch('vertices', data.(fields{1}).nodes, 'faces', data.(fields{1}).mesh, 'FaceVertexCData', labels(:,1), 'Facecolor','flat','edgecolor', 'none') ;
    %p.Facecolor = 'Flat';
    colorbar
    axis equal
end