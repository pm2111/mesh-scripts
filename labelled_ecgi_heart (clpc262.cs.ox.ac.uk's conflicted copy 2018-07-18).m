%atlas of the heart regions
% data = load('/data/ChrisReconstructions/ARVC2/ARVC2_Baseline.mat');
% plot_ventricle_regions(data)
% 
addpath('/users/petnov/Dropbox/shared/IO');
addpath('/users/petnov/Dropbox/mesh scripts/');
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

%% import our ventricular surfaces

mesh = read_VTK('/home/scratch/meshes ready/patient02/HEART_0.04.vtk');

region = csvread('/home/scratch/meshes ready/patient02/apex rv lv.csv',1);
region_names = dir('/home/scratch/meshes ready/patient02/regions/');
regions = cell(1);
%rearrange all regions into one matrix
for i=1:size(region_names,1)-2
    regions{i} =  csvread(strcat('/home/scratch/meshes ready/patient02/regions/',region_names(i+2).name),1);
end



for i=1:size(mesh.xyz,1)
    
    %check if its in labelled data
    dummy=0;
    for j=1:size(region_names,1)-2
        idx = find(regions{j}(:,1) == i);
        if idx ~= []
            break
        end
    end
    
    
end



find(norm(mesh.xyz - region(2,:))<1)

