%% body surface potentials plot
addpath('/users/petnov/Dropbox/shared/ ');
addpath('/data/Dphil/arcus/');

electrode_locations = load('/data/body surface potentials/ARVC_electrode_locations.mat');

electrodes = load('/data/meshes/ARVC_005/mpp/ECG_ELECTRODES.mat');

torso = load('/data/meshes/ARVC_005/mpp/fittedVEST.mat');

bodysp = load('/data/Data_ARVC_raw_electrodes_vest/Data_ARVC/ARVC2/ARVC2_Baseline.mat');

heart = read_VTK('/data/meshes/ARVC_005/HEART.vtk');

figure()
patch('vertices', torso.VEST.xyz, 'faces', torso.VEST.tri,'edgecolor','none')
camlight
alpha(0.4)                % set all patches transparency to 0.3
hold all
scatter3(electrode_locations.arvc2_electrodes(:,1),electrode_locations.arvc2_electrodes(:,2),electrode_locations.arvc2_electrodes(:,3),40,'filled','r')


figure()
patch('vertices', torso.BODY.xyz, 'faces', torso.BODY.tri,'edgecolor','none')
camlight
alpha(0.3)                % set all patches transparency to 0.3
hold all
scatter3(electrode_locations.arvc2_electrodes(:,1),electrode_locations.arvc2_electrodes(:,2),electrode_locations.arvc2_electrodes(:,3),40,'filled','r')

ele_names =fieldnames(electrodes.ELECTRODES);


figure()
patch('vertices', torso.BODY.xyz, 'faces', torso.BODY.tri,'edgecolor','none')
patch('vertices', heart.xyz, 'faces', heart.tri,'edgecolor','none')
%patch('vertices', torso.BODY.xyz, 'faces', torso.BODY.tri)

camlight('headlight')
alpha(0.4)                % set all patches transparency to 0.3
hold all
for i=1:12
scatter3(electrodes.ELECTRODES.(ele_names{i})(1),electrodes.ELECTRODES.(ele_names{i})(2),electrodes.ELECTRODES.(ele_names{i})(3),65,'filled','b')
end
colorbar
axis equal
axis off
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
