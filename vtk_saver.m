addpath('/users/petnov/Dropbox/shared/IO');
addpath('/users/petnov/Dropbox/mesh scripts/');

filepath = '/data/ChrisReconstructions/ARVC1/ARVC1_Baseline.mat';
data = load(filepath);
fields = fieldnames(data);

%isolines = csvread('/data/ecgi/isolines/iso01.csv',2);
[LS,LD,c] = isolines(data.arvc1_baseline.nodes,data.arvc1_baseline.mesh,data.arvc1_baseline.nodeparam.timestamps.activationTimes-data.arvc1_baseline.nodeparam.timestamps.activationStart,30.0);
[h,v] = IsoLine({data.(fields{1}).mesh, data.(fields{1}).nodes},data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart,20,'k');

%finding the breakthrough points 
activation = data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart;
[maxima, maxima_ind] = findpeaks((data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart),'Npeaks',10);
data_inv = 1.01*( data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart) - (data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart);
[minima,minima_ind] = findpeaks(data_inv,'Npeaks',10);
minima = activation(minima_ind);

[TF,P] = islocalmin(activation,'MaxNumExtrema',5,'MinProminence',30 );


%gradients of activation:
gradient_activation=diff(activation);
gradient_activation(1002) = 0;
% 
% M = struct;
% 
% M.xyz = data.arvc1_baseline.nodes;
% M.tri = data.arvc1_baseline.mesh;
% M.xyzATTRIBUTE =  data.arvc1_baseline.nodeparam.timestamps.activationTimes-data.arvc1_baseline.nodeparam.timestamps.activationStart;
% write_VTK_POLYDATA(M,'pat1_baseline.vtk')
% 
% V = zeros(2,2,size(data.arvc1_baseline.nodes,1));
% V(1,1,:) = data.arvc1_baseline.nodeparam.timestamps.activationTimes-data.arvc1_baseline.nodeparam.timestamps.activationStart;
% fv = isosurface(data.arvc1_baseline.nodes(:,1),data.arvc1_baseline.nodes(:,2),data.arvc1_baseline.nodes(:,3),V,15);
% 


figure()

%patch('vertices', data.arvc2_baseline.nodes, 'faces', data.arvc2_baseline.mesh, 'facecolor', 'interp', 'cdata', labels(:,1), 'edgecolor', 'none') ;
p = patch('vertices', data.(fields{1}).nodes, 'faces', data.(fields{1}).mesh,'facecolor', 'interp', 'cdata', data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart,'edgecolor', 'none') ;
%p.Facecolor = 'Flat';
cb = colorbar;
set(cb,'position',[0.75 .35 .03 .3])
caxis([0 70])
axis equal
hold all
view([50,-75])
%plot3(data.(fields{1}).nodes(TF,1),data.(fields{1}).nodes(TF,2),data.(fields{1}).nodes(TF,3),'r*');
%plot3(data.(fields{1}).nodes(a,1),data.(fields{1}).nodes(a,2),data.(fields{1}).nodes(a,3),'r*');
%plot3(iso_mat(:,1),iso_mat(:,2),iso_mat(:,3),'-')
%plot([LS(:,1) LD(:,1)]',[LS(:,2) LD(:,2)]','Color','k','LineWidth',2);
%[h,v] = IsoLine({data.(fields{1}).mesh, data.(fields{1}).nodes},data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart,10,'k');
colorbar off
axis off
saveas(gcf,'pat1.png')

image = imread('/data/ecgi/pat1.png');

min = imregionalmin(image);


%gradients of activation
figure()

%patch('vertices', data.arvc2_baseline.nodes, 'faces', data.arvc2_baseline.mesh, 'facecolor', 'interp', 'cdata', labels(:,1), 'edgecolor', 'none') ;
p = patch('vertices', data.(fields{1}).nodes, 'faces', data.(fields{1}).mesh,'facecolor', 'interp', 'cdata', gradient_activation,'edgecolor', 'none') ;
%p.Facecolor = 'Flat';
cb = colorbar;
set(cb,'position',[0.75 .35 .03 .3])
caxis([-40 10])
axis equal
hold all
view([50,-75])
%plot3(iso_mat(:,1),iso_mat(:,2),iso_mat(:,3),'-')
%plot([LS(:,1) LD(:,1)]',[LS(:,2) LD(:,2)]','Color','k','LineWidth',2);
%[h,v] = IsoLine({data.(fields{1}).mesh, data.(fields{1}).nodes},data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart,10,'k');

