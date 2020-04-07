addpath('/users/petnov/Dropbox/shared/IO');
addpath('/users/petnov/Dropbox/mesh scripts/');

filepath = '/data/ChrisReconstructions/ARVC10/ARVC10_Baseline.mat';
data = load(filepath);
fields = fieldnames(data);

%isolines = csvread('/data/ecgi/isolines/iso01.csv',2);
[h,v] = IsoLine({data.(fields{1}).mesh, data.(fields{1}).nodes},data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart,20,'k');

%finding the breakthrough points 
activation = data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart;
[maxima, maxima_ind] = findpeaks((data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart),'Npeaks',10);
data_inv = 1.01*( data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart) - (data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart);
[minima,minima_ind] = findpeaks(data_inv,'Npeaks',10);
minima = activation(minima_ind);

[TF,P] = islocalmin(activation,'MaxNumExtrema',5,'MinProminence',30 );


%do a regional based analysis
region_names=fieldnames(data.(fields{1}).labels);
for i=1:21
   activation_region = zeros(1002,1)+10000;
   activation_region(data.(fields{1}).labels.(region_names{i}))= activation( data.(fields{1}).labels.(region_names{i}));
   [value,id] = min(activation_region);
   region_minima_ind(i)=id;
   region_minima_value(i) = value;
end

%find top 3 earliest activations
region_dummy = region_minima_value;
for i=1:4
    [value,idx] = min(region_dummy);
    earliest_activation_idx(i)= region_minima_ind(idx);
    region_dummy(idx) = 1000;
    name_dummy(i) = idx;
end


%gradients of activation:
gradient_activation=diff(activation);
gradient_activation(1002) = 0;

% 
% 
% M = struct;
% 
% M.xyz = data.(fields{1}).nodes;
% M.tri = data.(fields{1}).mesh;
% M.xyzATTRIBUTE =  data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart;
% write_VTK_POLYDATA(M,'pat1_baseline.vtk')
% 
V = zeros(2,2,size(data.arvc1_baseline.nodes,1));
V(1,1,:) = data.arvc1_baseline.nodeparam.timestamps.activationTimes-data.arvc1_baseline.nodeparam.timestamps.activationStart;
fv = isosurface(data.arvc1_baseline.nodes(:,1),data.arvc1_baseline.nodes(:,2),data.arvc1_baseline.nodes(:,3),V,15);



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
%plot3(data.(fields{1}).nodes(region_minima_ind,1),data.(fields{1}).nodes(region_minima_ind,2),data.(fields{1}).nodes(region_minima_ind,3),'r*');
plot3(data.(fields{1}).nodes(earliest_activation_idx,1),data.(fields{1}).nodes(earliest_activation_idx,2),data.(fields{1}).nodes(earliest_activation_idx,3),'r*','markers',16);
%plot3(iso_mat(:,1),iso_mat(:,2),iso_mat(:,3),'-')
%plot([LS(:,1) LD(:,1)]',[LS(:,2) LD(:,2)]','Color','k','LineWidth',2);
%[h,v] = IsoLine({data.(fields{1}).mesh, data.(fields{1}).nodes},data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart,10,'k');
axis off



grad = csvread('/data/ecgi/gradients/grad10.csv',1);
grad_mag = sqrt( power(grad(:,2),2) + power(grad(:,3),2) + power(grad(:,4),2)); %find the magnitude of the powers



figure()

    %patch('vertices', data.arvc2_baseline.nodes, 'faces', data.arvc2_baseline.mesh, 'facecolor', 'interp', 'cdata', labels(:,1), 'edgecolor', 'none') ;
    p = patch('vertices', data.(fields{1}).nodes, 'faces', data.(fields{1}).mesh,'FaceVertexCData',grad_mag,'facecolor', 'interp','edgecolor', 'none') ;
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
    %[h,v] = IsoLine({data.(fields{1}).mesh, data.(fields{1}).nodes},data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart,10,'k');


%do a regional based analysis for the max gradients
region_names=fieldnames(data.(fields{1}).labels);
for i=1:21
   %grad_region = zeros(1002,1)+10000;
   grad_region= grad_mag( data.(fields{1}).labels.(region_names{i}));
   [value,id] = max(grad_region);   
   region_minima_value(i) = value;
end


%find top 3 earliest activations
region_dummy = region_minima_value;
for i=1:4
    [value,idx] = min(region_dummy);
    earliest_activation_idx(i)= region_minima_ind(idx);
    region_dummy(idx) = 1000;
    name_dummy(i) = idx;
end


