addpath('/users/petnov/Dropbox/shared/IO');
addpath('/users/petnov/Dropbox/mesh scripts/');

%% import our ventricular surfaces

mesh = read_VTK('/home/scratch/meshes ready/patient02/HEART_0.04.vtk');

region = csvread('/home/scratch/meshes ready/patient02/apex rv lv.csv',1);
region_names = dir('/home/scratch/meshes ready/patient02/regions/');
regions = cell(1);
regions_complete = cell(1);
%rearrange all regions into one matrix
for i=1:size(region_names,1)-2
    regions{i} =  csvread(strcat('/home/scratch/meshes ready/patient02/regions/',region_names(i+2).name),1);
    regions_complete{i} = csvread(strcat('/home/scratch/meshes ready/patient02/regions/',region_names(i+2).name),1);

end
tic

total_dummy=0;
for i=1:size(mesh.xyz,1)
    
    %check if its in labelled data
    dummy=0;
    for j=1:size(region_names,1)-2
        idx = find(regions{j}(:,1) == i);
        if isempty(idx)==0
            dummy =1;
            break
        end
    end
    total_dummy = total_dummy +dummy;

    %if its not in labelled data, compute distances
    if dummy==0
        for k=1:size(region_names,1)-2
            dist_dummy = [];
            for l=1:size(regions{k},1)
               
                dist_dummy(l) = norm( mesh.xyz(i,1:3) - regions{k}(l,2:4));
            end
            
            min_distances(k)=min(dist_dummy); %no need to keep the whole array of distances
            
            
        end
            [value,region_id] = min(min_distances);
            regions_complete{1,region_id}(end+1,1) = i;
            regions_complete{1,region_id}(end,2) = mesh.xyz(i,1);
            regions_complete{1,region_id}(end,3) = mesh.xyz(i,2);
            regions_complete{1,region_id}(end,4) = mesh.xyz(i,3);

            
    end
    
end
toc
