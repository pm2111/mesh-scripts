addpath(genpath('C:\Users\petnov\Dropbox\shared\'));
addpath(genpath('C:\Users\petnov\Dropbox\mesh scripts'))
enableVTK;

patient_sim_folder =  'D:\ARVC meshing automatic\patients\patient18\';

%H0 = read_VTK( fullfile( TORSO_MODEL_DIR , 'HEART_with_UPS_louie.vtk' ) ); %this is louie's heart!
H_mesh = read_CHASTE( 'D:\ARVC meshing automatic\patients\patient18\chaste_gmsh\HEART');
H_mesh.xyz = H_mesh.xyz*10;
%H0 = read_VTK('D:\ARVC meshing automatic\patients\patient09\ARVC09_ecgi_mesh.vtk'); %use the original ecgi mesh, not the clipped one!!
%H0 = read_VTK('D:\ARVC meshing automatic\patients\patient09\ecgi_clipped.vtk');
data = load('D:\ChrisReconstructions\ARVC18\ARVC18_Baseline.mat');
name = fieldnames(data);

H0 =struct;
H0.xyz = data.(name{1}).nodes;
H0.tri = data.(name{1}).mesh;

fields = fieldnames(data.(name{1}).labels);
    for i=1:size(fields,1)
        indices = data.(name{1}).labels.(fields{i,1});
        H0.xyzlabels([indices],1) = i;
    end

%chaste works in cm, our meshes are built in mm (BE CAREFUL WHEN LOCATING
%INDICES 



%%

H1 = read_VTK(  'D:\ARVC meshing automatic\patients\patient18\HEART.vtk'); %our target heart
[EPI1,LV1,RV1,~,MAT1] = HEARTparts( H1 );
% 
% figure()
% plotMESH( H0 ,'ne');
% hplotMESH( H1 , 'ne' ,'FaceColor',[1 1 1]*0.6)
% headlight; axis(objbounds);

%HD_loaded = load('Z:\Dropbox\patients\HD.mat');

%%

u = loadv( 'D:\ARVC meshing automatic\patients\patient18\mpp\Hmapping' , 'u' ); %transfer matrix u
%this is a function converting coordinates from H0 to H1
%it is generated by running the mpp_Heart_Mapping Function of Ernesto s
%pipeline!!


%check if the two meshes are in the same units of length by plotting side
%by side, if they are not, convert the ECGi or chaste mesh unit into mm
%H0.xyz = H0.xyz*10;

if 1
HD = H0;
HD.xyz = u( HD.xyz );
figure()
plotMESH( HD ,'ne','FaceColor',[0 0 1]*0.6);
hplotMESH( H1 , 'ne' ,'FaceColor',[1 0 0]*0.6)
headlight; axis(objbounds);
axis off
alpha(1)
end


EPI1 = Mesh(EPI1);

%%
%find the point in HD (ecgi mesh registered to the chaste mesh) corresponding to each of the nodes of thr chaste mesh 
closest_points = vtkClosestPoint(   Mesh(HD) ,EPI1.xyz);

epi_cell_property = H0.xyzlabels(closest_points);



%check that the epicardial surface contains the expected fibrotic and scar
%patterns
% 
 figure()
 plotMESH( EPI1, 'edgecolor','none' , 'cdata' , epi_cell_property' ,'facecolor','interp' )
% axis off
 headlight
% hold on
% plotMESH(LV1,'edgecolor','none' )
% plotMESH(RV1,'edgecolor','none' )
% alpha(0.3)
%remove the valve plane label
epi_cell_property(find(epi_cell_property ==1)) =0;

labelled_epi = find(epi_cell_property>0);
epi_nodes_not_labelled = find(epi_cell_property==0);

[idx,dist] = knnsearch(EPI1.xyz(labelled_epi,:),EPI1.xyz(epi_nodes_not_labelled,:)); %find the index in labelled_all which is closesst to point in uncategorised indices 

only_labelled = epi_cell_property(labelled_epi);

epi_cell_property(epi_nodes_not_labelled) = only_labelled(idx);

    fid = fopen(strcat(patient_sim_folder,'epi_cell_label','.txt'),'w');
    fprintf(fid,'%d \n',[epi_cell_property'].');
    fclose(fid)
    
figure()
plotMESH(EPI1,'edgecolor','none','cdata',epi_cell_property,'facecolor','interp')

%get rid of nodes on the valve plain
closest_points = vtkClosestPoint(   Mesh(HD) ,EPI1.xyz);


%give all of the nodes which are not part of the epicardium a 0 label

mesh_cell_property = zeros(size(H_mesh.xyz,1),1);
mesh_cell_property(1:size(epi_cell_property,1)) =epi_cell_property;


%find the nodes in the RV1 (endocardium) which are 6.5mm or further from the
%epicardium, call these septal nodes
closest_indices_in_epi=vtkClosestPoint(Mesh(EPI1),RV1.xyz);
vectors_epi_rv = EPI1.xyz(closest_indices_in_epi,:)-RV1.xyz;
distances_rv_epi = sqrt(vectors_epi_rv(:,1).*vectors_epi_rv(:,1)+vectors_epi_rv(:,2).*vectors_epi_rv(:,2)+vectors_epi_rv(:,3).*vectors_epi_rv(:,3));
septal_indices_rv = find(distances_rv_epi >7.5);

rv_labels = zeros(size(RV1.xyz,1),1);
rv_labels(septal_indices_rv)=1;

figure()
headlight
hold on
plotMESH(RV1,'edgecolor','none','cdata',rv_labels,'facecolor','interp')



septal_indices_full_mesh = vtkClosestPoint(Mesh(H_mesh),RV1.xyz(septal_indices_rv,:));
mesh_cell_property(septal_indices_full_mesh) =1;

H_mesh.xyzlabels = mesh_cell_property;
%now label the remaining nodes in the midmyocardium 
labelled_indices = find(mesh_cell_property>0);

nodes_not_labelled = find(H_mesh.xyzlabels ==0);


[idx,dist] = knnsearch(H_mesh.xyz(labelled_indices,:),H_mesh.xyz(nodes_not_labelled,:)); %find the index in labelled_all which is closesst to point in uncategorised indices 


H_mesh.xyzlabels(nodes_not_labelled) = mesh_cell_property(labelled_indices(idx));

mesh_cell_property(nodes_not_labelled)=  mesh_cell_property(labelled_indices(idx));

write_VTK(H_mesh,'D:\ARVC meshing automatic\patients\patient18\pat18_labelled_ecgi.vtk','b');


[~,idx] = sort(mesh_cell_property);

%save the coords of the different regions in separate files
H_mesh.xyz = H_mesh.xyz/10;
fields{1} = 'septum';
mkdir 'D:\ARVC meshing automatic\patients\patient18\regions';
% for i=1:size(unique(mesh_cell_property))-1
%     ind=find(mesh_cell_property ==i-1);
%     fid = fopen(strcat('D:\ARVC meshing automatic\patients\patient09\regions\',num2str(fields{i}),'.txt'),'w');
%     fprintf(fid,'%d %d %d\n',[H_mesh.xyz(ind,:)].');
%     fclose(fid)
% end
% 
%  fid = fopen(strcat('D:\ARVC meshing automatic\patients\patient18\','regions_pat05_in_mesh_vtk_ordering','.txt'),'w');
%     fprintf(fid,'%5.0f \n',[mesh_cell_property].');
%     fclose(fid)



%now save values for the conductivities scaling factors, assuming a linear
%rel between LGE intensity and conductivity of tissue (max CV reduction at
%50% LGEwhere CV is halved (conductivity is divided by 4).Assume
%conductivities change equally in the 3 fiber directions
cond =@(x) -0.009*x+1.0;
startRow = 2;

[GadLV,GadRV] = read_csv('D:\sim3d\late_gad.csv', startRow);
 GadRV = str2double(GadRV);
% GadRV(92) = 73;
% GadRV(94) = 76;
% GadRV(96) = 70;
% GadRV(99) = 73;

patient_nr=18;
regional_cond_lv =[];
sodium_cells = ones(size(H_mesh.xyz,1),1);
    for i=1:20
        j=i+20*(patient_nr-1);
        if isnan(str2double(GadLV{j})) ==0
             regional_cond_lv(j-20*(patient_nr-1)+1) = cond(str2double(GadLV{j}));
        elseif GadRV(j) >0
            regional_cond_lv(j-20*(patient_nr-1)+1) = cond(GadRV(j)*100/2-10);
        else
              regional_cond_lv(j-20*(patient_nr-1)+1)=1.000;
        end
        if isnan( regional_cond_lv(j-20*(patient_nr-1))+1) ==1
            regional_cond_lv(j-20*(patient_nr-1)+1)=1.0;
        end
        regional_cond_lv(1) =1; %septum late gad info is not present
%         fid = fopen(strcat('D:\ARVC meshing automatic\patients\patient05\conductivities\',num2str(Label{j}),'.txt'),'w');
%         fprintf(fid,'%d',regional_cond_lv(i));
%         fclose(fid)
    end


%save the coords of regions in the same file and save a number
%corresponding to the region number next to it
    fid = fopen(strcat('D:\ARVC meshing automatic\patients\patient18\','regions_full','.txt'),'w');
    fprintf(fid,'%d %d %d %d \n',[regional_cond_lv(mesh_cell_property)'].');
    fclose(fid)
    
        fid = fopen(strcat('D:\ARVC meshing automatic\patients\patient18\','mesh_cell_property','.txt'),'w');
    fprintf(fid,'%d \n',[mesh_cell_property].');
    fclose(fid)
    

H_mesh.xyzlabels = mesh_cell_property;
%write_VTK(H_mesh,'D:\ARVC meshing automatic\patients\patient05\sim_mesh_labelled.vtk','binary')
