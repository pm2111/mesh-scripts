
addpath(genpath('C:\Users\petnov\Dropbox\shared - copy\'));
addpath(genpath('C:\Users\petnov\Dropbox\mesh scripts'))
enableVTK;
pats={'09'};
    patient_nr=pats{1};
    %mkdir 'D:\ARVC meshing automatic\patients\patient01\regions';
H1 = read_VTK(  strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\HEART.vtk')); %our target heart
    [EPI1,LV1,RV1,~,MAT1] = HEARTparts( H1 );
    patient_sim_folder =  strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\');

    %H0 = read_VTK( fullfile( TORSO_MODEL_DIR , 'HEART_with_UPS_louie.vtk' ) ); %this is louie's heart!
    H_mesh = read_CHASTE( strcat(patient_sim_folder,    'chaste',patient_nr,'\HEART'));
    H_mesh.xyz = H_mesh.xyz*10;
    %H0 = read_VTK('D:\ARVC meshing automatic\patients\patient09\ARVC09_ecgi_mesh.vtk'); %use the original ecgi mesh, not the clipped one!!
    %\H0 = read_VTK('D:\ARVC meshing automatic\patients\patient09\ecgi_clipped.vtk');
    if patient_nr(1)=='0'
        data = load(strcat('D:\ChrisReconstructions\ARVC',patient_nr(2),'\ARVC',patient_nr(2),'_Baseline.mat'));
    else
        data = load(strcat('D:\ChrisReconstructions\ARVC',patient_nr,'\ARVC',patient_nr,'_Baseline.mat'));
    end
        name = fieldnames(data);
    H0 =struct;
    H0.xyz = data.(name{1}).nodes;
    H0.tri = data.(name{1}).mesh;
   mapping=[0,7,8,5,6,3,4,1,2,10,9,17,18,15,16,13,14,11,12,20,19]+1;
    fields = fieldnames(data.(name{1}).labels);
        for i=1:size(fields,1)
            indices = data.(name{1}).labels.(fields{i,1});
                H0.xyzlabels([indices],1) = mapping(i);
        end
       
        
       
   
    %regions= load(strcat(patient_sim_folder,'regions_labelled.txt'));
    epi_cell_property=dlmread(strcat(patient_sim_folder,'epi_cell_property.txt'));

    %give all of the nodes which are not part of the epicardium a 0 label

     mesh_cell_property = zeros(size(H_mesh.xyz,1),1);
    mesh_cell_property(1:size(epi_cell_property,1)) =epi_cell_property;

    H_mesh.xyzlabels = mesh_cell_property;
    %now label the remaining nodes in the midmyocardium 
    labelled_indices = find(mesh_cell_property>0);

    nodes_not_labelled = find(H_mesh.xyzlabels ==0);


    [idx,dist] = knnsearch(H_mesh.xyz(labelled_indices,:),H_mesh.xyz(nodes_not_labelled,:)); %find the index in labelled_all which is closesst to point in uncategorised indices 


    H_mesh.xyzlabels(nodes_not_labelled) = mesh_cell_property(labelled_indices(idx));

    mesh_cell_property(nodes_not_labelled)=  mesh_cell_property(labelled_indices(idx));

    
    figure()
    headlight
    hold on
    patch('vertices',EPI1.xyz,'faces',EPI1.tri,'edgecolor','none','FaceVertexCData',epi_cell_property,'facecolor','interp')



    %find the nodes in the RV1 (endocardium) which are 6.5mm or further from the
    %epicardium, call these septal nodes
%     closest_indices_in_epi=vtkClosestPoint(Mesh(EPI1),RV1.xyz);
%     vectors_epi_rv = EPI1.xyz(closest_indices_in_epi,:)-RV1.xyz;
%     distances_rv_epi = sqrt(vectors_epi_rv(:,1).*vectors_epi_rv(:,1)+vectors_epi_rv(:,2).*vectors_epi_rv(:,2)+vectors_epi_rv(:,3).*vectors_epi_rv(:,3));
%     septal_indices_rv = find(distances_rv_epi >6.5);
        distance= @(vec,vec1) ((vec(:,1)-vec1(:,1)).^2+(vec(:,2)-vec1(:,2)).^2+(vec(:,3)-vec1(:,3)).^2).^0.5;
% 
%         colors=zeros(size(RV1.xyz,1),1);
%         colors(septal_indices_rv)=1;
%         figure()
%         patch('vertices',RV1.xyz,'faces',RV1.tri,'edgecolor','none','FaceVertexCData',colors,'facecolor','interp')
%         headlight
%         
%         septal_indices_rv_local=septal_indices_rv-size(EPI1.xyz,1)-size(LV1.xyz,1);
%     %call all nodes further than 6.5 mm from epi surface septal nodes
%     closest_sept_rv_point=vtkClosestPoint(Mesh(RV1),H_mesh.xyz);
%     
%     distances_sept_rv_to_every_point=distance(H_mesh.xyz,RV1.xyz(closest_sept_rv_point,:));
%     
%     septal_nodes_rv=find(distances_sept_rv_to_every_point<0.8);
%     %all nodes further than 6.5 mm from epi surf are septal
      closest_idx_in_epi=vtkClosestPoint(Mesh(EPI1),H_mesh.xyz);
    distances=distance(EPI1.xyz(closest_idx_in_epi,:),H_mesh.xyz);
     septal_indices = find(distances >12.5);

     mesh_cell_property(septal_indices)=1;
     figure()
    headlight
    hold on
    patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',mesh_cell_property,'facecolor','interp')
    colorbar
    camlight 
    axis off
    
    
    
    %septum stuff below, not necessary for every patient
    septal_nodes=find(mesh_cell_property==1);
    nearestlv=vtkClosestPoint(Mesh(LV1),H_mesh.xyz(septal_nodes,:));
    nearestrv=vtkClosestPoint(Mesh(RV1),H_mesh.xyz(septal_nodes,:));
        
        nearestEPI=vtkClosestPoint(Mesh(EPI1),H_mesh.xyz(septal_nodes,:));

distlv=distance(LV1.xyz(nearestlv,:),H_mesh.xyz(septal_nodes,:));
distrv=distance(RV1.xyz(nearestrv,:),H_mesh.xyz(septal_nodes,:));
distepi=distance(EPI1.xyz(nearestEPI,:),H_mesh.xyz(septal_nodes,:));
for i=1:size(distlv,1)
    if distrv(i) >= distepi(i) && distrv(i)*2>= distlv(i)
        region_old(i)=1;
    
    else % distlv(i) >= distepi(i) && distlv(i) >=distrv(i)
                region_old(i)=2;

    end
   if distepi(i)*2 >= distlv(i) && distepi(i) >=distrv(i)

        if distlv(i) < 2/3*(distlv(i)+distrv(i))
            region_old(i)=3;
        else
            region_old(i)=4;
        end
    end

end

    

    septal_lv=septal_nodes(region_old==3);
    septal_rv=septal_nodes(region_old==4);
    
    mesh_cell_property(septal_lv)=22;
        mesh_cell_property(septal_rv)=23;

        dlmwrite(strcat(patient_sim_folder,'regions_septal','.csv'),mesh_cell_property)

%     mesh_cell_property(septal_lv)=22;
%         mesh_cell_property(septal_rv)=23;
% 
%     H_mesh.xyzlabels(septal_lv)=22;
%         H_mesh.xyzlabels(septal_rv)=23;
    
    
    
        H_mesh.xyzregions=regions;
        H_mesh.triORTHO=0;
    write_VTK(H_mesh,strcat(patient_sim_folder,'results\regions_labelled_sept.vtk'))


    figure()
    headlight
    hold on
    patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',regions,'facecolor','interp')
    colorbar

    
         mesh_cell_property=dlmread(strcat(patient_sim_folder,'regions_labelled_septal.csv'));
       cond =@(x) -0.011*x+1.0;
        cond5= @(x) -0.0125*x+1.0;
    cond6= @(x) -0.014*x+1.0;
    cond7=@(x) -0.01485*x+1.0;
    startRow = 2;

    [GadLV,GadRV] = read_csv('D:\sim3d\late_gad.csv', startRow);
     GadRV = str2double(GadRV);
    % GadRV(92) = 73;
    % GadRV(94) = 76;
    % GadRV(96) = 70;
    % GadRV(99) = 73;
    
    patient_nr=pats{1};
    patient_nr_str=patient_nr;
    if patient_nr(1) =='0' ==0
        patient_nr=str2double(patient_nr);
    else
        patient_nr= str2double(patient_nr(2));

    end
    regional_cond_lv =[];
    CV_scaling=ones(23,1);
              CV_scaling_int=ones(23,1);
  CV_scaling_strong=ones(23,1);
  CV_scaling_vstrong=ones(23,1);

    sodium_cells = ones(size(H_mesh.xyz,1),1);
    regional_cond_lv =[];

        for i=1:20
            j=i+20*(patient_nr-1);
            if isnan(str2double(GadLV{j})) ==0
                 if(GadLV{j}==0)
                    CV_scaling(j-20*(patient_nr-1)+1) = 1;
                                        CV_scaling_int(j-20*(patient_nr-1)+1) = 1;
                                        CV_scaling_strong(j-20*(patient_nr-1)+1) = 1;
                                        CV_scaling_vstrong(j-20*(patient_nr-1)+1) = 1;
                                        regional_cond_lv(j-20*(patient_nr-1)+1) = 1;

                 else

                      CV_scaling(j-20*(patient_nr-1)+1) = sqrt(cond(str2double(GadLV{j})));
                      CV_scaling_int(j-20*(patient_nr-1)+1) = sqrt(cond5(str2double(GadLV{j})));
                      CV_scaling_strong(j-20*(patient_nr-1)+1) = sqrt(cond6(str2double(GadLV{j})));
                      CV_scaling_vstrong(j-20*(patient_nr-1)+1) = sqrt(cond7(str2double(GadLV{j})));
                         regional_cond_lv(j-20*(patient_nr-1)+1) = cond5(str2double(GadLV(j)));

                 end

            elseif GadRV(j) >0
                regional_cond_lv(j-20*(patient_nr-1)+1) = cond5(GadRV(j)*100/3);

                CV_scaling(j-20*(patient_nr-1)+1) =sqrt(cond(GadRV(j)*100/3));
                CV_scaling_int(j-20*(patient_nr-1)+1) =sqrt(cond5(GadRV(j)*100/3));
                CV_scaling_strong(j-20*(patient_nr-1)+1) =sqrt(cond6(GadRV(j)*100/3));
                CV_scaling_vstrong(j-20*(patient_nr-1)+1) =sqrt(cond7(GadRV(j)*100/3));

            else
                                  regional_cond_lv(j-20*(patient_nr-1)+1)=1.000;

                  CV_scaling(j-20*(patient_nr-1)+1)=1.000;
                  CV_scaling_int(j-20*(patient_nr-1)+1)=1.000;
                  CV_scaling_strong(j-20*(patient_nr-1)+1)=1.000;
                  CV_scaling_vstrong(j-20*(patient_nr-1)+1)=1.000;

            end
           
            CV_scaling(1) =1; %septum late gad info is not present
            CV_scaling_int(1) =1; %septum late gad info is not present
            CV_scaling_strong(1) =1; %septum late gad info is not present
            CV_scaling_vstrong(1) =1; %septum late gad info is not present

    %         fid = fopen(strcat('D:\ARVC meshing automatic\patients\patient05\conductivities\',num2str(Label{j}),'.txt'),'w');
    %         fprintf(fid,'%d',regional_cond_lv(i));
    %         fclose(fid)
        end
% %   

        CV_scaling=CV_scaling([0,7,8,5,6,3,4,1,2,10,9,17,18,15,16,13,14,11,12,20,19]+1);
        CV_scaling_int=CV_scaling_int([0,7,8,5,6,3,4,1,2,10,9,17,18,15,16,13,14,11,12,20,19]+1);
        CV_scaling_strong=CV_scaling_strong([0,7,8,5,6,3,4,1,2,10,9,17,18,15,16,13,14,11,12,20,19]+1);
        CV_scaling_vstrong=CV_scaling_vstrong([0,7,8,5,6,3,4,1,2,10,9,17,18,15,16,13,14,11,12,20,19]+1);


%   regional_cond_lv(22)=0.5;
%     regional_cond_lv(23)=0.5;
% 
%         CV_scaling(find(strcmp(fields,'septum')))=1.0;
%         CV_scaling(22)=0.55;
%         CV_scaling(23)=0.55;

   % regions=csvread(strcat(patient_sim_folder,'regions_labelled.txt'));
    %fib=dlmread(strcat(patient_sim_folder,'septal_fib.csv'),',',[1 0 0 0]);
    ts_scaling=1.0./CV_scaling;
    ts_scaling_int=1.0./CV_scaling_int;
    ts_scaling_strong=1.0./CV_scaling_strong;
    ts_scaling_vstrong=1.0./real(CV_scaling_vstrong);
    
    
   % ts_scaling_septal=ts_scaling_strong;
    %ts_scaling_septal(23)=4;
    
    dlmwrite(strcat(patient_sim_folder,'eikonal09_fine_nodeScaling_1f','.csv'),[ts_scaling(mesh_cell_property(:,1))'].');
    dlmwrite(strcat(patient_sim_folder,'eikonal09_fine_nodeScaling_2f','.csv'),[ts_scaling_int(mesh_cell_property(:,1))'].');
    dlmwrite(strcat(patient_sim_folder,'eikonal09_fine_nodeScaling_3f','.csv'),[ts_scaling_strong(mesh_cell_property(:,1))'].');
    dlmwrite(strcat(patient_sim_folder,'eikonal09_fine_nodeScaling_4f','.csv'),[ts_scaling_vstrong(mesh_cell_property(:,1))'].');
    dlmwrite(strcat(patient_sim_folder,'eikonal09_fine_nodeScaling_3f_septal_and_rv','.csv'),[ts_scaling_septal(mesh_cell_property(:,1))'].');

    
   dummy= dlmread(strcat(patient_sim_folder,'eikonal09_fine_nodeScaling_3f_septal_5f','.csv'),',');
       dummy1= dlmread(strcat(patient_sim_folder,'eikonal09_fine_nodeScaling_3f','.csv'),',');

    roots=csvread(strcat(patient_sim_folder,'eikonal09_fine_rootNodes','.csv'));
         
   roots=[roots(1:3),roots([5:7])];
   
   dlmwrite(strcat(patient_sim_folder,'eikonal06_fine_rootNodes_fib','.csv'),roots)
    figure()
    headlight
    hold on
    patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',[ts_scaling_strong(mesh_cell_property)],'facecolor','interp')
    colorbar
    camlight 
    caxis( [0 5])
    axis off
    scatter3(H_mesh.xyz(roots,1),H_mesh.xyz(roots,2),H_mesh.xyz(roots,3),105,'filled','r')
    axis off
    
    
        fid = fopen(strcat(patient_sim_folder,'regions_full_light','.txt'),'w');
        fprintf(fid,'%d %d %d \n',[CV_scaling(cell_regions(:,1))' CV_scaling(cell_regions(:,1)) regions].');
        fclose(fid)