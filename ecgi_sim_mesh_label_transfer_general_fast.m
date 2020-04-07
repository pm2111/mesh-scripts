
addpath(genpath('C:\Users\petnov\Dropbox\shared - copy\'));
addpath(genpath('C:\Users\petnov\Dropbox\mesh scripts'))
addpath('C:\Users\petnov\Dropbox\shared - copy\MESHES\');

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
        
        fields_rearr=fields(mapping);
        

    regions= load(strcat(patient_sim_folder,'regions_labelled.txt'));
    
%     figure()
%       patch('vertices',data.arvc6_baseline.nodes,'faces',data.arvc6_baseline.mesh,'edgecolor','none','FaceVertexCData',H0.xyzlabels,'facecolor','interp')

      %smooth the regions
      regions_ecgi=H0.xyzlabels;
      %if a number of neighboring triangles are of another region, convert to that region
      
      %how do we find the neighbouring triangles?
      %can search in tris to see where indices appear
    
     
%        figure()
%       patch('vertices',data.arvc6_baseline.nodes,'faces',data.arvc6_baseline.mesh,'edgecolor','none','FaceVertexCData',regions_ecgi,'facecolor','interp')
%       hold on
%       patch('vertices',data.arvc6_baseline.nodes,'faces',data.arvc6_baseline.mesh(uncertain_tris,:),'edgecolor','none','facecolor','r')
%         
      %colorbar
      
      %assign uncertain triangles based on geometric proximity
      %find first 5 neareast neighbors
      %check if the nearest neighbors are in uncertain list
%       to_be_removed=[];
%       for i=1:size(nearest_neighbors_uncertain,1)
%           if isempty(find(H0.tri(uncertain_tris,:)==nearest_neighbors_uncertain(i)))==0
%               to_be_removed=vertcat(to_be_removed,find(H0.tri(uncertain_tris,:)==nearest_neighbors_uncertain(i)));
%           end
%       end
        x=1;
            
%        
%              figure()
%       patch('vertices',data.arvc6_baseline.nodes,'faces',data.arvc6_baseline.mesh,'edgecolor','none','FaceVertexCData',regions_ecgi,'facecolor','interp')
%       hold on
      
      u = loadv( 'D:\ARVC meshing automatic\patients\patient09\mpp\Hmapping_ecgi' , 'u' ); %transfer matrix u

      HD = H0;
HD.xyz = u( HD.xyz );
% figure()
% plotMESH( HD ,'ne','FaceColor',[0 0 1]*0.6);
% hplotMESH( H1 , 'ne' ,'FaceColor',[1 0 0]*0.6)
% headlight; axis(objbounds);


EPI1 = Mesh(EPI1);

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
%       closest_idx_in_epi=vtkClosestPoint(Mesh(EPI1),H_mesh.xyz*10);
%     distances=distance(EPI1.xyz(closest_idx_in_epi,:),H_mesh.xyz*10);
    septal_nodes=find(regions==1);
    
    
    
    %septum stuff below, not necessary for every patient
%     
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
    
    regions(septal_lv)=22;
        regions(septal_rv)=23;

        
    regions=csvread(strcat(patient_sim_folder,'regions_labelled_septal.csv'));
          nearest_neighbors=knnsearch(H_mesh.xyz,H_mesh.xyz,'K',200);      

        regions_experimental=regions;
                       
        all_neighbors_all_nodes=nearest_neighbors(1:size(H_mesh.xyz,1),2:200);

            for n=1:100
                 all_regions_neightbors_nodes=regions_experimental(all_neighbors_all_nodes);
            
                    hists_regions_neightbors_nodes=histc(all_regions_neightbors_nodes,[1:23],2); %every 4 rows concern 1 tri
                    [~,id]=max(hists_regions_neightbors_nodes,[],2);
                    regions_experimental=id;
                    disp(['changed' num2str(sum(regions~=regions_experimental))]);
            end
            
            
            
            dlmwrite(strcat(patient_sim_folder,'regions_labelled_septal_smooth.csv'),regions_experimental,'precision',10);

%     
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
    
    
    
    figure()
    headlight
    hold on
    patch('vertices',H1.xyz,'faces',H1.tri,'edgecolor','none','FaceVertexCData',regions,'facecolor','interp')
    colorbar
    
    
        regions=csvread(strcat(patient_sim_folder,'regions_labelled_septal_smooth.csv'));


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
%   regional_cond_lv(22)=0.5;
%     regional_cond_lv(23)=0.5;
% 
%         CV_scaling(find(strcmp(fields,'septum')))=1.0;
%         CV_scaling(22)=0.55;
%         CV_scaling(23)=0.55;
   % dlmwrite(strcat(patient_sim_folder,'regions_labelled_septal.csv'),regions);

    regions=csvread(strcat(patient_sim_folder,'regions_labelled_septal.csv'));
          nearest_neighbors=knnsearch(H_mesh.xyz,H_mesh.xyz,'K',20);      

        regions_experimental=regions;
        
            all_neighbors_all_nodes_in_tri=nearest_neighbors(H_mesh.tri(:,:),2:20);
            all_regions_neightbors_tri=regions_experimental(all_neighbors_all_nodes_in_tri);
            all_neighbors_all_nodes_in_tri=[];
            
             hists=histc(all_regions_neightbors_tri(:,:),[1:23],2);

         
       
%     
%        for n=1:100
%            regions_tet=regions(H_mesh.tri);
%            nr_regions=unique(regions_tet,'rows');
%            if size(unique(regions_tet),1)==2
%                for k=1:4
%                    idx= find(regions_tet==regions_tet(k));
%                     if size(idx,1)>2
%                         regions(H_mesh.tri(i,:))=regions_tet(k);
%                         break;
%                     end
%                end
%            end
%       end
%       
         figure()
    headlight
    hold on
    patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',regions_experimental,'facecolor','interp')
    colorbar
    camlight 
    axis off
    
    
    
    ts_scaling=1.0./CV_scaling;
    ts_scaling_int=1.0./CV_scaling_int;
    ts_scaling_strong=1.0./CV_scaling_strong;
    ts_scaling_vstrong=1.0./real(CV_scaling_vstrong);
    
 %   ts_scaling_int(21)=4;
%     
%     ts_scaling_strong=ones(1,23);
%     ts_scaling_strong(20)=3;
%         ts_scaling_strong(21)=3;
%     ts_scaling_strong(23)=3;

    
   % for i=1:3
   
%      ts_scaling(22)=1+i;
%     ts_scaling_int(22)=1+i;
%     ts_scaling_strong(22)=1+i;
    
    
    dlmwrite(strcat(patient_sim_folder,'eikonal06_fine_nodeScaling_',num2str(1),'f','.csv'),[ts_scaling(regions(:,1))'].');
    dlmwrite(strcat(patient_sim_folder,'eikonal06_fine_nodeScaling_',num2str(2),'f','.csv'),[ts_scaling_int(regions(:,1))'].');
    dlmwrite(strcat(patient_sim_folder,'eikonal06_fine_nodeScaling_',num2str(3),'f','.csv'),[ts_scaling_strong(regions(:,1))'].');
   % end
      
        fid = fopen(strcat(patient_sim_folder,'regions_full_light','.txt'),'w');
        fprintf(fid,'%d %d %d \n',[regional_cond_lv(regions(:,1))' regional_cond_lv(regions(:,1))' regions].');
        fclose(fid)
    
    
    roots=csvread(strcat(patient_sim_folder,'eikonal06_fine_rootNodes','.csv'));
         

   
   root_pos=[3.26,-9.915,-0.9604];
   LV1_new_root = vtkClosestPoint( Mesh(H_mesh) , root_pos );

   
      rootsNew=[roots(1:4),LV1_new_root];

   dlmwrite(strcat(patient_sim_folder,'eikonal06_fine_rootNodes_rv_ant_middle','.csv'),rootsNew)
   
   acts=dlmread(strcat(patient_sim_folder,'results\Eikonal\exp22\10.0x_3f_0.5fiber_roots_9_no_rv_pur_0\10.0x_3f_0.5fiber_roots_9_no_rv_pur_0_NEW_pred_ATMap.csv'));
   
    figure()
    headlight
    hold on
    patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',[ts_scaling_int(regions(:,1))']','facecolor','interp')
    colorbar
    camlight 
    axis off
    scatter3(H_mesh.xyz(roots,1),H_mesh.xyz(roots,2),H_mesh.xyz(roots,3),105,'filled','r')
    axis off
    view([-109 18])
    
     figure()
    headlight
    hold on
    patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',acts,'facecolor','interp')
    colorbar
    camlight 
    axis off
    scatter3(H_mesh.xyz(rootsNew,1),H_mesh.xyz(rootsNew,2),H_mesh.xyz(rootsNew,3),105,'filled','r')
    axis off
    view([-109 18])