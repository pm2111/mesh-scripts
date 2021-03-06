
addpath(genpath('C:\Users\petnov\Dropbox\shared - copy\'));
addpath(genpath('C:\Users\petnov\Dropbox\mesh scripts'))
addpath('C:\Users\petnov\Dropbox\shared - copy\MESHES\');

enableVTK;
pats={'06'};
    patient_nr=pats{1};
    %mkdir 'D:\ARVC meshing automatic\patients\patient01\regions';
H1 = read_VTK(  strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\HEART_0.70.vtk')); %our target heart
    [EPI1,LV1,RV1,~,MAT1] = HEARTparts( H1 );
    patient_sim_folder =  strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\');

    %H0 = read_VTK( fullfile( TORSO_MODEL_DIR , 'HEART_with_UPS_louie.vtk' ) ); %this is louie's heart!
    H_mesh = read_CHASTE( strcat(patient_sim_folder,    'chaste06_coarse','\HEART'));
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
        
%         fields_rearr=fields(mapping);
%         
% 
%     regions= load(strcat(patient_sim_folder,'regions_labelled.txt'));
%     
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
%       colorbar
      
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
            
       
             figure()
      patch('vertices',data.arvc6_baseline.nodes,'faces',data.arvc6_baseline.mesh,'edgecolor','none','FaceVertexCData',regions_ecgi,'facecolor','interp')
      hold on
      
      u = loadv( 'D:\ARVC meshing automatic\patients\patient06\mpp\Hmapping_ecgi_coarse' , 'u' ); %transfer matrix u

      
      figure()
      HD = H0;
HD.xyz = u( HD.xyz );
figure()
plotMESH( HD ,'ne','FaceColor',[0 0 1]*0.6);
hplotMESH( H1 , 'ne' ,'FaceColor',[1 0 0]*0.6)
headlight; axis(objbounds);

    EPI1 = Mesh(EPI1);

    %%
    %find the point in HD (ecgi mesh registered to the chaste mesh) corresponding to each of the nodes of thr chaste mesh 
    closest_points = vtkClosestPoint(   Mesh(H0) ,EPI1.xyz);

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
    colorbar
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
    septal_indices_rv = find(distances_rv_epi >6.5);
        distance= @(vec,vec1) ((vec(:,1)-vec1(:,1)).^2+(vec(:,2)-vec1(:,2)).^2+(vec(:,3)-vec1(:,3)).^2).^0.5;

    %call all nodes further than 6.5 mm from epi surface septal nodes
    closest_idx_in_epi=vtkClosestPoint(Mesh(EPI1),H_mesh.xyz*10);
    distances=distance(EPI1.xyz(closest_idx_in_epi,:),H_mesh.xyz*10);
    septal_nodes=find(distances>11.0);
    
    
    rv_labels = zeros(size(RV1.xyz,1),1);
    rv_labels(septal_indices_rv)=1;

    
    
    figure()
    headlight
    hold on
    patch('vertices',RV1.xyz,'faces',RV1.tri,'edgecolor','none','FaceVertexCData',rv_labels,'facecolor','interp')


    
    septal_indices_full_mesh = vtkClosestPoint(Mesh(H_mesh),RV1.xyz(septal_indices_rv,:));
    mesh_cell_property(septal_indices_full_mesh) =1;

    H_mesh.xyzlabels = mesh_cell_property;
    %now label the remaining nodes in the midmyocardium 
    labelled_indices = find(mesh_cell_property>0);

    nodes_not_labelled = find(H_mesh.xyzlabels ==0);


    [idx,dist] = knnsearch(H_mesh.xyz(labelled_indices,:),H_mesh.xyz(nodes_not_labelled,:)); %find the index in labelled_all which is closesst to point in uncategorised indices 


    H_mesh.xyzlabels(nodes_not_labelled) = mesh_cell_property(labelled_indices(idx));

    mesh_cell_property(nodes_not_labelled)=  mesh_cell_property(labelled_indices(idx));

    %write_VTK(H_mesh,'D:\ARVC meshing automatic\patients\patient05\pat05_labelled_ecgi.vtk','b');

    [~,idx] = sort(mesh_cell_property);
    
    
        H_mesh.xyzregions=mesh_cell_property;
        H_mesh.triORTHO=0;
    write_VTK(H_mesh,strcat(patient_sim_folder,'results\regions_labelled_sept_coarse.vtk'))
  
   % H_recov=read_VTK(strcat(patient_sim_folder,'results\regions_labelled_sept_coarse.vtk'));
  
    regions= mesh_cell_property;

    %reassign regions of the LV which enter the RV
    %region 2 => region 12 if the distance to the RV endo is smaller to the
    %lv endo
    pairs=[2 3 8 9 10 11 ;12 13 18 19 20 21 ]';
    
    mesh_cell_property2=mesh_cell_property;
    distance=@(x,y) sqrt((x(:,1)-y(:,1)).^2+(x(:,2)-y(:,2)).^2+(x(:,3)-y(:,3)).^2);
    k=1;
    for i = 1:size(pairs,1)
    j=pairs(i,1);
    indices_currently_region2=find(mesh_cell_property==j);

    closest_lv=vtkClosestPoint(LV1,H_mesh.xyz(indices_currently_region2,:));
        closest_rv=vtkClosestPoint(RV1,H_mesh.xyz(indices_currently_region2,:));

    
    dists_lv=distance(LV1.xyz(closest_lv,:),H_mesh.xyz(indices_currently_region2,:));
        dists_rv=distance(RV1.xyz(closest_rv,:),H_mesh.xyz(indices_currently_region2,:));

        swap = find(dists_rv<dists_lv);
        
        mesh_cell_property2(indices_currently_region2(swap))=pairs(k,2);
        k=k+1;
    end
        
        figure()
    headlight
    hold on
    patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',mesh_cell_property2,'facecolor','interp')
    colorbar
    camlight 
    axis off
    
        
        
     %   regions=csvread(strcat(patient_sim_folder,'regions_labelled_septal_smooth.csv'));


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

%    regions=csvread(strcat(patient_sim_folder,'regions_labelled_septal_smooth.csv'));
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
%          figure()
%     headlight
%     hold on
%     patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',regions_experimental,'facecolor','interp')
%     colorbar
%     camlight 
%     axis off
    
    
    %fibrotic scaling factors 
  
    RV_scal=[1,0.8,0.6,0.4];
    LV_scal=[1,0.8,0.6,0.4];

   % ts_scaling_vstrong=1.0./real(CV_scaling_vstrong);
%     
%        for j=1:4
            for i=1:5
           
             ts_scaling=1.0./CV_scaling;
             ts_scaling_int=1.0./CV_scaling_int;
             ts_scaling_strong=1.0./CV_scaling_strong;
             ts_scaling_no_fib=ones(23,1);

             %make large LGE patch homogenious outside
             ts_scaling(find(ts_scaling==max(ts_scaling)))=ts_scaling(end-2);
             ts_scaling_int(find(ts_scaling_int==max(ts_scaling_int)))=ts_scaling_int(end-2);
             ts_scaling_strong(find(ts_scaling_strong==max(ts_scaling_strong)))=ts_scaling_strong(end-2);

             
             scar_dists=distance(H_mesh.xyz,[51.26,-93.27,-10.02]);
             scar_ids=find(scar_dists<5); %semicircle of radius 1cm
             border_ids=find(scar_dists>5 & scar_dists <7.5);

%            ts_scaling_no_fib(12:21)= ts_scaling_no_fib(12:21)*1/RV_scal(j);
%            ts_scaling_no_fib(1:11)= ts_scaling_no_fib(1:11)*1/LV_scal(i);
% 
%            ts_scaling([12:21])= ts_scaling([12:21])*1/RV_scal(j);
%                       ts_scaling([1:11])= ts_scaling([1:11])*1/LV_scal(i);
% 
%            ts_scaling_int([12:21])= ts_scaling_int([12:21])*1/RV_scal(j);
%                       ts_scaling_int([1:11])= ts_scaling_int([1:11])*1/LV_scal(i);
% 
%            ts_scaling_strong([12:21])= ts_scaling_strong([12:21])*1/RV_scal(j);
%            ts_scaling_strong([1:11])= ts_scaling_strong([1:11])*1/LV_scal(i);
% 
%            
           no_lge= ts_scaling_no_fib(mesh_cell_property2);
            light_lge=ts_scaling(mesh_cell_property2);
           int_lge=ts_scaling_int(mesh_cell_property2);
          adv_lge=ts_scaling_strong(mesh_cell_property2);
          
          
            no_lge(scar_ids)=2*i;
            
            light_lge(scar_ids)=4*max(light_lge);
            light_lge(border_ids)=20*i;

             int_lge(scar_ids)=4*max(int_lge);
            adv_lge(scar_ids)=4*max(adv_lge);

%            no_lge(H_mesh.rv)=1;
%            light_lge(H_mesh.rv)=1;
%            int_lge(H_mesh.rv)=1;
%            adv_lge(H_mesh.rv)=1;
%            
%            no_lge(H_mesh.lv)=1;
%            light_lge(H_mesh.lv)=1;
%            int_lge(H_mesh.lv)=1;
%            adv_lge(H_mesh.lv)=1;
        
           
           % dlmwrite(strcat(patient_sim_folder,'scalings\eikonal06_coarse_fib_',num2str(0),'.csv'), no_lge,'precision',10);
            dlmwrite(strcat(patient_sim_folder,'scalings\eikonal06_coarse_fib_',num2str(1),'_scar_',num2str(2),'_border_',num2str(i),'.csv'), light_lge,'precision',10);
         %   dlmwrite(strcat(patient_sim_folder,'scalings\eikonal06_coarse_fib_',num2str(2),'_scar.csv'),  int_lge,'precision',10);
          %  dlmwrite(strcat(patient_sim_folder,'scalings\eikonal06_coarse_fib_',num2str(3),'_scar.csv'),  adv_lge,'precision',10);
           %end
            end
   
   
   
       
 figure()
 headlight
 hold on
 patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',light_lge,'facecolor','interp')
colorbar 
 
[x,y,z] = ginput();


% and similarly for ZData in 3D plots
       
   ts_scalingRV=ones(21,1);
      ts_scalingRV_int=ones(21,1);
         ts_scalingRV_strong=ones(21,1);
         
    ts_scalingRV([12:21])=2;
    ts_scalingRV_int([12:21])=3;
    ts_scalingRV_strong([12:21])=4;
    
    lightRV=ts_scalingRV(mesh_cell_property);
    intRV=ts_scalingRV_int(mesh_cell_property);
    strongRV=ts_scalingRV_strong(mesh_cell_property);
    
    lightRV(H_mesh.rv)=1;
        lightRV(H_mesh.lv)=1;

    intRV(H_mesh.rv)=1;
        intRV(H_mesh.lv)=1;

    strongRV(H_mesh.rv)=1;
        strongRV(H_mesh.lv)=1;

    lightRV(EPI1.tri)=1;
    intRV(EPI1.tri)=1;
    strongRV(EPI1.tri)=1;

    
 figure()
 headlight
 hold on
 patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',intRV,'facecolor','interp')
colorbar 
axis off
title('conduction slowing factor','Fontsize',21)
    

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
    
      %  csvwrite(strcat(patient_sim_folder,'eikonal06_coarse_nodeScaling_',num2str(1),'f','.csv'),[ts_scaling(regions(:,1))']);

    dlmwrite(strcat(patient_sim_folder,'eikonal06_coarse_RV_sub_nodeScaling_',num2str(1),'f','.csv'), lightRV,'precision',10);
    dlmwrite(strcat(patient_sim_folder,'eikonal06_coarse_RV_sub_nodeScaling_',num2str(2),'f','.csv'),intRV,'precision',10);
    dlmwrite(strcat(patient_sim_folder,'eikonal06_coarse_RV_sub_nodeScaling_',num2str(3),'f','.csv'),strongRV,'precision',10);
   % end
      
        fid = fopen(strcat(patient_sim_folder,'regions_full_light_coarse','.txt'),'w');
        fprintf(fid,'%d %d %d \n',[regional_cond_lv(regions(:,1))' regional_cond_lv(regions(:,1))' regions].');
        fclose(fid)
    
    
   % roots=csvread(strcat(patient_sim_folder,'eikonal06_fine_rootNodes','.csv'));
         

%    
%    root_pos=[3.26,-9.915,-0.9604];
%    LV1_new_root = vtkClosestPoint( Mesh(H_mesh) , root_pos );
% 
%    
%       rootsNew=[roots(1:4),LV1_new_root];
% 
%    dlmwrite(strcat(patient_sim_folder,'eikonal06_fine_rootNodes_rv_ant_middle','.csv'),rootsNew)
%    
%    acts=dlmread(strcat(patient_sim_folder,'results\Eikonal\exp22\10.0x_3f_0.5fiber_roots_9_no_rv_pur_0\10.0x_3f_0.5fiber_roots_9_no_rv_pur_0_NEW_pred_ATMap.csv'));
%    
    figure()
    headlight
    hold on
        patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',mesh_cell_property(:,1),'facecolor','interp')
     %patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none')
    caxis([1 2])
    
    
    colorbar
    figure()
     plotMESH( MakeMesh( H_mesh , H_mesh.face ) , 'r','edgecolor','none' );
     headlight
     
     figure()
     plot(1:100,sqrt(cond(1:100)),'r','LineWidth',4)
     hold on
     plot(1:100,sqrt(cond5(1:100)),'b','LineWidth',4)
     plot(1:100,sqrt(cond6(1:100)),'m','LineWidth',4)
     set(gca,'Fontsize',45)
     xlabel('LGE Intensity [%]','Fontsize',45)
     ylabel('CV scaling factor ','Fontsize',45)
     legend({'low density fibrosis','intermediate density fibrosis','high density fibrosis'},'Fontsize',45)
%     camlight 
%     axis off
%     scatter3(H_mesh.xyz(roots,1),H_mesh.xyz(roots,2),H_mesh.xyz(roots,3),105,'filled','r')
%     axis off
%     view([-109 18])
%     
%      figure()
%     headlight
%     hold on
%     patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',acts,'facecolor','interp')
%     colorbar
%     camlight 
%     axis off
%     scatter3(H_mesh.xyz(rootsNew,1),H_mesh.xyz(rootsNew,2),H_mesh.xyz(rootsNew,3),105,'filled','r')
%     axis off
%     view([-109 18])