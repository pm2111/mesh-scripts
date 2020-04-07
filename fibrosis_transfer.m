
addpath(genpath('C:\Users\petnov\Dropbox\shared - copy\'));
addpath(genpath('C:\Users\petnov\Dropbox\mesh scripts'))
enableVTK;
pats={'06'};
    patient_nr=pats{1};
    %mkdir 'D:\ARVC meshing automatic\patients\patient01\regions';

    patient_sim_folder =  strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\');

    %H0 = read_VTK( fullfile( TORSO_MODEL_DIR , 'HEART_with_UPS_louie.vtk' ) ); %this is louie's heart!
    H_mesh = read_CHASTE( strcat(patient_sim_folder,    'chaste',patient_nr,'\HEART'));
    
    cell_regions=load(strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\','regions_labelled','.txt'));
    
      %septum stuff below, not necessary for every patient
          H1 = read_VTK(  strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\HEART.vtk')); %our target heart

        [EPI1,LV1,RV1,~,MAT1] = HEARTparts( H1 );
    septal_nodes=find(cell_regions==1);
    nearestlv=vtkClosestPoint(Mesh(LV1),H_mesh.xyz(septal_nodes,:)*10);
        nearestrv=vtkClosestPoint(Mesh(RV1),H_mesh.xyz(septal_nodes,:)*10);
        
        nearestEPI=vtkClosestPoint(Mesh(EPI1),H_mesh.xyz(septal_nodes,:)*10);
        distance= @(vec,vec1) ((vec(:,1)-vec1(:,1)).^2+(vec(:,2)-vec1(:,2)).^2+(vec(:,3)-vec1(:,3)).^2).^0.5;

distlv=distance(LV1.xyz(nearestlv,:),H_mesh.xyz(septal_nodes,:)*10);
distrv=distance(RV1.xyz(nearestrv,:),H_mesh.xyz(septal_nodes,:)*10);
distepi=distance(EPI1.xyz(nearestEPI,:),H_mesh.xyz(septal_nodes,:)*10);
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

% 
     septal_lv=septal_nodes(region_old==3);
     septal_rv=septal_nodes(region_old==4);
    cell_regions(septal_lv)=22;
        cell_regions(septal_rv)=23;

    figure()
    headlight
    hold on
    patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',cell_regions,'facecolor','interp')
    colorbar
    
    cond =@(x) -0.011*x+1.0;
        cond5= @(x) -0.0125*x+1.0;
    cond6= @(x) -0.014*x+1.0;
    cond7=@(x) -0.016*x+1.0;
    startRow = 2;

    [GadLV,GadRV] = read_csv('D:\sim3d\late_gad.csv', startRow);
     GadRV = str2double(GadRV);
     
     
    patient_nr_str=patient_nr;
    if patient_nr(1) =='0' ==0
        patient_nr=str2double(patient_nr);
    else
        patient_nr= str2double(patient_nr(2));

    end
    regional_cond_lv =[];
    CV_scaling=ones(23,1);
    sodium_cells = ones(size(H_mesh.xyz,1),1);
        for i=1:20
            j=i+20*(patient_nr-1);
            if isnan(str2double(GadLV{j})) ==0
                 regional_cond_lv(j-20*(patient_nr-1)+1) = cond5(str2double(GadLV{j}));
                 if(GadLV{j}==0)
                    CV_scaling(j-20*(patient_nr-1)+1) = 1;
                 else
                      CV_scaling(j-20*(patient_nr-1)+1) = sqrt(cond6(str2double(GadLV{j})));
                 end

            elseif GadRV(j) >0
                regional_cond_lv(j-20*(patient_nr-1)+1) = cond5(GadRV(j)*100/3);
                CV_scaling(j-20*(patient_nr-1)+1) =sqrt(cond6(GadRV(j)*100/3));

            else
                  regional_cond_lv(j-20*(patient_nr-1)+1)=1.000;
                  CV_scaling(j-20*(patient_nr-1)+1)=1.000;

            end
            if isnan( regional_cond_lv(j-20*(patient_nr-1))+1) ==1
                regional_cond_lv(j-20*(patient_nr-1)+1)=1.0;
                CV_scaling(j-20*(patient_nr-1)+1)=1.0;
                
            end
            regional_cond_lv(1) =1; %septum late gad info is not present
            CV_scaling(1) =1; %septum late gad info is not present

    %         fid = fopen(strcat('D:\ARVC meshing automatic\patients\patient05\conductivities\',num2str(Label{j}),'.txt'),'w');
    %         fprintf(fid,'%d',regional_cond_lv(i));
    %         fclose(fid)
        end
            patient_nr=pats{1};

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
    fields=fields(mapping);
   % CV_scaling(find(strcmp(fields,'mid_inferolateral_rv')))=0.5;
    %    CV_scaling(find(strcmp(fields,'basal_inferolateral_rv')))=0.5;

%       CV_scaling(22)=0.55;
         CV_scaling(23)=1.0;%rv septum
         CV_scaling(22)=0.50;%lv septum

    ts_scaling=1.0./CV_scaling;
    
dlmwrite(strcat(patient_sim_folder,'eikonal06_fine_nodeScaling_1f','.csv'),[ts_scaling(cell_regions(:,1))'].');
 figure()
    headlight
    hold on
    patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',ts_scaling(cell_regions(:,1)),'facecolor','interp')
    colorbar
    
     figure()
    headlight
    hold on
    patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',cell_regions(:,1),'facecolor','interp')
    colorbar
    