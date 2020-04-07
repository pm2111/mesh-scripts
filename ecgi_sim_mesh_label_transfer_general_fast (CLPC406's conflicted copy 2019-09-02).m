
addpath(genpath('C:\Users\petnov\Dropbox\shared - copy\'));
addpath(genpath('C:\Users\petnov\Dropbox\mesh scripts'))
enableVTK;
pats={'06'};
    patient_nr=pats{1};
    %mkdir 'D:\ARVC meshing automatic\patients\patient01\regions';

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
        
    regions= load(strcat(patient_sim_folder,'regions_labelled.txt'));
    
    figure()
    headlight
    hold on
    patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',regions,'facecolor','interp')
    colorbar

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

                 else
                      CV_scaling(j-20*(patient_nr-1)+1) = sqrt(cond(str2double(GadLV{j})));
                      CV_scaling_int(j-20*(patient_nr-1)+1) = sqrt(cond5(str2double(GadLV{j})));
                      CV_scaling_strong(j-20*(patient_nr-1)+1) = sqrt(cond6(str2double(GadLV{j})));
                      CV_scaling_vstrong(j-20*(patient_nr-1)+1) = sqrt(cond7(str2double(GadLV{j})));

                 end

            elseif GadRV(j) >0
                CV_scaling(j-20*(patient_nr-1)+1) =sqrt(cond(GadRV(j)*100/3));
                CV_scaling_int(j-20*(patient_nr-1)+1) =sqrt(cond5(GadRV(j)*100/3));
                CV_scaling_strong(j-20*(patient_nr-1)+1) =sqrt(cond6(GadRV(j)*100/3));
                CV_scaling_vstrong(j-20*(patient_nr-1)+1) =sqrt(cond7(GadRV(j)*100/3));

            else
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

    ts_scaling=1.0./CV_scaling;
    ts_scaling_int=1.0./CV_scaling_int;
    ts_scaling_strong=1.0./CV_scaling_strong;
    ts_scaling_vstrong=1.0./real(CV_scaling_vstrong);

    dlmwrite(strcat(patient_sim_folder,'eikonal06_fine_nodeScaling_1f','.csv'),[ts_scaling(regions(:,1))'].');
    dlmwrite(strcat(patient_sim_folder,'eikonal06_fine_nodeScaling_2f','.csv'),[ts_scaling_int(regions(:,1))'].');
    dlmwrite(strcat(patient_sim_folder,'eikonal06_fine_nodeScaling_3f','.csv'),[ts_scaling_strong(regions(:,1))'].');
    dlmwrite(strcat(patient_sim_folder,'eikonal06_fine_nodeScaling_4f','.csv'),[ts_scaling_vstrong(regions(:,1))'].');

    roots=csvread(strcat(patient_sim_folder,'eikonal06_fine_rootNodes','.csv'));
         
   roots=[roots(1:3),roots([6])];
   
   dlmwrite(strcat(patient_sim_folder,'eikonal06_fine_rootNodes_fib','.csv'),roots)
    figure()
    headlight
    hold on
    patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',ts_scaling(regions(:,1)),'facecolor','interp')
    colorbar
    scatter3(H_mesh.xyz(roots,1),H_mesh.xyz(roots,2),H_mesh.xyz(roots,3),105,'filled','r')
    axis off