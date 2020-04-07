addpath(genpath('C:\users\petnov\Dropbox\'))
mesh_cell_property= load('D:\ARVC meshing automatic\patients\patient06\regions_full_light_coarse.txt');
H_mesh = read_CHASTE('D:\ARVC meshing automatic\patients\patient06\chaste06_coarse\HEART');
patient_folder ='D:\ARVC meshing automatic\patients\patient06\';
%now save values for the conductivities scaling factors, assuming a linear
%rel between LGE intensity and conductivity of tissue (max CV reduction at
%50% LGEwhere CV is halved (conductivity is divided by 4).Assume
%conductivities change equally in the 3 fiber directions

       cond =@(x) -0.011*x+1.0;
        cond2= @(x) -0.0125*x+1.0;
    cond3= @(x) -0.014*x+1.0;


[GadLV,GadRV] = read_csv('D:\sim3d\late_gad.csv', startRow);
 GadRV = str2double(GadRV);
 GadRV(92) = 73;
 GadRV(94) = 76;
 GadRV(96) = 70;
 GadRV(99) = 73;

patient_nr=6;
regional_cond_lv =[];
    for i=1:20
        j=i+20*(patient_nr-1);
        if isnan(str2double(GadLV{j})) ==0
             regional_cond_lv(j-20*(patient_nr-1)+1) = cond3(str2double(GadLV{j}));
        elseif isnan(GadRV(j)) ==0
            regional_cond_lv(j-20*(patient_nr-1)+1) = cond3(GadRV(j));
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
    
    
    strong_LGE= regional_cond_lv(mesh_cell_property(:,3));
    
      figure()
    headlight
    hold on
    patch('vertices',H_mesh.xyz,'faces',H_mesh.tri,'edgecolor','none','FaceVertexCData',strong_LGE','facecolor','interp')
    colorbar
    camlight 
    axis equal

%save the coords of regions in the same file and save a number
%corresponding to the region number next to it
    fid = fopen(strcat(patient_folder,'regions_full_strong_small','.txt'),'w');
        fprintf(fid,'%d %d \n',[regional_cond_lv(mesh_cell_property(:,3))' regional_cond_lv(mesh_cell_property(:,3))'].');
    fclose(fid)
    
