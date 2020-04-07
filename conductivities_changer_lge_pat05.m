addpath(genpath('C:\users\petnov\Dropbox\'))
mesh_cell_property= load('D:\ARVC meshing automatic\patients\patient05\regions_pat05_in_mesh_vtk_ordering.txt');
H_mesh = read_CHASTE('D:\ARVC meshing automatic\patients\patient05\chaste_gmsh\HEART');
%now save values for the conductivities scaling factors, assuming a linear
%rel between LGE intensity and conductivity of tissue (max CV reduction at
%50% LGEwhere CV is halved (conductivity is divided by 4).Assume
%conductivities change equally in the 3 fiber directions
cond =@(x) -0.0095*x+1.0;

%this second function was calibrated to explain the conduction delay in the
%LV 
cond2= @(x) -0.0055*x+1.0;

cond3= @(x) -0.011*x+1.0;

startRow = 2;

[GadLV,GadRV] = read_csv('D:\sim3d\late_gad.csv', startRow);
 GadRV = str2double(GadRV);
 GadRV(92) = 67;
 GadRV(94) = 70;
 GadRV(96) = 63;
 GadRV(99) = 67;

patient_nr=5;
regional_cond_lv =[];
lge_percentage = [];
    for i=1:20
        j=i+20*(patient_nr-1);
        if isnan(str2double(GadLV{j})) ==0
             regional_cond_lv(j-20*(patient_nr-1)+1) = cond3(str2double(GadLV{j}));
             lge_percentage(j-20*(patient_nr-1)+1) = str2double(GadLV{j});
        elseif isnan(GadRV(j)) ==0
            regional_cond_lv(j-20*(patient_nr-1)+1) = cond3(GadRV(j));
            lge_percentage(j-20*(patient_nr-1)+1) = GadRV(j);
               
        else
              regional_cond_lv(j-20*(patient_nr-1)+1)=1.000;
              lge_percentage(j-20*(patient_nr-1)+1) =0;

        end
        if isnan( regional_cond_lv(j-20*(patient_nr-1))+1) ==1
            regional_cond_lv(j-20*(patient_nr-1)+1)=1.0;
        end
        regional_cond_lv(1) =1; %septum late gad info is not present
%         fid = fopen(strcat('D:\ARVC meshing automatic\patients\patient05\conductivities\',num2str(Label{j}),'.txt'),'w');
%         fprintf(fid,'%d',regional_cond_lv(i));
%         fclose(fid)
    end

regional_ina_scaling = ones(size(regional_cond_lv,2),1);
regional_ina_scaling(3) =0.8;
regional_ina_scaling(5) =0.8;
regional_ina_scaling(10) =0.8;
regional_ina_scaling(11) =0.8;
regional_ina_scaling(12:end) = 0.8;
regional_ina_scaling = regional_ina_scaling';
    
%save the coords of regions in the same file and save a number
%corresponding to the region number next to it
    fid = fopen(strcat('D:\ARVC meshing automatic\patients\patient05\','regions_full','.txt'),'w');
    fprintf(fid,'%d %d \n',[regional_cond_lv(mesh_cell_property)' regional_ina_scaling(mesh_cell_property)'].');
    fclose(fid)
    
    H_mesh.xyzFibrosis = lge_percentage(mesh_cell_property)';
    
    write_VTK(H_mesh,'D:\ARVC meshing automatic\patients\patient05\fibrosis_lge_sim_mesh.vtk','binary')
    
    
