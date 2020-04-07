addpath(genpath('C:\users\petnov\Dropbox\shared - copy'));
enableVTK;
cd D:\ARVC' meshing automatic'\patients

pat_nr='06';
sim_name='TestHeart06Disrupted';
control_sim= 'TestHeart01Control';

%mesh_vtk = read_VTK('/data/sim3d/louie/HEART_with_UPS.vtk'); 
mesh = read_CHASTE(strcat('patient',pat_nr,'\chaste',pat_nr,'\HEART'));

    fileinfo = hdf5info(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',sim_name,'\results.h5'));
    permutation = load(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',sim_name,'\permutation.txt'));
    permutation = permutation+1;
    try 
            upstroke = csvread(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',sim_name,'\eikonal09_fine_NEW_pred_ATMap_150_1x_3f_NormalRootsLightFibLeftSeptumAnteriorIntermediate.csv'));
    catch
        
        upstroke = h5read(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',sim_name,'\results.h5'),fileinfo.GroupHierarchy.Datasets(3).Name);
    end
    fib=load(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\regions_full_light_lv_sept_int.txt'));
    fib2=csvread(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\eikonal09_fine_nodeScaling_2f.csv'));
   
    csvwrite(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\','Eikonal','\ActivationTimesStrongLGEBid.csv'),upstroke(1,permutation,1));

    csvwrite(  strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',sim_name,'\eikonal06_fine_true_ATMap_150_1x_3f.csv'),upstroke(1,permutation,1));
    fileinfo2 = hdf5info(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',control_sim,'\results_maps.h5'));

    upstroke_control = h5read(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',control_sim,'\results_maps.h5'),fileinfo2.GroupHierarchy.Datasets(1).Name);
   permutation_control = load(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',control_sim,'\permutation.txt'));
    permutation_control=permutation_control+1;
       
    upstroke_heart1 = upstroke(1,permutation,1); %select only the heart activation times!! 
    
    rms_error=sum((acts-upstroke_heart1').^2)/size(acts,1);
    
    
    
    upstroke_heart2 = upstroke(1,permutation,2); %select only the heart activation times!! 

    upstroke_control1=upstroke_control(1,permutation_control,1); %select only the heart activation times!! 

    %upstroke_control2=upstroke_control(1,permutation,2); %select only the heart activation times!! 

    N=size(mesh.xyz,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    mesh.xyzUpstroke = transpose(upstroke_heart1);
    mesh.xyzUpstroke2 = transpose(upstroke_heart2);

    %mesh.xyzUpstrokeControl = transpose(upstroke_control2);

    %compare results with control values for this patient
    % control_res = 'D:\ARVC meshing automatic\patients\patient05\results\TestHeart05FullBeatandActivation\';
    % upstroke_control = h5read(strcat(control_res,'results.h5'),fileinfo.GroupHierarchy.Datasets(3).Name);
    % upstroke_control = upstroke_control(1,permutation,1);
     mesh.xyzFib=fib(:,1);
          mesh.xyzFib2=fib2;

    % mesh.xyzUpstroke2 = transpose(upstroke_heart2);
    % 
    % mesh.xyzUpstrokeDifference = transpose(upstroke_heart- upstroke_control);
    % upstroke_diff= transpose(upstroke_heart- upstroke_control);
    mesh.celltype =10;
    mesh.xyzeikonal=upstroke;
    write_VTK(mesh,strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',sim_name,'\results.vtk'),'binary');
