
addpath(genpath('C:\Users\petnov\Dropbox\shared - copy\'));
addpath(genpath('C:\Users\petnov\Dropbox\mesh scripts'))
addpath(genpath('C:\Users\petnov\Dropbox\shared\IO\'));

enableVTK;
pats={'09'};
    patient_nr=pats{1};
    patient_sim_folder=strcat('D:\ARVC meshing automatic\patients\patient',pats{1},'\');
        H_mesh = read_CHASTE( strcat(patient_sim_folder,    'chaste',patient_nr,'\HEART'));
        
        
sim_names= {'TestHeart09Control','TestHeart09LGELight','TestHeart09LGEIntermediate','TestHeart09LGEStrong'};
    permutation = load(strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\results\',sim_names{1},'\permutation.txt'));
    permutation = permutation+1;
 
upstroke={};
for i =1:4
    try
    fileinfo = hdf5info(strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\results\',sim_names{i},'\results_maps.h5'));

    upstroke{i} = h5read(strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\results\',sim_names{i},'\results_maps.h5'),fileinfo.GroupHierarchy.Datasets(1).Name);
        permutation = load(strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\results\',sim_names{i},'\permutation.txt'));
    permutation = permutation+1;
        upstroke{i}=upstroke{i}(permutation);
    catch
       fileinfo = hdf5info(strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\results\',sim_names{i},'\results.h5'));

    upstroke{i} = h5read(strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\results\',sim_names{i},'\results.h5'),fileinfo.GroupHierarchy.Datasets(3).Name);
   permutation = load(strcat('D:\ARVC meshing automatic\patients\patient',patient_nr,'\results\',sim_names{i},'\permutation.txt'));
    permutation = permutation+1;
        upstroke{i}=upstroke{i}(permutation);
    end
     
end

    figure()

    labels={'no LGE', 'LGE max CV scaling factor 0.58','LGE max CV scaling factor 0.5','LGE max CV scaling factor 0.4'};
for i=1:4
subplot_tight(2,2,i)
    headlight
    hold on
    patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',[upstroke{i}],'facecolor','interp')
    axis off
    caxis([0 60])
    title(strcat('activation times [ms] ', labels{i}),'FontSize',25)
    colormap(flipud(jet))
end
colorbar('southoutside')
set(gca,'Fontsize',25)