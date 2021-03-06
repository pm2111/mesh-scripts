%explore the fibre generator!

%LINUX filepaths
addpath(genpath('/users/petnov/Dropbox/'));

torso_louie = read_CHASTE('/data/sim3d/louie/data/nejib_20131010/FullMeshVolume');


ortho = dlmread('/data/sim3d/louie/data/nejib_20131010/FullMeshVolume.ortho', ' ', 1, 0);

torso_peter =  read_CHASTE('/data/sim3d/patient02/TORSO');


%windows path

addpath(genpath('D:/shared/'));

torso = read_CHASTE('D:\ARVC meshing automatic\patients\patient02\chaste_ernesto\TORSO');

ortho = dlmread('D:\ARVC meshing automatic\patients\patient02\chaste_ernesto\HEART.ortho',' ', [1 0 ]);


x = randi(10,3);

x = magic(3);
dlmwrite('D:\sim3d\write\myfile.txt',x,'delimiter',' ', 'newline','pc');

fid = fopen('D:\ARVC meshing automatic\patients\patient02\chaste_ernesto\HEART.ortho');

nr_fiber_rows = fgetl(fid); %this is a useful function to read a file line by line! 

fid1 = fopen('/data/sim3d/patient02/TORSO.ortho');

nr_torso_rows = fgetl(fid1);

fid2 = fopen('/data/sim3d/patient02/ortho_from_heart');

nr_torso_rows = fgetl(fid2);

nr_heart_fibre_lines = str2double(nr_torso_rows);

fid2 = fopen('/data/sim3d/patient02/TORSO.ortho');
indices_to_check = [];
for i=1:60422161 %nr_heart_fibre_lines
    line = fgetl(fid2);
    if size(line) < 160
        indices_to_check = [indices_to_check,i]
    end
end
        
    


%need to assign a fiber orientation for each element in the torso

total_nodes = size(torso.tri,1);
%need to append the difference between the total number of elements and the
%number of fiber rows currently generated. We are only generating for the
%heart atm.

rows_to_append = total_nodes- str2num(nr_fiber_rows);
26959113
what_to_append = [1, 0, 0 , 0, 1, 0, 0, 0, 1];
dlmwrite('Z:\write\what to append.txt',what_to_append,'delimiter',' ', 'newline','pc','roffset',0,'precision','%.6f')
60422160

%use this following line to append to the chaste.ortho file (need as many
%lines as there are elements (tetras)!)
dlmwrite('Z:\write\chaste_original.ortho',what_to_append,'-append','delimiter',' ', 'newline','pc','roffset',0,'precision','%.6f')

for i=1:rows_to_append
    dlmwrite('D:\ARVC meshing automatic\patients\patient02\chaste_ernesto\HEARTappended.ortho',what_to_append,'-append','delimiter',' ', 'newline','pc','roffset',0,'precision','%.6f')
end
