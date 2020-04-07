addpath(genpath('C:\Users\petnov\Dropbox\shared\'));
addpath(genpath('C:\Users\petnov\Dropbox\'));

TORSO_MODEL_DIR = 'D:\ARVC meshing automatic\patients\patient05\chaste_ernesto\';

%H0 = read_VTK( fullfile( TORSO_MODEL_DIR , 'HEART_with_UPS_louie.vtk' ) ); %this is louie's heart!

H0 = read_VTK('D:\ARVC meshing automatic\patients\patient05\ecgi_clipped.vtk');

%chaste works in cm, our meshes are built in mm (BE CAREFUL WHEN LOCATING
%INDICES )

LV0_roots = [ 157258 , 129516 , 140114 , 131341 ];
RV0_roots = [ 198857 , 209309 , 186671 ];



figure()
plotMESH( H0 , 'ne' );
hplot3d( H0.xyz( LV0_roots ,:) , 'o1kb20' );
hplot3d( H0.xyz( RV0_roots ,:) , 'o1kr20' );
headlight


%%

H1 = read_VTK(  'D:\ARVC meshing automatic\patients\patient05\HEART.vtk'); %our target heart

plotMESH( H0 ,'ne');
hplotMESH( H1 , 'ne' ,'FaceColor',[1 1 1]*0.6)
headlight; axis(objbounds);

%HD_loaded = load('Z:\Dropbox\patients\HD.mat');

%%

u = loadv( 'D:\ARVC meshing automatic\patients\patient05\mpp\Hmapping' , 'u' ); %transfer matrix u
%this is a function converting coordinates from H0 to H1
%it is generated by running the mpp_Heart_Mapping Function of Ernesto s
%pipeline!!

if 0
HD = H0;
HD.xyz = u( HD.xyz );

plotMESH( HD ,'ne');
hplotMESH( H1 , 'ne' ,'FaceColor',[1 1 1]*0.6)
headlight; axis(objbounds);
axis off
end


%%

LV1_roots = vtkClosestPoint( Mesh(H1) , u( H0.xyz( LV0_roots ,:) ) );
RV1_roots = vtkClosestPoint( Mesh(H1) , u( H0.xyz( RV0_roots ,:) ) );


LV1_roots_chaste  = LV1_roots-1; %convert to chaste indexing
RV1_roots_chaste = RV1_roots-1;

fid = fopen('D:\ARVC meshing automatic\patients\patient05\root nodes\start.txt','w');
fprintf(fid,'%u\t',vertcat(LV1_roots_chaste,RV1_roots_chaste))

%check that these nodes lie on the endocardium:
fileID = fopen('D:\ARVC meshing automatic\patients\patient02\chaste_ernesto\rv.edges','r');
rv_nodes = fread(fileID);
%check that root nodes lie on endo surface
root_liesRV = [];
for i=1:3
    root_liesRV(i) = find(rv_edge_nodes == RV1_roots(i));
end

root_liesLV = [];
for i=1:4
    root_liesLV(i) = find(lv_edge_nodes == LV1_roots(i));
end


joined_roots = vertcat(LV1_roots,RV1_roots);



permutation = load('D:\sim3d\patient02\permutation.txt');
%the rows in the permutation file correspond to the index in chaste and t
%the number in the cell corresponds to the vtk ordering!

for i=1:size(joined_roots)
    roots_chaste(i) = find(permutation ==joined_roots(i));
end

%roots_chaste = [2021939,2559829,2006872,2210285,2930394,2649272,73922];
figure()
plotMESH( H1 , 'ne' );
hplot3d( H1.xyz( LV1_roots ,:) , 'o1kb20' );
hplot3d( H1.xyz( RV1_roots ,:) , 'o1kr20' );
headlight;