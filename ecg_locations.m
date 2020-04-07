addpath(genpath('C:/Users/petnov/Dropbox/shared - copy/'));
addpath(genpath('C:/Users/petnov/Dropbox/'));

enableVTK;
sim_name='TestControlECGPat09';
pat_nr='06';

FBfine =read_CHASTE(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\chaste',pat_nr,'\HEART'));

H.face      = MeshBoundary( FBfine.tri );
H.tri =FBfine.tri;
H.xyz=FBfine.xyz;


electrodes = load(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\mpp\ECG_ELECTRODES.mat'));

E = [ electrodes.ELECTRODES.LA ; electrodes.ELECTRODES.RA ; electrodes.ELECTRODES.LL ; electrodes.ELECTRODES.RL ; electrodes.ELECTRODES.V1 ; electrodes.ELECTRODES.V2 ; electrodes.ELECTRODES.V3 ; electrodes.ELECTRODES.V4 ; electrodes.ELECTRODES.V5 ; electrodes.ELECTRODES.V6 ];
E=E/10; %chaste is in cm


%body surface electrode locs

locs=load('D:\Data_ARVC_raw_electrodes_vest\Data_ARVC\ARVC_electrode_locations');
locs=locs.arvc6_electrodes/10;
for i=1:10
    idx(i)=knnsearch(locs,E(i,:));
end
 117    32   107    22    60    92    80    69   112   102
idx_2=[ 117    21   107    11    60    92    80    69   112   102];

figure()
%plotMESH( MakeMesh( H_mesh , H_mesh.face ) , 'r' ); hplotMESH( T.VEST , 'nf' );
 patch('vertices',H_mesh.xyz,'faces',H_mesh.face,'edgecolor','none','FaceVertexCData',[upstroke{1}],'facecolor','interp')

hplot3d( E , '*m' ,'LineWidth',4,'MarkerSize',10 );
%hplot3d( locs(idx,:) , '*g' ,'LineWidth',4,'MarkerSize',10 );
headlight
axis equal
caxis([0 18])
hold on
%scatter3(E(idx,1),locs(idx,2),locs(idx,3),'filled','r')