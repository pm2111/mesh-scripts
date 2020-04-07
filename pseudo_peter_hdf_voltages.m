clear all;
close all;
addpath(genpath('C:/Users/petnov/Dropbox/shared - copy'));

enableVTK;
sim_name='TestHeart06NewRootsOldFibers';
pat_nr='06';

FBfine =transform(read_CHASTE(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\chaste',pat_nr,'\HEART')),'s',1/1);
H.face      = MeshBoundary( FBfine.tri );
H.triCENTER = meshFacesCenter( FBfine );
H.triVOL    = meshQuality( FBfine , 'volume' );
H.xyzIDS    = vtkClosestPoint( Mesh( FBfine ) , FBfine.xyz );
H.Gop       = meshGradient( FBfine );
H.tri =FBfine.tri;
H.xyz=FBfine.xyz;

mesh=struct();
mesh.xyz=H.xyz;
mesh.tri=H.tri;

%%
electrodes = load(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\mpp\ECG_ELECTRODES.mat'));


E = [ electrodes.ELECTRODES.LA ; electrodes.ELECTRODES.RA ; electrodes.ELECTRODES.LL ; electrodes.ELECTRODES.RL ; electrodes.ELECTRODES.V1 ; electrodes.ELECTRODES.V2 ; electrodes.ELECTRODES.V3 ; electrodes.ELECTRODES.V4 ; electrodes.ELECTRODES.V5 ; electrodes.ELECTRODES.V6 ];
E=E/10; %chaste is in cm

%plotMESH( MakeMesh( H , H.face ) , 'r' ); hplotMESH( T , 'nf' );
%hplot3d( E , '*m' ,'LineWidth',4,'MarkerSize',10 );
%axis(objbounds)
%instead of using positions from torso use positions of body surface pot
%recordings:
locs=load('D:\Data_ARVC_raw_electrodes_vest\Data_ARVC\ARVC_electrode_locations');
locs=locs.arvc6_electrodes/10;
idx=[];
for i=1:10
    idx(i)=knnsearch(locs,E(i,:));
end

E_compute_ecg=locs(idx,:);


E_left=locs(idx,:);
E_left(:,1)=locs(idx,1)-2.5;

E_compute_ecg=locs(idx,:);

%%save the voltages into 1ms bins

%fileinfo = hdf5info('D:/ARVC meshing automatic/patients/patient18/results/TestHeart18ControlECG/results.h5');
permutation = load(strcat('D:/ARVC meshing automatic/patients/patient',pat_nr,'/results/',sim_name,'/permutation.txt'));
permutation = permutation+1;
N=size(FBfine.xyz,1);

ecg_full=[];
for n=0:399

           % voltages = h5read('TestHeart05FullBeatandActivation\results.h5',fileinfo.GroupHierarchy.Datasets(1).Name,[1,1,t] ,[1,N,1]);
            fid=fopen(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',sim_name,'\voltages\voltages.bin',int2str(800+n)),'rb');
           voltages = fread(fid,'double');

            voltages_heart = voltages(permutation); %select only the heart activation times!! 
          disp(n);
          fclose(fid);

          %gradV(t}=reshape( H.Gop * getv( voltages_heart' , H.xyzIDS ) ,[],3);
          H.triG(:,:,1) = reshape( H.Gop * getv( voltages_heart , H.xyzIDS ) ,[],3);
    
    [~,I,II,III,aVR,aVL,aVF,V1,V2,V3,V4,V5,V6] = computePECG(  H , E );
     int=[I(end),II(end),III(end),aVR(end),aVL(end),aVF(end),V1(end),V2(end),V3(end),V4(end),V5(end),V6(end)];

    ecg_full = [ecg_full,int];
    H=rmfield(H,'triG');
end


save(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',sim_name,'\ecg.mat'),'ecg_full');




%%
