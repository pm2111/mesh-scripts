close all;
addpath(genpath('C:/Users/petnov/Dropbox/shared - copy'));

enableVTK;
sim_name='controlAT';
pat_nr='09';
folders=dir(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,',pat_nr',sim_name,'\'));

FBfine =transform(read_CHASTE(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\chaste',pat_nr,'\HEART')),'s',1/1);
H.face      = MeshBoundary( FBfine.tri );
H.triCENTER = meshFacesCenter( FBfine );
H.triVOL    = meshQuality( FBfine , 'volume' );
H.xyzIDS    = vtkClosestPoint( Mesh( FBfine ) , FBfine.xyz );
H.Gop       = meshGradient( FBfine );
H.tri =FBfine.tri;
H.xyz=FBfine.xyz;

FBfine=struct;
%%
electrodes = load(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\mpp\ECG_ELECTRODES.mat'));


E = [ electrodes.ELECTRODES.LA ; electrodes.ELECTRODES.RA ; electrodes.ELECTRODES.LL ; electrodes.ELECTRODES.RL ; electrodes.ELECTRODES.V1 ; electrodes.ELECTRODES.V2 ; electrodes.ELECTRODES.V3 ; electrodes.ELECTRODES.V4 ; electrodes.ELECTRODES.V5 ; electrodes.ELECTRODES.V6 ];
E=E/10; %chaste is in cm

%plotMESH( MakeMesh( H , H.face ) , 'r' ); hplotMESH( T , 'nf' );
%hplot3d( E , '*m' ,'LineWidth',4,'MarkerSize',10 );
%axis(objbounds)
%instead of using positions from torso use positions of body surface pot
%recordings:
% locs=load('D:\Data_ARVC_raw_electrodes_vest\Data_ARVC\ARVC_electrode_locations');
% locs=locs.arvc9_electrodes/10;
% idx=[];
% for i=1:10
%     idx(i)=knnsearch(locs,E(i,:));
% end
% 
% E_compute_ecg=locs(idx,:);
% 
% 
% E_left=locs(idx,:);
% E_left(:,1)=locs(idx,1)-2.5;
% 
% E_compute_ecg=locs(idx,:);

%%save the voltages into 1ms bins

%fileinfo = hdf5info('D:/ARVC meshing automatic/patients/patient18/results/TestHeart18ControlECG/results.h5');
permutation = load(strcat('D:/ARVC meshing automatic/patients/patient',pat_nr,'/results/',sim_name,'/permutation.txt'));
permutation = permutation+1;
N=size(H.xyz,1);
ecg_full=[zeros(1,12)];
batch_size=5;
t_end=600;
nr_batches=int16(t_end/batch_size)+1;
for i=1:nr_batches
    for n=1:batch_size
          id=(i-1)*batch_size+n;
           % voltages = h5read('TestHeart05FullBeatandActivation\results.h5',fileinfo.GroupHierarchy.Datasets(1).Name,[1,1,t] ,[1,N,1]);
                  % inf = h5info(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',sim_name,'\results.h5'));

           voltages = h5read(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',sim_name,'\results.h5'),'/Data',[1, 1, double(id)], [1, N ,1] );

            voltages_heart = voltages(permutation); %select only the heart activation times!! 
          disp(id);

          %gradV(t}=reshape( H.Gop * getv( voltages_heart' , H.xyzIDS ) ,[],3);
          H.triG(:,:,i) = reshape( H.Gop * getv( voltages_heart' , H.xyzIDS ) ,[],3);

          %gradV(t}=reshape( H.Gop * getv( voltages_heart' , H.xyzIDS ) ,[],3);
    end
    [~,I,II,III,aVR,aVL,aVF,V1,V2,V3,V4,V5,V6] = computePECG(  H , E );
    ecg_full = vertcat(ecg_full,I,II,III,aVR,aVL,aVF,V1,V2,V3,V4,V5,V6);
    H=rmfield(H,'triG');
end
I_real =[];
II_real= [];
III_real= [];
aVR_real= [];
aVL_real= [];
aVF_real=[];
V1_real= [];
V2_real= [];
V3_real=[];
V4_real=[];
V5_real=[];
V6_real= [];
for i=1:nr_batches
    I_real= vertcat(I_real,ecg_full(:,1+(i-1)*12));
    II_real= vertcat(II_real,ecg_full(:,2+(i-1)*12));
    III_real= vertcat(III_real,ecg_full(:,3+(i-1)*12));
    aVR_real= vertcat(aVR_real,ecg_full(:,4+(i-1)*12));
    aVL_real= vertcat(aVL_real,ecg_full(:,5+(i-1)*12));
    aVF_real= vertcat(aVF_real,ecg_full(:,6+(i-1)*12));
    V1_real= vertcat(V1_real,ecg_full(:,7+(i-1)*12));
    V2_real= vertcat(V2_real,ecg_full(:,8+(i-1)*12));
    V3_real= vertcat(V3_real,ecg_full(:,9+(i-1)*12));
    V4_real= vertcat(V4_real,ecg_full(:,10+(i-1)*12));
    V5_real= vertcat(V5_real,ecg_full(:,11+(i-1)*12));
    V6_real= vertcat(V6_real,ecg_full(:,12+(i-1)*12));
end


save(strcat('D:\ARVC meshing automatic\patients\patient',pat_nr,'\results\',sim_name,'\ecg.mat'),'I_real','II_real','III_real','aVR_real','aVL_real','aVF_real','V1_real','V2_real','V3_real','V4_real','V5_real','V6_real');




%%
