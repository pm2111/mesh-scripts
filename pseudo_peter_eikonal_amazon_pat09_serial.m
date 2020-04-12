addpath(genpath('C:/Users/petnov/Dropbox/shared - copy'));
addpath(genpath('C:/Users/petnov/Dropbox/mesh scripts/'));

mesh_path='D:\ARVC meshing automatic\patients\patient09\eikonal09\';


sim_path='D:\peter_results\root_calibration\';
ecg_path='D:\peter_results\root_calibration_ECGs\';


pat_nr='09';

FBfine=struct;
FBfine.xyz=csvread(strcat(mesh_path,'eikonal09_fine_xyz.csv'));
FBfine.tri=csvread(strcat(mesh_path,'eikonal09_fine_tri.csv'));
FBfine.triCENTER = meshFacesCenter( FBfine );
FBfine.triVOL    = meshQuality( FBfine , 'volume' );
FBfine.xyzIDS    = 1:size(FBfine.xyz,1);
FBfine.Gop       = meshGradient( FBfine );


   disp(strcat("i have loaded the mesh"));

% rundir='/home/ubuntu/Eikonal_peter/';
% runDir='/home/ubuntu/Eikonal_peter/';


%%
electrodes = load(strcat('D:\ARVC meshing automatic\patients\patient09\mpp\','ECG_ELECTRODES.mat'));


E = [ electrodes.ELECTRODES.LA ; electrodes.ELECTRODES.RA ; electrodes.ELECTRODES.LL ; electrodes.ELECTRODES.RL ; electrodes.ELECTRODES.V1 ; electrodes.ELECTRODES.V2 ; electrodes.ELECTRODES.V3 ; electrodes.ELECTRODES.V4 ; electrodes.ELECTRODES.V5 ; electrodes.ELECTRODES.V6 ];
E=E/10; %chaste is in cm


%body surface electrode locs



idx=[ 117    21   107    11    60    92    80    69   112   102];
N=size(FBfine.xyz,1);

folders=dir(sim_path);
%gcp=parpool(1);
%addAttachedFiles(gcp,{'parsave.m','computePECG.m'});

for j=3:size(folders,1)
    names=dir(strcat(sim_path,folders(j).name));
    
    acts=dlmread(strcat(sim_path,folders(j).name,'\',names(3).name));
    batch_size=1;
   % Gop=[];
   % triCENTER=[];
  %  triVOL=[];
 %   xyzIDS=[];    
    t_end=max(acts);
    tot_batches= int16(t_end/batch_size)+1;
    ecg_full=[zeros(1,12)];
    voltages_heart=ones(size(acts))*-89.0;
    for i=1:tot_batches
        for n=1:batch_size
              id=(i-1)*batch_size+n;
              idx=find(acts==id);
              voltages_heart(idx)=40.0;
              disp(id);
              FBfine.triG(:,:,n) = reshape( FBfine.Gop * getv( voltages_heart , FBfine.xyzIDS ) ,[],3);
        end
        [~,I,II,III,aVR,aVL,aVF,V1,V2,V3,V4,V5,V6] = computePECG(  FBfine ,E );
        int=[I,II,III,aVR,aVL,aVF,V1,V2,V3,V4,V5,V6];
        ecg_full = vertcat(ecg_full,int);
        Hbatch.triG=[];
    end
    save(strcat(ecg_path,names(3).name,'pseudo_ECG.mat'),'ecg_full');
   disp(strcat("i am moving onto ECG #",num2str(j-2)));

end

