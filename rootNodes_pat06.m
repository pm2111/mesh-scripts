addpath(genpath('C:\Users\petnov\Dropbox\shared - copy'));
enableVTK;
TORSO_MODEL_DIR = 'D:\ARVC meshing automatic\patients\patient06\chaste06\';
SUBJECT_DIR = 'D:\ARVC meshing automatic\patients\patient06\';
%H0 = read_VTK( fullfile( TORSO_MODEL_DIR , 'HEART_with_UPS_louie.vtk' ) ); %this is louie's heart!
vest=load(strcat(SUBJECT_DIR,'mpp\','fittedVEST.mat'));
H0 = read_VTK('C:\Users\petnov\Dropbox\shared - copy\MeshPersonalizationPipeline\TORSO\HEART.vtk');

%chaste works in cm, our meshes are built in mm (BE CAREFUL WHEN LOCATING
%INDICES OF ROOTS, always compare the same type of file eg vtk generated by
%pipeline or convert units

LV0_roots = [ 157258 , 129516 , 140114 , 131341 ];
RV0_roots = [ 198857 , 209309 , 186671 ];



figure()
plotMESH( H0 , 'ne' );
hplot3d( H0.xyz( LV0_roots ,:) , 'o1kb20' );
hplot3d( H0.xyz( RV0_roots ,:) , 'o1kr20' );
headlight


%%

H1 = read_VTK(  'D:\ARVC meshing automatic\patients\patient06\HEART.vtk'); %our target heart

roots=load('D:\ARVC meshing automatic\patients\patient06\start.txt');





figure()
plotMESH( H1 , 'ne' );
hplot3d( H1.xyz( vertcat(roots(1:4), roots(6:7)) ,:) , 'o1kb20' );
headlight

pos_new=[29.52 -59.82 -19.88];
distance=@(v1,v2) sqrt((v1(:,1)-v2(:,1)).^2+(v1(:,2)-v2(:,2)).^2+(v1(:,3)-v2(:,3)).^2);
dist=distance(pos_new,H1.xyz);
[~,id]=min(dist);

pos_new2=[40.53 -83.21 21.36];
dist=distance(pos_new2,H1.xyz);
[~,id2]=min(dist);

roots_new= vertcat(roots(1:5),roots(7) );

dlmwrite('D:\ARVC meshing automatic\patients\patient06\eikonal06_fine_rootNodesRootRVSeptUpAnteriorForward.csv',roots_new,'delimiter',',','precision',10);


figure()
%plotMESH( H0 ,'ne');
hplotMESH( H1 , 'ne' ,'FaceColor',[1 1 1]*0.6)
hplot3d( H1.xyz( vertcat(roots_new+1),:) , 'o1kb20' );
headlight; axis(objbounds);

%HD_loaded = load('Z:\Dropbox\patients\HD.mat');

%%
H1_full = read_CHASTE(  'D:\ARVC meshing automatic\patients\patient06\chaste06\HEART'); %our target heart

u = loadv( strcat(SUBJECT_DIR, 'mpp\Hmapping') , 'u' ); %transfer matrix u
%this is a function converting coordinates from H0 to H1
%it is generated by running the mpp_Heart_Mapping Function of Ernesto s
%pipeline!!

if 0
% idx=vtkClosestPoint(Mesh(EPI1),EPI0.xyz(1,:));
% vec = EPI1.xyz(idx,:) - EPI0.xyz(1,:);
% distance = sqrt(vec * vec');
% 
% column = ones(size(H0.xyz,1),1);
% vec_repeated(:,1) = column.*vec(1);
% vec_repeated(:,2) = column.*vec(2);
% vec_repeated(:,3) = column.*vec(3);
% 
% if distance > 100
%     H0.xyz = H0.xyz+vec_repeated;
% end

HD = H0;
HD.xyz = u( HD.xyz );
figure()
plotMESH( HD ,'ne');
hplotMESH( H1 , 'ne' ,'FaceColor',[1 1 1]*0.6)
headlight; axis(objbounds);
axis off
end


%%

LV1_roots = vtkClosestPoint( Mesh(H1) , u( H0.xyz( LV0_roots ,:) ) );
RV1_roots = vtkClosestPoint( Mesh(H1) , u( H0.xyz( RV0_roots ,:) ) );

try 
    root_nodes =load(strcat(SUBJECT_DIR,'\start.txt'));
    LV1_roots=[];
    RV1_roots =[];
    LV1_roots =root_nodes(1:4)+1;
    RV1_roots = root_nodes(5:end)+1;
catch
    printf('no root node file found');
end

figure()
plotMESH( H1 , 'ne' );
hplot3d( H1.xyz( LV1_roots(1) ,:) , 'o1kb20' );
hplot3d( H1.xyz( LV1_roots(2) ,:) , 'o1kg20' );
hplot3d( H1.xyz( LV1_roots(3) ,:) , 'o1km20' );
hplot3d( H1.xyz( LV1_roots(4) ,:) , 'o1ky20' );
headlight
hplot3d( H1.xyz( RV1_roots(1) ,:) , 'o1kr20' );
hplot3d( H1.xyz( RV1_roots(2) ,:) , 'o1kb20' );
hplot3d( H1.xyz( RV1_roots(3) ,:) , 'o1kg20' );

dlmwrite(strcat(SUBJECT_DIR,'start_pos.txt'),[vertcat(H1.xyz(LV1_roots,:),H1.xyz(RV1_roots,:))]/10.0)




LV1_roots_chaste  = LV1_roots-1; %convert to chaste indexing
RV1_roots_chaste = RV1_roots-1;

fid = fopen(strcat(SUBJECT_DIR,'start.txt'),'w');
fprintf(fid,'%u \n',[vertcat(LV1_roots_chaste,RV1_roots_chaste)]')
fclose(fid);



%now implement a technique to move root nodes around 
%first transpose the ventricles to align the apexbase with the z axis
apexbase = load(strcat(TORSO_MODEL_DIR,'apexbase'));
apexbase_vec = apexbase(2,:)-apexbase(1,:);
apexbase_dist = sqrt(apexbase_vec * apexbase_vec')*10;
z_vec = [0,0,1];
subject= SUBJECT_DIR;
iZ = loadv( fullfile(strcat(subject,'mpp\'),'HM'),'iZ');
H1iz=transform(H1,iZ);
figure()
plotMESH( H1iz , 'ne' );
headlight

LV=H1_full.lv;

lies_on_rv_endo = zeros(size(H1iz.xyz,1),1);
lims=range(H1_full.rv);
limsLV=range(H1_full.lv);

lies_on_rv_endo(lims(1):lims(2)) =1;
lies_on_lv_endo = zeros(size(H1iz.xyz,1),1);
lies_on_lv_endo(limsLV(1):limsLV(2))=1;


%now select a rootnode
%1st root is the one on the LV anterior base
%move it a set distance, say 48.5% of the apexbase vector magnitude
distances = [repmat(H1iz.xyz(LV1_roots(1),: ),size(H1iz.xyz,1),1)]-H1iz.xyz(:,:);
distances = sqrt(distances(:,1).*distances(:,1)+distances(:,2).*distances(:,2)+distances(:,3).*distances(:,3));

%criteria for selection
% 1) root node must be at the same z-height as original node
% 2) it must be moving around the LV endocardium along the lateral wall,
% hence the distance between it and a RV rootnode 3 must be increasing

dist = 0.48*apexbase_dist<distances &  distances < .49*apexbase_dist;
dist_ind = find(dist);

z_height = H1iz.xyz(LV1_roots(1),3)*0.98> H1iz.xyz(:,3) &  H1iz.xyz(:,3)> H1iz.xyz(LV1_roots(1),3)*1.02;

%the new root node has to be closer to from LV node 2

dist_node2 = [repmat(H1iz.xyz(LV1_roots(2),: ),size(H1iz.xyz,1),1)]-H1iz.xyz(:,:);
dist_node2 = sqrt(dist_node2(:,1).*dist_node2(:,1)+dist_node2(:,2).*dist_node2(:,2)+dist_node2(:,3).*dist_node2(:,3));

prev_dist_LV1_LV2_nodes = H1iz.xyz(LV1_roots(1))-H1iz.xyz(LV1_roots(2));
prev_dist_LV1_LV2_nodes = sqrt(prev_dist_LV1_LV2_nodes*prev_dist_LV1_LV2_nodes);

further_away = dist_node2<dist_node2(LV1_roots(1));

score = dist+z_height+further_away;
candidates = find(score ==3);

%pick the index closest to node 2
[~,id]= min(dist_node2(candidates));

LV1_roots_new(1) = candidates(id);

%now try to move root node in basal RV free wall (RV root node 2) radially towards the
%posterior base

%conditions
%1)Distances between the new root node and old root node need to be 0.2 x
%apexbase dist
%2) The z -height of the candidate root node must be the same as the
%original
%3) The new root node must be closer to the lv node 3 than the previous
%root node
%) it must lie on the endocardium
%move it a set distance, say 48.5% of the apexbase vector magnitude
distances = [repmat(H1iz.xyz(RV1_roots(2),: ),size(H1iz.xyz,1),1)]-H1iz.xyz(:,:);
distances = sqrt(distances(:,1).*distances(:,1)+distances(:,2).*distances(:,2)+distances(:,3).*distances(:,3));
dist = 0.45*apexbase_dist<distances &  distances < .46*apexbase_dist;

z_height = H1iz.xyz(RV1_roots(2),3)*0.98> H1iz.xyz(:,3) &  H1iz.xyz(:,3)> H1iz.xyz(RV1_roots(2),3)*1.02;
const_height=find(z_height ==1);
% 
% figure()
% plotMESH( H1iz , 'ne' );
% hplot3d(H1iz.xyz(const_height,:),'o1kb20');
% headlight
% hplot3d( H1iz.xyz( RV1_roots(2) ,:) , 'o1km20' );

%the new root node has to be closer to LV node3

dist_node2 = [repmat(H1iz.xyz(LV1_roots(3),: ),size(H1iz.xyz,1),1)]-H1iz.xyz(:,:);
dist_node2 = sqrt(dist_node2(:,1).*dist_node2(:,1)+dist_node2(:,2).*dist_node2(:,2)+dist_node2(:,3).*dist_node2(:,3));

prev_dist_RV1_RV2_nodes = H1iz.xyz(LV1_roots(3))-H1iz.xyz(RV1_roots(2));
prev_dist_RV1_RV2_nodes = sqrt(prev_dist_RV1_RV2_nodes*prev_dist_RV1_RV2_nodes);

further_away = dist_node2<dist_node2(RV1_roots(2));

score = dist+z_height+further_away+lies_on_rv_endo;
candidates = find(score ==4);

%pick the index closest to node 2
[~,id]= min(dist_node2(candidates));
RV1_roots_new(2) = candidates(id);


%now try to move root node in basal RV anterior  (RV root node 1) vertically towards the
%apex
%conditions

%1) The z -height of the candidate root node must be 1.485 x z-height of
%original RV root node 1
%2) The root node must have a similar xy value to the original
%3) new node needs to be lie on the endocardial RV
z_height = H1iz.xyz(RV1_roots(1),3)-0.3*apexbase_dist> H1iz.xyz(:,3) &  H1iz.xyz(:,3)> H1iz.xyz(RV1_roots(1),3)-0.31*apexbase_dist;




y_change = abs(H1iz.xyz(:,2)-repmat(H1iz.xyz(RV1_roots(1),2),size(H1iz.xyz,1),1));
x_change = abs(H1iz.xyz(:,1)-repmat(H1iz.xyz(RV1_roots(1),1),size(H1iz.xyz,1),1));
xy_change = x_change < 30 &y_change <30; 

score = z_height+xy_change+lies_on_rv_endo;
candidates = find(score ==3);

%pick the index closest to upper z boundary
[~,id]= max(H1iz.xyz(candidates,3));
RV1_roots_new(1) = candidates(id);

%now lets move node RV 3

%conditions:
%1) its going to move vertically upwards by 41%of the apexbase distance
%2) the x and y coords are similar to the previous root node 3
%3) the new root node is on the endocardium

z_height = H1iz.xyz(RV1_roots(3),3)+0.41*apexbase_dist> H1iz.xyz(:,3) &  H1iz.xyz(:,3)> H1iz.xyz(RV1_roots(3),3)+apexbase_dist*0.40;
if 0
figure()
plotMESH( H1 , 'ne' );
hplot3d( H1.xyz( find(z_height==1) ,:) , 'o1kb20' );
end
y_change = abs(H1iz.xyz(:,2)-repmat(H1iz.xyz(RV1_roots(3),2),size(H1iz.xyz,1),1));
x_change = abs(H1iz.xyz(:,1)-repmat(H1iz.xyz(RV1_roots(3),1),size(H1iz.xyz,1),1));

total_change = x_change+y_change;

xy_change = x_change < 5 &y_change <5; 

score = z_height+xy_change+lies_on_rv_endo;
candidates = find(score ==3);

%pick the index with largest z-coord (closest to the z=0 line)
[~,id]= min(total_change(candidates));
RV1_roots_new(3) = candidates(id);

%lets move root node LV3 upwards
%conditions:
%1) its going to move vertically upwards by 30%of the apexbase distance
%2) the x and y coords are similar to the previous root node 3
%3) the new root node is on the endocardium

z_height = H1iz.xyz(LV1_roots(3),3)+0.31*apexbase_dist> H1iz.xyz(:,3) &  H1iz.xyz(:,3)> H1iz.xyz(LV1_roots(3),3)+apexbase_dist*0.30;

y_change = abs(H1iz.xyz(:,2)-repmat(H1iz.xyz(LV1_roots(3),2),size(H1iz.xyz,1),1));
x_change = abs(H1iz.xyz(:,1)-repmat(H1iz.xyz(LV1_roots(3),1),size(H1iz.xyz,1),1));


xy_change = x_change < 20 &y_change <20; 

score = z_height+xy_change+lies_on_lv_endo;
candidates = find(score ==3);

[~,id]=min(total_change(candidates));
%pick the index with largest z-coord (closest to the z=0 line)
LV1_roots_new(3) = candidates(id);

%now lets move the septal LV root node 4 up
%conditions:
%1) its going to move vertically upwards by 40%of the apexbase distance
%2) the x and y coords are similar to the previous root node 3
%3) the new root node is on the endocardium

z_height = H1iz.xyz(LV1_roots(4),3)+0.41*apexbase_dist> H1iz.xyz(:,3) &  H1iz.xyz(:,3)> H1iz.xyz(LV1_roots(4),3)+apexbase_dist*0.40;
y_change = abs(H1iz.xyz(:,2)-repmat(H1iz.xyz(LV1_roots(4),2),size(H1iz.xyz,1),1));
x_change = abs(H1iz.xyz(:,1)-repmat(H1iz.xyz(LV1_roots(4),1),size(H1iz.xyz,1),1));

total_change = x_change +y_change;
xy_change = x_change < 20 &y_change <20; 
score = z_height+xy_change+lies_on_lv_endo;
candidates = find(score ==3);

[~,id]=min(total_change(candidates));
%pick the index with largest z-coord (closest to the z=0 line)
LV1_roots_new(4) = candidates(id);

%now lets move the apical LV root node 2 radially towards the anterior
%conditions:
%1) its going to move radially towards the anterior apex by 10%of the apexbase distance
    % ie the distance between the old node and the new node is going to be
    % 10% of apexbase dist
%2) the z coords are similar to the previous root node 2
%3) the new root node is on the endocardium
%4) pick root node so that its closest to original LV1 root node

distances = [repmat(H1iz.xyz(LV1_roots(2),: ),size(H1iz.xyz,1),1)]-H1iz.xyz(:,:);
distances = sqrt(distances(:,1).*distances(:,1)+distances(:,2).*distances(:,2)+distances(:,3).*distances(:,3));
dist = 0.3*apexbase_dist<distances &  distances < .31*apexbase_dist;

z_height = H1iz.xyz(LV1_roots(2),3)*0.98> H1iz.xyz(:,3) &  H1iz.xyz(:,3)> H1iz.xyz(LV1_roots(2),3)*1.02;

score = z_height+dist+lies_on_lv_endo;
candidates = find(score ==3);

distances = [repmat(H1iz.xyz(LV1_roots(1),: ),size(candidates,1),1)]-H1iz.xyz(candidates,:);
distances = sqrt(distances(:,1).*distances(:,1)+distances(:,2).*distances(:,2)+distances(:,3).*distances(:,3));

[~,id] = min(distances);

LV1_roots_new(2) = candidates(id);


fid = fopen(strcat(SUBJECT_DIR,'permuted.txt'),'w');
fprintf(fid,'%5.0f \n',[LV1_roots_new-1 RV1_roots_new-1]);%chaste ordering starts from 0
fclose(fid)

fid = fopen(strcat(SUBJECT_DIR,'permuted_nodes_lv.txt'),'w');
fprintf(fid,'%5.0f \n',[LV1_roots_new-1 ]);%chaste ordering starts from 0
fclose(fid)

fid = fopen(strcat(SUBJECT_DIR,'permuted_nodes_rv.txt'),'w');
fprintf(fid,'%5.0f \n',[LV1_roots_new-1 ]);%chaste ordering starts from 0
fclose(fid)


figure()
plotMESH( H1 , 'ne' );
hplot3d( H1.xyz( LV1_roots(1) ,:) , 'o1kb20' );
hplot3d( H1.xyz( LV1_roots(2) ,:) , 'o1kg20' );
hplot3d( H1.xyz( LV1_roots(3) ,:) , 'o1km20' );
hplot3d( H1.xyz( LV1_roots(4) ,:) , 'o1ky20' );
hplot3d( H1.xyz( LV1_roots_new(1) ,:) , 'o1kb20' );
hplot3d( H1.xyz( LV1_roots_new(2) ,:) , 'o1kg20' );
hplot3d( H1.xyz( LV1_roots_new(3) ,:) , 'o1km20' );
hplot3d( H1.xyz( LV1_roots_new(4) ,:) , 'o1ky20' );

headlight
hplot3d( H1.xyz( RV1_roots(1) ,:) , 'o1kr20' );
hplot3d( H1.xyz( RV1_roots(2) ,:) , 'o1kb20' );
hplot3d( H1.xyz( RV1_roots(3) ,:) , 'o1kg20' );
hplot3d( H1.xyz( RV1_roots_new(2) ,:) , 'o1kb20' );
hplot3d( H1.xyz( RV1_roots_new(1) ,:) , 'o1kr20' );
hplot3d( H1.xyz( RV1_roots_new(3) ,:) , 'o1kg20' );


figure()
plotMESH( H1 , 'ne' );
hplot3d( H1.xyz( LV1_roots ,:) , 'o1kb20' );
hplot3d( H1.xyz( RV1_roots ,:) , 'o1kr20' );
headlight;