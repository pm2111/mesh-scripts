addpath(genpath('C:/Users/petnov/Dropbox/'));
addpath(genpath('C:/Users/petnov/Dropbox/qrs/'));
addpath(genpath('C:/Users/petnov/Dropbox/sim 3d scripts/'));
addpath(genpath('C:/Users/petnov/Dropbox/shared - copy/mesh scripts/'));

addpath('C:\Users\petnov\Dropbox\sim3d scripts\');
addpath(genpath('C:\Users\petnov\Dropbox\mesh scripts\'));
addpath(genpath('C:\Users\petnov\Dropbox\qrs\'));
addpath(genpath('C:\Users\petnov\Dropbox\tools\'));

addpath('C:\Users\petnov\Dropbox\shared - copy\IO\');
 patient_sim_folder =  strcat('D:\ARVC meshing automatic\patients\DTI003\');

    %H0 = read_VTK( fullfile( TORSO_MODEL_DIR , 'HEART_with_UPS_louie.vtk' ) ); %this is louie's heart
H_mesh=struct;
H_mesh.xyz=csvread('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_xyz.csv');
H_mesh.tri=csvread('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_tri.csv');
H_mesh.rv=csvread('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_rvface.csv');
H_mesh.lv=csvread('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_lvface.csv');
H_mesh.epi=csvread('D:\ARVC meshing automatic\patients\DTI003\DTI003_mesh\DTI003_coarse_epiface.csv');
H_mesh.face=vertcat(H_mesh.epi,H_mesh.lv,H_mesh.rv);

dists = dlmread('/data/meshes/nejib/dists.txt',' ');

H_mesh.xyzDistLV=dists(:,2);
H_mesh.xyzDistRV=dists(:,3);
H_mesh.xyzDistEpi=dists(:,1);

distlv=dists(:,2);
distrv=dists(:,3);
distepi=dists(:,1);


mesh=struct;
mesh.xyz=H_mesh.xyz;
% mesh.xyzDistLV=dists(:,2);
% mesh.xyzDistRV=dists(:,3);
% mesh.xyzDistEpi=dists(:,1);
mesh.tri=H_mesh.tri(H_mesh.triATTS==0,:);
mesh.triFibers=H_mesh.triORTHO(H_mesh.triATTS==0,1:3);
mesh.triX=H_mesh.triORTHO(H_mesh.triATTS==0,1);
mesh.triY=H_mesh.triORTHO(H_mesh.triATTS==0,2);
mesh.triZ=H_mesh.triORTHO(H_mesh.triATTS==0,3);

write_VTK(mesh,'/data/meshes/nejib/mesh.vtk','binary');

edges=dlmread('/data/meshes/patient06/chaste06/eikonal06/HEART_fine_edges.csv');

%analyse the distances chaste uses to determine what area of the heart each
%node belongs to

for i=1:size(dists,1)
    if distrv(i) >= distepi(i) && distrv(i)*2>= distlv(i)
        region_old(i)=1;
    
    else % distlv(i) >= distepi(i) && distlv(i) >=distrv(i)
                region_old(i)=2;

    end
    if distepi(i)*2 >= distlv(i) && distepi(i) >=distrv(i)

        if distlv(i) < 2/3*(distlv(i)+distrv(i))
            region_old(i)=3;
        else
            region_old(i)=4;
        end
    end

end


figure()
hold on
scatter3(mesh.xyz(idx,1),mesh.xyz(idx,2),mesh.xyz(idx,3),20,region_old(idx),'filled')
%quiver3(mesh.xyz(idx,1),mesh.xyz(idx,2),mesh.xyz(idx,3),fibres(:,1),fibres(:,2),fibres(:,3),0.5)
colorbar()


% region=region_old;
% sz= max(distrv(region_old==4)) ;
% mn=mean(distlv(region_old==4));
% st=std(distlv(region_old==4));
% for i=1:size(dists,1)
%     if region_old(i)==2 && distrv(i)>= distlv(i)*0.1
%         region(i)=1;
%     end
% 
% end

figure()
hold on
scatter3(mesh.xyz(idx,1),mesh.xyz(idx,2),mesh.xyz(idx,3),20,region(idx),'filled')
%quiver3(mesh.xyz(idx,1),mesh.xyz(idx,2),mesh.xyz(idx,3),fibres(:,1),fibres(:,2),fibres(:,3),0.5)
colorbar()


dist_epi_combined=distlv(H_mesh.epi)+distrv(H_mesh.epi);
[~,idx_insertion]=min(dist_epi_combined);

vec_insertion_epi = mesh.xyz(H_mesh.epi,:)-mesh.xyz(H_mesh.epi(idx_insertion),:);
dist_insertion_epi=norm(vec_insertion_epi(:,1),vec_insertion_epi(:,2),vec_insertion_epi(:,3));

idx_far_insertion=find(dist_insertion_epi>5);

[~,idx3]=min(dist_epi_combined(idx_far_insertion));

idx_insertion_epi2=idx_far_insertion(idx3);

figure()
hold on
scatter3(mesh.xyz(idx,1),mesh.xyz(idx,2),mesh.xyz(idx,3),20,region(idx),'filled')
%quiver3(mesh.xyz(idx,1),mesh.xyz(idx,2),mesh.xyz(idx,3),fibres(:,1),fibres(:,2),fibres(:,3),0.5)
colorbar()
scatter3(H_mesh.xyz(H_mesh.epi(idx_insertion),1),H_mesh.xyz(H_mesh.epi(idx_insertion),2),6,30, 'r','filled')
scatter3(H_mesh.xyz(H_mesh.epi(idx_insertion_epi2),1),H_mesh.xyz(H_mesh.epi(idx_insertion_epi2),2),6,30, 'r','filled')


H_mesh.triORTHO=H_mesh.triORTHO(:,:);
mesh=struct;
mesh.tri=H_mesh.tri;
mesh.xyz=H_mesh.xyz;
mesh.triFIBERS=H_mesh.triORTHO(:,:);
mesh.celltype=10;
mesh.xyzregions=region_old';
write_VTK(mesh,'/data/meshes/chaste03/mesh_new_sep.vtk','binary');

H_mesh.celltype =10;

epicardium = H_mesh.epi;
rv= H_mesh.rv;
lv= H_mesh.lv;

%coords = Hsmall.xyz(Hsmall.tri(:,1),:);


H_mesh.xyzORTHO= H_mesh.triORTHO;
mesh=struct;
mesh.xyz =H_mesh.xyz;
mesh.triORTHO=H_mesh.triORTHO(:,:,1);
mesh.tri=H_mesh.tri;
mesh.Celltype=10;


plane = @(x,y,z) -5.5+z;
norm = @(x,y,z) sqrt(x.^2+y.^2+z.^2);
d=3;
dist = plane(mesh.xyz(:,1),mesh.xyz(:,2),mesh.xyz(:,3));

idx = find(abs(dist) <0.0001);
id_tri=[];
for i=1:size(idx,1)
    dummy= find(mesh.tri(:,1)==idx(i));
    try
        id_tri(i)=dummy(1,1);
    catch
                dummy= find(mesh.tri(:,2)==idx(i));
                if size(dummy,1)==0
                     dummy= find(mesh.tri(:,3)==idx(i));

                end
                  if size(dummy,1)==0
                     dummy= find(mesh.tri(:,4)==idx(i));

                  end
            if size(dummy,1)>0
                id_tri(i)=dummy(1,1);
            end
    end

end



fibres=H_mesh.triORTHO(id_tri',1:3,1);

in_plane = [mesh.xyz(idx(1),1)-mesh.xyz(idx(5),1) mesh.xyz(idx(1),2)-mesh.xyz(idx(5),2) mesh.xyz(idx(1),3)-mesh.xyz(idx(5),3)];

apexbase= load('/data/meshes/nejib/apexbase');
apexbase_vec = -(apexbase(1,:)-apexbase(2,:));
angle =[];
for i=1:size(fibres,1)
    angle(i) = 180/pi*acos(dot(apexbase_vec,fibres(i,:))/(norm(apexbase_vec(1),apexbase_vec(2),apexbase_vec(3))*norm(fibres(i,1),fibres(i,2),fibres(i,3))));


end
 


figure()
hold on
scatter3(mesh.xyz(idx,1),mesh.xyz(idx,2),mesh.xyz(idx,3),20,angle(idx),'filled')
%quiver3(mesh.xyz(idx,1),mesh.xyz(idx,2),mesh.xyz(idx,3),fibres(:,1),fibres(:,2),fibres(:,3),0.5)
colorbar()


figure()
scatter3(mesh.xyz(idx,1),mesh.xyz(idx,2),mesh.xyz(idx,3),20,angle_in_plane','filled')

%streamline(mesh.xyz(idx,1),mesh.xyz(idx,2),mesh.xyz(idx,3),mesh.triORTHO(idx,1),mesh.triORTHO(idx,2),mesh.triORTHO(idx,3))
colorbar
hold on
title("Angle of fibre diretion with respect to apex-base vector")
grid off
scatter3(mesh.xyz(idx(1),1),mesh.xyz(idx(1),2),mesh.xyz(idx(1),3),50,'r','filled')
scatter3(mesh.xyz(idx(5),1),mesh.xyz(idx(5),2),mesh.xyz(idx(5),3),50,'g','filled')

coords = H_mesh.xyz(H_mesh.tri(:,1),:);
%coords=H_mesh.xyz;

lv_coords= H_mesh.xyz(lv+1,:);
rv_coords= H_mesh.xyz(rv+1,:);
epi_coords= H_mesh.xyz(epicardium+1,:);

xs=coords(1:2,:);
indices =[];
for i=1:3000
    rnd=randi([1 max(rv)]);
    x=coords(rnd,:);
     while (isempty(find(lv_coords(:,1)==x(1)))==1 && isempty(find(lv_coords(:,2)==x(2)))==1 &&isempty(find(lv_coords(:,3)==x(3)))==1)% || (isempty(find(rv_coords(:,1)==x(1)))==1 && isempty(find(rv_coords(:,2)==x(2)))==1 &&isempty(find(rv_coords(:,3)==x(3)))==1) || (isempty(find(epi_coords(:,1)==x(1)))==1 && isempty(find(epi_coords(:,2)==x(2)))==1 &&isempty(find(epi_coords(:,3)==x(3)))==1)
         ind2=randi([1,size(coords,1)]);
         x=coords(ind2,:);
     end
    dummy= ones(size(xs,1),3);
    dummy(:,1)=dummy(:,1)*x(1);
    dummy(:,2)=dummy(:,1)*x(2);
    dummy(:,3)=dummy(:,1)*x(3);

    
    dist = ((dummy(:,1)-xs(:,1)).^2 + (dummy(:,2)-xs(:,2)).^2+(dummy(:,3)-xs(:,3)).^2).^0.5 ;

  %  if isempty(find(dist<30))==1
        indices =[indices,ind2];
        xs(end+1,1)=coords(ind2,1);
        xs(end,2)=coords(ind2,2);
        xs(end,3)=coords(ind2,3);
   % end

end

for i=1:3000
    rnd=randi([1 max(rv)]);
    x=coords(rnd,:);
     while  (isempty(find(rv_coords(:,1)==x(1)))==1 && isempty(find(rv_coords(:,2)==x(2)))==1 &&isempty(find(rv_coords(:,3)==x(3)))==1) %|| (isempty(find(epi_coords(:,1)==x(1)))==1 && isempty(find(epi_coords(:,2)==x(2)))==1 &&isempty(find(epi_coords(:,3)==x(3)))==1)
        ind2=randi([1,size(coords,1)]);
         x=coords(ind2,:);
     end
    dummy= ones(size(xs,1),3);
    dummy(:,1)=dummy(:,1)*x(1);
    dummy(:,2)=dummy(:,1)*x(2);
    dummy(:,3)=dummy(:,1)*x(3);

    
    dist = ((dummy(:,1)-xs(:,1)).^2 + (dummy(:,2)-xs(:,2)).^2+(dummy(:,3)-xs(:,3)).^2).^0.5 ;

   % if isempty(find(dist<30))==1
        indices =[indices,ind2];
        xs(end+1,1)=coords(ind2,1);
        xs(end,2)=coords(ind2,2);
        xs(end,3)=coords(ind2,3);
  %  end

end

for i=1:3000
    rnd=randi([1 max(rv)]);
    x=coords(rnd,:);
     while  (isempty(find(epi_coords(:,1)==x(1)))==1 && isempty(find(epi_coords(:,2)==x(2)))==1 &&isempty(find(epi_coords(:,3)==x(3)))==1)
        ind2=randi([1,size(coords,1)]);
         x=coords(ind2,:);
     end
    dummy= ones(size(xs,1),3);
    dummy(:,1)=dummy(:,1)*x(1);
    dummy(:,2)=dummy(:,1)*x(2);
    dummy(:,3)=dummy(:,1)*x(3);

    
    dist = ((dummy(:,1)-xs(:,1)).^2 + (dummy(:,2)-xs(:,2)).^2+(dummy(:,3)-xs(:,3)).^2).^0.5 ;

   % if isempty(find(dist<30))==1
        indices =[indices,ind2];
        xs(end+1,1)=coords(ind2,1);
        xs(end,2)=coords(ind2,2);
        xs(end,3)=coords(ind2,3);
    %end

end

x = H_mesh.xyz(,1); 
y = H_mesh.xyz(:,2); 
z = H_mesh.xyz(:,3); 

vx = 30 + 30*cos((Y-100)*pi/100); 
vy = 5*cos(X/10).*cos((Y-100)*pi/100); 
figure
pcolor(X,Y,hypot(vx,vy))
shading interp
N = 75; 
xstart = max(x)*rand(N,1); 
ystart = max(y)*rand(N,1); 
h=streamline(X,Y,vx,vy,xstart,ystart);
set(h,'color','red')

H_mesh.xyzFIBERS=fibers(inds,1:3);

write_VTK(H_mesh,'D:\ARVC meshing automatic\patients\patient05\fibers_small.vtk','binary')

y = H_mesh.xyz(:,2); 
clearvars
    inds = [];
    for i=1:size(H_mesh.xyz,1)
        [row,col]=find(H_mesh.tri==i);
         inds=[inds,row(1)];
    end
        
    mag_V = sqrt(fibers(:,1).^2 + fibers(:,2).^2 + fibers(:,3).^2);
    % Filtering out zero velocity walls.
    index = mag_V > 0.1;
    xx = H_mesh.xyz(:,1);
    yy = H_mesh.xyz(:,2);
    zz = H_mesh.xyz(:,3);
    Vxx = fibers(inds,1);
    Vyy = fibers(inds,2);
    Vzz = fibers(inds,3);
    FVx = scatteredInterpolant(xx,yy,zz,Vxx,'linear','none');
    FVy = scatteredInterpolant(xx,yy,zz,Vyy,'linear','none');
    FVz = scatteredInterpolant(xx,yy,zz,Vzz,'linear','none');
    % This is a very ugly grid for your problem but it works.
    % A cylindrical grid would probably be more resource efficient
    elements = 40; % how many subdivisions per dimensions.
    % Mind you, this is 3D. Total memory gets large fast
    [X3, Y3, Z3] = meshgrid(linspace(min(xx),max(xx),elements),...
        linspace(min(yy),max(yy),elements),...
        linspace(min(zz),max(zz),elements));
    V3x = FVx(X3,Y3,Z3);
    V3y = FVy(X3,Y3,Z3);
    V3z = FVz(X3,Y3,Z3);
    mag_V3 = sqrt(V3x.^2 + V3y.^2 + V3z.^2); % velocity vector magnitude
    starting_x = min(xx)*ones(5,5); % starting points for streamlines
    [starting_y, starting_z] = meshgrid(linspace(min(yy),max(yy),5), ...
        linspace(min(zz),max(zz),5));
    figure()
    clf
   plot3(xx,yy,zz,'.k','MarkerSize',4);
    hold on
    %quiver3(xx,yy,zz, V3x,V3y,V3z);
   % h = slice(X3,Y3,Z3,mag_V3,[],0,-0.02:.02:0);
   % set(h,'EdgeColor','w','FaceAlpha',0.75');
    h = streamline(X3,Y3,Z3,V3x,V3y,V3z,starting_x,starting_y,starting_z);
    set(h,'Color','r','LineWidth',4);
    axis equal
    grid on
    hold off

figure()


hold on
hplotMESH(Hsmall,'ne')
streamline(fibers(:,1:3),fibers(:,4:6),2,'LineWidth',4)
headlight
%q=quiver3(coords(1:100000:end,1),coords(1:100000:end,2),coords(1:100000:end,3),fibers(1:100000:end,4),fibers(1:100000:end,5),fibers(1:100000:end,6),2,'LineWidth',2.0)
%q.Color='r';
alpha(0.6)
axis off