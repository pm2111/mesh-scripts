addpath(genpath('C:\Users\petnov\Dropbox\shared - copy\'));
addpath(genpath('C:\Users\petnov\Dropbox\mesh scripts'))
addpath(genpath('C:\Users\peter\Dropbox\mesh scripts'))
addpath(genpath('C:\Users\peter\Dropbox\shared - copy'))

z =@(x,y,a,b,c)  c*(sqrt((1-x.^2/a.^2 -y.^2/b.^2)));

[X,Y]=meshgrid(-30:0.5:30,-30:0.5:30)



Z = z(X,Y,-15,-15,-15);

Z2 = z(X,Y,-17,-17,-17);


figure()
surf(X,Y,real(Z))


figure()
hold on
surf(X,Y,real(Z2))
surf(X,Y,real(Z))
alpha =0.5

keep = Z<-10;
keep2=Z2<-10;
figure()
scatter3(X(keep),Y(keep),Z(keep))

%now make a mesh by joining all nodes to nearest 2 nodes

nodes_inner = [X(keep),Y(keep),Z(keep)];
nodes_outer = [X(keep2),Y(keep2),Z2(keep2)];

[t]=MyCrustOpen(nodes_inner);
[t2]=MyCrustOpen(nodes_outer);
mesh.tri=t;
mesh.xyz=nodes_inner;
mesh2.tri=t2;
mesh2.xyz=nodes_outer;
figure()
hplotMESH(mesh);
hplotMESH(mesh2);
axis equal

idx=find(abs(nodes_inner(:,3)+10) <0.35);

figure()
hold on
hplotMESH(mesh);
scatter3(mesh.xyz(idx,1),mesh.xyz(idx,2),mesh.xyz(idx,3),'filled')

idx2=find(abs(nodes_outer(:,3)+10) <0.35);

figure()
hold on
hplotMESH(mesh2);
scatter3(mesh2.xyz(idx2,1),mesh2.xyz(idx2,2),mesh2.xyz(idx2,3),'filled')

mesh3.xyz=vertcat(mesh.xyz(idx,:),mesh2.xyz(idx2,:));
mesh3.tri= MyCrustOpen(mesh3.xyz);

figure()
hplotMESH(mesh3)
mesh_combined=struct;
mesh_combined.xyz=vertcat(mesh.xyz,mesh2.xyz);
mesh_combined.tri=vertcat(mesh.tri,mesh2.tri+size(mesh.xyz,1));
figure()
hplotMESH(mesh_combined);

idx = find(abs(mesh_combined.xyz(:,3)+10) <0.35);
mesh_combined.tri=vertcat(mesh_combined.tri,idx(t4));
figure()
hplotMESH(mesh_combined);
EDGE_LENGTH=0.4;
  [~,M] = jigsaw_surface_2_volume( mesh_combined , 'delfront' , 'absolute','geom_feat',true,'hfun_hmax',EDGE_LENGTH,'hfun_hmin',EDGE_LENGTH*0.9);


  figure()
  hplotMESH(M);
  
Htetgen = gmshfn(M,['v']);
figure()
hplotMESH(Htetgen)

write_VTK(Htetgen,'D:\meshes\ellipsoid.vtk','binary')



[Cnodes,Celems,fib]=benchmark_ellipse_linear();

figure
plot3(Cnodes(:,1),Cnodes(:,2),Cnodes(:,3),'+')
grid on
axis equal

ellipsoid=struct;
ellipsoid.xyz=Cnodes;
ellipsoid.tri=Celems;

%centred at (0,0,z)
ellipsoid.epi=[];
ellipsoid.endo=[];
un=unique(Cnodes(:,3));

rs = @(x,y) sqrt(x.^2+y.^2);
radii = rs(Cnodes(:,1),Cnodes(:,2));

for i=1:size(un,1)
    idx = Cnodes(:,3)==un(i);
    radii_unique=unique(radii(idx));
    a= min(radii_unique);
    b=max(radii_unique);
    if i<4
        idx=find(Cnodes(:,3)==un(i));
        ellipsoid.epi=vertcat(ellipsoid.epi,idx);
    
    elseif 3<i && i<6
         idx=find(Cnodes(:,3)==un(i));   
         ellipsoid.endo=vertcat(ellipsoid.endo,idx);
    
    elseif abs(a-b) <0.1 && b > max(radii(ellipsoid.epi))
                 idx=find(Cnodes(:,3)==un(i));   

                ellipsoid.epi=vertcat(ellipsoid.epi,idx);
    elseif i==size(un,1)
        idx =find(Cnodes(:,3)==un(i) & radii<a+0.5);
                         ellipsoid.endo=vertcat(ellipsoid.endo,idx);
                         
                idx =find(Cnodes(:,3)==un(i) & radii>a+0.5);
        ellipsoid.epi=vertcat(ellipsoid.epi,idx);

    else
                         idx=find(Cnodes(:,3)==un(i));   

                 ellipsoid.endo=vertcat(ellipsoid.endo,idx);
    end

end

figure()
hold on
scatter3(ellipsoid.xyz(ellipsoid.epi,1),ellipsoid.xyz(ellipsoid.epi,2),ellipsoid.xyz(ellipsoid.epi,3),'r','filled')
scatter3(ellipsoid.xyz(ellipsoid.endo,1),ellipsoid.xyz(ellipsoid.endo,2),ellipsoid.xyz(ellipsoid.endo,3),'b','filled')


mesh_epi= struct;
mesh_epi.xyz=Cnodes(ellipsoid.epi,:);
mesh_epi.tri=MyCrustOpen(mesh_epi.xyz);

mesh_endo=struct;
mesh_endo.xyz=Cnodes(ellipsoid.endo,:);
mesh_endo.tri=MyCrustOpen(mesh_endo.xyz);

rim_epi = find(mesh_epi.xyz(:,3)> max(mesh_epi.xyz(:,3))-0.1);
rim_endo= find(mesh_endo.xyz(:,3)> max(mesh_endo.xyz(:,3))-0.1);

figure()
hplotMESH(mesh_endo)
hold on
hplotMESH(mesh_epi)
scatter3(mesh_endo.xyz(rim_endo,1),mesh_endo.xyz(rim_epi,2),mesh_endo.xyz(rim_epi,3));

mesh_lid = struct;
mesh_lid.xyz=vertcat(mesh_epi.xyz(rim_epi,:),mesh_endo.xyz(rim_endo,:));

figure()
scatter3(mesh_lid.xyz(:,1),mesh_lid.xyz(:,2),mesh_lid.xyz(:,3))
mesh_lid.tri=delaunayTriangulation(mesh_lid.xyz(:,1:2));

mesh_combined=struct;
mesh_combined.xyz=vertcat(mesh_epi.xyz,mesh_endo.xyz,mesh_lid.xyz);
mesh_combined.tri=vertcat(mesh_epi.tri,mesh_endo.tri+size(mesh_epi.xyz,1));

idx_top_unclosed=find(abs(mesh_combined.xyz(:,3)-max(mesh_combined.xyz(:,3)))<0.1);





%look through tris and elim any tri that has any of the indices
%in eliminate ids



idx=  find(abs(mesh_combined.xyz(:,3)-max(mesh_combined.xyz(:,3)))<0.1);

figure()
hplotMESH(mesh_combined)
hold on
scatter(mesh_combined.xyz(idx,1),mesh_combined.xyz(idx,2),mesh_combined.xyz(idx,3))


lid=struct;
lid.xyz=mesh_combined.xyz(idx,:);
[p,tri]=MyCrustOpen2(mesh_combined.xyz(idx,:));


%only keep nodes ad tri of existing elements
pkeep=p(1:size(idx,1))';
sz=size(idx,1);
tri_keep=[];
for i=1:size(tri,1)
    keep_row=tri(i,1)<sz && tri(i,2)<sz && tri(i,3)<sz ;
    if keep_row==1
        tri_keep=vertcat(tri_keep,tri(i,:));
    end
end
lid.tri=tri_keep;



figure()
hplotMESH(lid)
%eliminate inner tris

radii = rs(mesh_combined.xyz(:,1),mesh_combined.xyz(:,2));
outer = find(radii(idx) >min(radii(idx))+0.1);
inner=  find(radii(idx) <min(radii(idx))+0.1);

tri_keep2=[];
for i=1:size(tri_keep,1)
    keep_row(i)=isempty(find(tri_keep(i,1)==inner)) + isempty(find(tri_keep(i,2)==inner)) + isempty(find(tri_keep(i,3)==inner)) ;
   
    if keep_row(i)>0
        tri_keep2=vertcat(tri_keep2,tri_keep(i,:));
    end
end

lid.tri=tri_keep2;

figure()
hplotMESH(lid)

mesh_final=mesh_combined;

mesh_final.tri=vertcat(mesh_combined.tri,idx(lid.tri));

figure()
plotMESH(mesh_final)
EDGE_LENGTH=0.4;
  [~,M] = jigsaw_surface_2_volume( mesh_final , 'delfront' , 'absolute','geom_feat',true,'hfun_hmax',EDGE_LENGTH,'hfun_hmin',EDGE_LENGTH*0.9);

  figure()
plotMESH(M)

  
Htetgen = gmshfn(M,['v']);
figure()
hplotMESH(Htetgen)


final=MeshFixCellOrientation(Htetgen);


MeshQuality()

write_VTK(Htetgen,'D:\meshes\ellipsoid.vtk','binary')

