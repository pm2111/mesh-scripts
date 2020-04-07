try
addpath(genpath('/home/scratch/petnov/Dropbox/shared - copy/IO/'));
addpath(genpath('/home/scratch/petnov/Dropbox/shared - copy/'));
addpath(genpath('/home/scratch/petnov/Dropbox/mesh scripts/'));
catch
    addpath(genpath('C:\Users\petnov\Dropbox\shared - copy\'));
addpath(genpath('C:\Users\petnov\Dropbox\mesh scripts'))
addpath(genpath('C:\Users\peter\Dropbox\mesh scripts'))
addpath(genpath('C:\Users\peter\Dropbox\shared - copy'))
end
coords=[];
for i=0:50
    for j=0:20
        
            coords=vertcat(coords,[ones(21,1)*0.04*i,ones(21,1)*0.04*j,linspace(0,0.8,21)']);
        
    end
end

%surf = meshgrid(coords(:,1),coords(:,2),coords(:,3));

coords_surf1=[];
i=0;

    for j=0:19
        
            coords_surf1=vertcat(coords_surf1,[ones(20,1)*0.04*i,ones(20,1)*0.04*j,linspace(0,0.8,20)']);
        
    end
    

    
mesh1=struct;
mesh1.xyz=coords_surf1;

tris=delaunay(coords_surf1(:,2:3));
mesh1.tri=tris;

figure()
hplotMESH(mesh1)
    coords_surf2=[];
i=49;

    for j=0:19
        
            coords_surf2=vertcat(coords_surf2,[ones(20,1)*0.04*i,ones(20,1)*0.04*j,linspace(0,0.8,20)']);
        
    end
    
mesh2=struct;
mesh2.xyz=coords_surf2;
mesh2.tri=delaunay(coords_surf2(:,2:3));
figure()
hplotMESH(mesh2)


coords_surf3=[];
for i=0:49
    for j=0:19
        
            coords_surf3=vertcat(coords_surf3,[0.04*i,0.04*j,0]);
        
    end
end

mesh3=struct;
mesh3.xyz=coords_surf3;
mesh3.tri=delaunay(coords_surf3(:,1:2));
figure()
hplotMESH(mesh3)

coords_surf4=[];
for i=0:49
    for j=0:19
        
            coords_surf4=vertcat(coords_surf4,[0.04*i,0.04*j,1*0.8]);
        
    end
end

mesh4=struct;
mesh4.xyz=coords_surf4;
mesh4.tri=delaunay(coords_surf4(:,1:2));
figure()
hplotMESH(mesh4)



coords_surf5=[];
j=0;
for i=0:49
        
            coords_surf5=vertcat(coords_surf5,[ones(20,1)*0.04*i,ones(20,1)*0.04*j,linspace(0,0.8,20)']);
        
    
end

mesh5=struct;
mesh5.xyz=coords_surf5;
mesh5.tri=delaunay(coords_surf5(:,1),coords_surf5(:,3));
figure()
hplotMESH(mesh5)


coords_surf6=[];
j=19;
for i=0:49
        
            coords_surf6=vertcat(coords_surf6,[ones(20,1)*0.04*i,ones(20,1)*0.04*j,linspace(0,0.8,20)']);
        
    
end

mesh6=struct;
mesh6.xyz=coords_surf6;
mesh6.tri=delaunay(coords_surf6(:,1),coords_surf6(:,3));
figure()
hplotMESH(mesh6)




sz1=size(mesh1.xyz,1);
sz2=size(mesh2.xyz,1);
sz3=size(mesh3.xyz,1);
sz4=size(mesh4.xyz,1);
sz5=size(mesh5.xyz,1);
sz6=size(mesh6.xyz,1);

surf=struct;
surf.xyz=vertcat(mesh1.xyz,mesh2.xyz,mesh3.xyz,mesh4.xyz,mesh5.xyz,mesh6.xyz);
surf.tri=vertcat(mesh1.tri,mesh2.tri+sz1,mesh3.tri+sz1+sz2,mesh4.tri+sz1+sz2+sz3,mesh5.tri+sz1+sz2+sz3+sz4,mesh6.tri+sz1+sz2+sz3+sz4+sz5);

figure()
hplotMESH(surf)



Htetgen = gmshfn(surf,['v']);



figure()
scatter3(coords_surf1(:,1),coords_surf1(:,2),coords_surf1(:,3))


mesh=struct;
mesh.xyz=coords;
mesh.tri=delaunay(coords);

model=createpde(1);
geometryFromMesh(model,mesh.xyz',mesh.tri');
generateMesh(model,'Hmax',0.04,'Hmin',0.03,'GeometricOrder','linear')

mesh=struct;
mesh.tri=model.Mesh.Elements';
mesh.xyz=model.Mesh.Nodes';
figure()
hplotMESH(mesh)

model = createpde(1);
model.Geometry=coords;


constraint=[];
for i=0:980
    constraint = [constraint,vertcat(1+(20*i)+0:20+(20*i), 2+(20*i)+0:21+(20*i),1+(20*(i+1))+0:20+(20*(i+1)))];
end
tri= delaunayTriangulation(coords,contraint');



%mesh =read_CHASTE('/data/meshes/mesh0.4mm/slab_0.4mm.ele');
sz = size(mesh.tri,1);
fid = fopen(strcat('/data/meshes/mesh0.4mm/','slab_0.4mm','.ortho'),'w');
fprintf(fid,'%d \n',sz);

orthos = vertcat(ones(1,sz),zeros(1,sz),zeros(1,sz),zeros(1,sz),ones(1,sz),zeros(1,sz),zeros(1,sz),zeros(1,sz),ones(1,sz));   

fprintf(fid,'%d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \n',orthos)

fclose(fid)

%find elements on faces
%first find nodes on faces

x= max(mesh.xyz);


face1 = mesh.xyz(:,1)==0;
face2 = mesh.xyz(:,1)==x(1);
face3= mesh.xyz(:,2)==0;
face4=mesh.xyz(:,2)==x(2);
face5=mesh.xyz(:,3)==0;
face6=mesh.xyz(:,3)==x(3);

%do each face individually


all_faces = face1;
faces_full =all_faces>0;

score = faces_full(mesh.tri(:,1))+faces_full(mesh.tri(:,2))+faces_full(mesh.tri(:,3))+faces_full(mesh.tri(:,4));
idx= [ faces_full(mesh.tri(:,1)), faces_full(mesh.tri(:,2)), faces_full(mesh.tri(:,3)), faces_full(mesh.tri(:,4))];
mesh_faces = score>2;
mesh_faces_line=mesh.tri(mesh_faces,:);
idx_faces=idx(mesh_faces,:);
mesh_faces_tri=[];
for i=1:size(mesh_faces_line,1)
    mesh_faces_tri = vertcat(mesh_faces_tri,mesh_faces_line(i,find(idx_faces(i,:)==1)));
end
meshtest=struct;
meshtest.xyz=mesh.xyz;
meshtest.tri=mesh_faces_tri;

figure()
hplotMESH(meshtest)



all_faces = face2;
faces_full =all_faces>0;

score = faces_full(mesh.tri(:,1))+faces_full(mesh.tri(:,2))+faces_full(mesh.tri(:,3))+faces_full(mesh.tri(:,4));
idx= [ faces_full(mesh.tri(:,1)), faces_full(mesh.tri(:,2)), faces_full(mesh.tri(:,3)), faces_full(mesh.tri(:,4))];
mesh_faces = score>2;
mesh_faces_line=mesh.tri(mesh_faces,:);
idx_faces=idx(mesh_faces,:);
for i=1:size(mesh_faces_line,1)
    mesh_faces_tri = vertcat(mesh_faces_tri,mesh_faces_line(i,find(idx_faces(i,:)==1)));
end

meshtest.tri=mesh_faces_tri;
figure()
hplotMESH(meshtest)

all_faces = face3;
faces_full =all_faces>0;

score = faces_full(mesh.tri(:,1))+faces_full(mesh.tri(:,2))+faces_full(mesh.tri(:,3))+faces_full(mesh.tri(:,4));
idx= [ faces_full(mesh.tri(:,1)), faces_full(mesh.tri(:,2)), faces_full(mesh.tri(:,3)), faces_full(mesh.tri(:,4))];
mesh_faces = score>2;
mesh_faces_line=mesh.tri(mesh_faces,:);
idx_faces=idx(mesh_faces,:);
for i=1:size(mesh_faces_line,1)
    mesh_faces_tri = vertcat(mesh_faces_tri,mesh_faces_line(i,find(idx_faces(i,:)==1)));
end

meshtest.tri=mesh_faces_tri;
figure()
hplotMESH(meshtest)

all_faces = face4;
faces_full =all_faces>0;

score = faces_full(mesh.tri(:,1))+faces_full(mesh.tri(:,2))+faces_full(mesh.tri(:,3))+faces_full(mesh.tri(:,4));
idx= [ faces_full(mesh.tri(:,1)), faces_full(mesh.tri(:,2)), faces_full(mesh.tri(:,3)), faces_full(mesh.tri(:,4))];
mesh_faces = score>2;
mesh_faces_line=mesh.tri(mesh_faces,:);
idx_faces=idx(mesh_faces,:);
for i=1:size(mesh_faces_line,1)
    mesh_faces_tri = vertcat(mesh_faces_tri,mesh_faces_line(i,find(idx_faces(i,:)==1)));
end


all_faces = face5;
faces_full =all_faces>0;

score = faces_full(mesh.tri(:,1))+faces_full(mesh.tri(:,2))+faces_full(mesh.tri(:,3))+faces_full(mesh.tri(:,4));
idx= [ faces_full(mesh.tri(:,1)), faces_full(mesh.tri(:,2)), faces_full(mesh.tri(:,3)), faces_full(mesh.tri(:,4))];
mesh_faces = score>2;
mesh_faces_line=mesh.tri(mesh_faces,:);
idx_faces=idx(mesh_faces,:);
for i=1:size(mesh_faces_line,1)
    mesh_faces_tri = vertcat(mesh_faces_tri,mesh_faces_line(i,find(idx_faces(i,:)==1)));
end


all_faces = face6;
faces_full =all_faces>0;mesh=struct;
mesh.xyz=p;
mesh.tri=tri;

score = faces_full(mesh.tri(:,1))+faces_full(mesh.tri(:,2))+faces_full(mesh.tri(:,3))+faces_full(mesh.tri(:,4));
idx= [ faces_full(mesh.tri(:,1)), faces_full(mesh.tri(:,2)), faces_full(mesh.tri(:,3)), faces_full(mesh.tri(:,4))];
mesh_faces = score>2;
mesh_faces_line=mesh.tri(mesh_faces,:);
idx_faces=idx(mesh_faces,:);
for i=1:size(mesh_faces_line,1)
    mesh_faces_tri = vertcat(mesh_faces_tri,mesh_faces_line(i,find(idx_faces(i,:)==1)));
end
meshf = struct;
meshf.xyz= coords;
meshf.tri=tri.ConnectivityList;
mesh.face=mesh_faces_tri;
meshf.face=mesh_faces_tri;

meshf.Celltype=10;
counter=0;
stim = [mesh.xyz(:,1)<0.5 , mesh.xyz(:,2)<0.5 ,mesh.xyz(:,3)<0.5] ;

meshtest=struct;
meshtest.xyz=mesh.xyz;
meshtest.tri=mesh_faces_tri;

figure()
hplotMESH(meshtest)

write_CHASTE(meshf,'/data/meshes/mesh0.4mm/slab_0.4mm','ascii')

fid = fopen(strcat('/data/meshes/mesh0.4mm/','slab_0.4mm','.face'),'w');
fprintf(fid,'%d %u \n',[size(mesh.face,1), 0])
fprintf(fid,'%d %d %d %d \n',vertcat([0:size(mesh.face,1)-1], mesh.face'-1))


test= read_CHASTE('/data/meshes/mesh0.2mm/slab_0.2mm.ele');
write_VTK(test,'/data/meshes/mesh0.2mm/test.vtk','binary')


test.tri=test.face;
figure()
hplotMESH(test)
