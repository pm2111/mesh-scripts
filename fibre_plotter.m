addpath(genpath('C:/users/petnov/Dropbox/'))

H_mesh= read_CHASTE('D:/ARVC meshing automatic/patients/patient01/chaste01/HEART');
H_mesh.celltype=10;

 fibres=dlmread('D:/ARVC meshing automatic/patients/patient01/chaste01/HEART.ortho','\t',1,0);
 mesh=struct;
mesh.celltype=10;
mesh.xyz=H_mesh.xyz;
mesh.tri=H_mesh.tri;
mesh.triORTHO=[fibres(:,1:3)];
H_mesh.xyzFibers=fibres(:,1:3);
write_VTK(mesh,'D:/ARVC meshing automatic/patients/patient01/chaste01/mesh.vtk','binary')

H_mesh.celltype =10;

Hsmall=read_VTK('D:\ARVC meshing automatic\patients\patient01\HEART_0.20.vtk');
Hsmall.xyz = Hsmall.xyz/10;

epicardium = load('D:/ARVC meshing automatic/patients/patient01/chaste01/HEART.epi');
rv= load('D:/ARVC meshing automatic/patients/patient01/chaste01/HEART.rv');
lv= load('D:/ARVC meshing automatic/patients/patient01/chaste01/HEART.lv');

%coords = Hsmall.xyz(Hsmall.tri(:,1),:);

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

fiber_nodes = ones(size(mesh.xyz,1),1);
fiber_nodes(mesh.tri(:,1)) = angle(H_mesh.tri(:,1));
fiber_nodes(mesh.tri(:,2)) = angle(mesh.tri(:,2));
fiber_nodes(mesh.tri(:,3)) = angle(mesh.tri(:,3));
fiber_nodes(mesh.tri(:,4)) = angle(mesh.tri(:,4));

%need 1 datapoint per vertex!! find a way to assign these efficiently!

%write_VTK(mesh_vtk,'/data/temp/louieORdATs_ordered.vtk','binary');
 figure()
 camlight
patch('Faces',mesh.face,'Vertices',mesh.xyz,'FaceVertexCData',fiber_nodes(P),'FaceColor','interp','EdgeColor','none');
%patch('Faces',mesh.face,'Vertices',mesh.xyz,'FaceColor',colors,'EdgeColor','none');