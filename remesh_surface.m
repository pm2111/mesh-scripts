addpath(genpath('C:\Users\petnov\Dropbox\shared\'));
M = read_CHASTE('D:\meshes\data\2Dmesh');
M.celltype =5;

%the mesh is currently in decimeters, convert to mm for remeshing
M.xyz = M.xyz*100;
M.xyz(:,3) = zeros(size(M.xyz,1),1);
L=M;
for EDGE_LENGTH = [ 10 5 3 1  0.8 0.4]
  
  fprintf('now EDGE_LENGTH = %g\n', EDGE_LENGTH );
  L=M;

  [~,V] = jigsaw_surface_2_volume2D( L , 'delfront' , 'absolute','geom_feat',true,'hfun_hmax',EDGE_LENGTH,'hfun_hmin',EDGE_LENGTH*0.9,'verbose');
  L=V;

end
% M.xyz =M.xyz/10; %convert distances to cm for chaste
% 
% M.celltype =3;
% 
% V= struct;
% V.celltype=3;
% V.xyz = M.xyz(:,1:2);
% V.tri = M.tri;

figure()
scatter(M.xyz(M.xyzISBOUNDARY,1)*10,M.xyz(M.xyzISBOUNDARY,2)*10)
hold on
plotMESH(V)
%find boundary nodes

bounday_points_upscaled=vtkClosestPoint(Mesh(V),M.xyz(M.xyzISBOUNDARY,1:2));
if 1
figure()
scatter(M.xyz(M.xyzISBOUNDARY,1)*10,M.xyz(M.xyzISBOUNDARY,2)*10)
hold on
scatter(V.xyz(bounday_points_upscaled,1),V.xyz(bounday_points_upscaled,2))
plotMESH(V)
end
Z=struct;
Z.xyz = V.xyz(:,1:2)/10;
Z.tri = V.tri;
write_CHASTE(Z,'D:\meshes\data\2Dmeshcm')

% %For LV
% E = meshEdges( M );
%   fid = fopen( Fullfile( 'chaste' , sprintf('%s.edges',E) ),'w');
%   fprintf( fid , '%u\n' , size(xv_edges_reindexed,1) );
%   fprintf(  fid , '%u\t%u \n', 0:size(xv_edges_reindexed,1)-1,xv_edges_reindexed'-1 ); % 0 based!!
%   fclose( fid );


E =meshEdges(Z);
  xv_edges = E( all(E,2) ,:); % Node indices for edges with both nodes in xv.
  xv_edge_nodes = unique(xv_edges); % List of nodes found in xV edges (less than original file!)
% 
%   fid = fopen( strcat('D:\meshes\data\2Dmeshcm' , sprintf('%s_edge_nodes',xv))  ,'w' );
%   fprintf( fid , '%u\n' , xv_edge_nodes-1 ); % 0 based!!
%   fclose( fid );

  % Make node indices contiguous
  [~, xv_edges_reindexed] = ismember( xv_edges , xv_edge_nodes );
a = [transpose(0:size(xv_edges_reindexed,1)-1),xv_edges_reindexed-1].'; %strange way for fprintf to print in the right order
  fid = fopen('D:\meshes\data\2Dmeshcm.edge'  ,'w');
  fprintf( fid , '%u %u\n' , [size(xv_edges_reindexed,1) 0] );
  fprintf(  fid , '%u %u %u \n', a ); % 0 based!!
  fclose( fid );
% 
%   EL = sqrt( sum( ( H.xyz( xv_edges(:,2) ,:) - H.xyz( xv_edges(:,1) ,:) ).^2 ,2) );
% 
%   fid = fopen( Fullfile( 'chaste' , sprintf('%s.weights',xv) ),'w');
%   fwrite( fid , mm2cm( EL ) , 'double');
%   fclose( fid );


