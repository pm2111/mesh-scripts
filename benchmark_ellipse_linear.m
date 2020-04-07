% Generates linear lagrange hexahedral meshes for the mechanics benchmark problems 2 and 3
% - cmgui .exnode/.exelem
% - simple text output:
%     .X for nodes
%     .T  for elements (zero-based node indices in xi-based order per element)
%     .F for fibers:  angle x y z per line in the same order as .X
% ".in" problem file for my own code, which may be useful to extract the dirichlet and pressure boundary conditions
%
% Input:
% h = approxite element size
% or h = [nab ncirc nr] where
% - nab:   number of elements apex-to-base
% - ncirc: number of elements circumferentially
% - nr:    number of elements transmurally
function [Cnodes,Celems,fib,naive_epi1endo0]=benchmark_ellipse_linear(h)

if nargin == 0
  h = 1;
end

if length(h)==3
  nab = h(1);
  ncirc = h(2);
  nr = h(3);
else
  assert(length(h)==1,'h must be one element or three, see function description');
  wt=0.7;
  top=2.0;
  r=2.9+wt;
  outer=6.13+0.7;
  ncirc = round(2*pi*(r-wt/2)/h);
  nab =  round(  ((outer/2+top)/(2*outer) * pi*sqrt((outer-wt)^2+r^2) ) / h);
  nr = round(wt/h);
end

output_tag = sprintf('ellipse_benchmark_lin_%d-%d-%d',nr,ncirc,nab);

% ellipse geometry
top = 5;

outer_long_axis = 20;
outer_short_axis = 10;
wall_thickness = 3;

% fiber angles
fiber_angle_epi  = -90;
fiber_angle_endo = 90;

% mesh generation: points to base nodes off

points_circ  = ncirc;
points_ab    = nab + 1;
points_trans = nr  + 1;

iz = 1;
coords = zeros(5,points_ab,points_circ+1,iz);

for long_r=(outer_long_axis-wall_thickness):(wall_thickness/(points_trans-1)):outer_long_axis
  short_r = outer_short_axis - (outer_long_axis-long_r);
  
  uu = -acos(top/long_r);
  u = -pi:(pi+uu)/(points_ab-1):uu;
  v = -pi:(2*pi/points_circ):pi;
  
  coords(1,:,:,iz) = short_r .* sin(u)' * cos(v);
  
  figure()
  coords(2,:,:,iz) = short_r .* sin(u)' * sin(v);
  coords(3,:,:,iz) = long_r .* cos(u)' * ones(1,length(v));
  coords(4,:,:,iz) = repmat(u',1,length(v));
  coords(5,:,:,iz) = repmat(v,length(u),1);
  
  iz = iz+1;
end
assert(points_trans == iz-1)

% mesh generation: number nodes naively first, then connect later

dof_per_element=8;
el_nodes = zeros(dof_per_element,5);
naive_elements = zeros(nab*ncirc*nr,dof_per_element);
naive_nodes    = zeros(dof_per_element*nab*ncirc*nr,5);
naive_epi1endo0 = zeros(dof_per_element*nab*ncirc*nr,1);

ei = 0;
for el_ab = 1:nab
  for el_circ = 1:ncirc
    for el_r = 1:nr
      
      ix_ab = (el_ab-1)   + (1:2);
      ix_c  = (el_circ-1) + (1:2);
      ix_r  = (el_r-1)    + (1:2);
      
      for xi3=1:2
        for xi2=1:2
          for xi1=1:2
            el_nodes(4*(xi3-1) + 2*(xi2-1) + xi1 ,:) = coords(:,ix_ab(xi2),ix_c(xi1),ix_r(xi3));
          end
        end
      end
      
      ei=ei+1;
      nn = (ei-1)*8 + (1:8);
      naive_elements(ei,:) = nn;
      naive_nodes(nn,:) = el_nodes;
      naive_epi1endo0(nn) = (el_r-1)/nr + (floor((0:7)/4))/nr;
    end
  end
end

% merge nodes to create mesh
[Celems,Cnodes,nodemap,revmap] = merge_nodes(naive_elements,naive_nodes(:,1:3));
UV = naive_nodes(revmap,4:5);

% determine angles in actual nodes and check they are consistent
epi1endo0 = zeros(size(Cnodes,1),1);
for i=1:size(Cnodes,1)
  f = naive_epi1endo0(nodemap==i);
  assert( max(f) - min(f) < 1e-10);
  epi1endo0(i) = f(1);
end

fib = zeros(size(Cnodes,1),3);
angle = zeros(size(Cnodes,1),1);
for ni=1:size(Cnodes,1)
  u = UV(ni,1); v = UV(ni,2);
  fiber_angle = (fiber_angle_endo + epi1endo0(ni) * (fiber_angle_epi - fiber_angle_endo)) * (pi/180);
  
  long_r  = outer_long_axis - (1-epi1endo0(ni))*wall_thickness;
  short_r = outer_short_axis - (outer_long_axis-long_r);
  
  deriv_dir = [sin(fiber_angle); cos(fiber_angle)];
  
  M = [short_r*cos(u)*cos(v) -short_r*sin(u)*sin(v);  % these are simply the d/du and d/dv in a matrix
    short_r*cos(u)*sin(v)  short_r*sin(u)*cos(v);
    -long_r*sin(u)          0                   ];
  M(:,1) = M(:,1) / norm(M(:,1)); % normalize directions
  M(:,2) = M(:,2) / norm(M(:,2));
  
  fib(ni,:) = M * deriv_dir;
  fib(ni,:) = fib(ni,:) / norm(fib(ni,:));
  
  angle(ni) = fiber_angle*(180/pi);
  if abs(u+pi)<1e-6 % apex - change this if you want another choice for apical fiber direction
    fib(ni,:) = 0;
    angle(ni) = 0;
  end
end

% ---------- MESH OUTPUT

fh =fopen(sprintf('%s.exelem',output_tag),'w');
fprintf(fh,'Group name: solution\nShape. Dimension=3, line*line*line\n #Scale factor sets=0\n #Nodes=8\n #Fields=1\n 1) coordinates, coordinate, rectangular cartesian, #Components=3\n');
for x='xyz'
  fprintf(fh,['  ' x '. l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.\n  #Nodes=8\n']);
  for i=1:8
    fprintf(fh,['  ' int2str(i) '. #Values=1\n     Value indices: 1\n     Scale factor indices: 0\n']);
  end
end

for e=1:size(Celems,1)
  fprintf(fh,'\nElement: %d 0 0\nNodes: %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d \n',e,Celems(e,:));
end
fclose(fh);

fh=fopen(sprintf('%s.exnode',output_tag),'w');
fprintf(fh,'Region: /\n #Fields=2\n 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
for x='xyz'
  fprintf(fh,['   ' x '.  Value index=1, #Derivatives=0, #Versions=1\n']);
end
fprintf(fh,'2) fibers, coordinate, rectangular cartesian, #Components=3\n')
for x='xyz'
  fprintf(fh,'   %s.  Value index=%d, #Derivatives=0, #Versions=1\n',x,x-'x'+4);
end

for ni=1:size(Cnodes,1)
  fprintf(fh,'\nNode: %d\n %.8f %.8f %.8f\t%.8f %.8f %.8f',ni,Cnodes(ni,:),fib(ni,:));
end
fclose(fh);
fprintf('Wrote exnode, exelem files for linear hexes to ./out/%s.*\n',output_tag);


% Merges identical nodes
% nn x 3 -> newnn x 3  . elements size irrelevant
function [elements,nodes,nodemap,reversenodemap] = merge_nodes(elements,nodes)
os = size(nodes,1);
nodes(:,[1 2 3]) = nodes(:,[3 2 1]);
[nodes,reversenodemap,nodemap] = unique(roundn(nodes,-10),'rows'); % with tolerence for rounding
nodes(:,[1 2 3]) = nodes(:,[3 2 1]);

elements = nodemap(elements);
fprintf('removed duplicates %d -> %d\n',os,size(nodes,1));

