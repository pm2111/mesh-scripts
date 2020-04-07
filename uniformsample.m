%Uniform density of points algorithm
nr_points =10000;
n = 5000; %final nr points required

%Generate the random sample of points in 2D
xs = -5 + (5+5)*rand(nr_points,1);
ys = -5 + (5+5)*rand(nr_points,1);


%parameters of the algorithm:
% specify the total initial number of points (nr_points)
% specify the points we want to end up with (n)

%compute the total volume (in this case area) of the cube

volume = (max(xs)-min(xs))*(max(ys)-min(ys)); %this is in 2d atm
%select the number of smaller cubes which will make the big cube (these
%need to be as close to equal size as possible)

nr_cubes = n; %important!! : the number of cubes may not be exactly equal
%to the desired amout ( error decreases with increasing number of points)
dx = sqrt(volume) / sqrt(nr_cubes);

% 
% dx = (max(xs)-min(xs))/nr_cubes;
% dy= (max(ys)-min(ys))/nr_cubes;


%compute the bottom left corner of the 1st cube edge 
cube_edges = [min(xs) min(ys)];
cube_centres = [min(xs)+dx/2  min(ys)+dx/2];

%compute the bottom left corner of the cube (square) in cartesian coordinates
for i=1:round(sqrt(nr_cubes))
    for j=1:round(sqrt(nr_cubes))
         cube_edges = vertcat(cube_edges, [cube_edges(1,1)+dx*(j-1) cube_edges(1,2)+dx*(i-1)]);
         cube_centres = vertcat(cube_centres, [cube_edges(1,1)+0.5*dx*(j-1) cube_edges(1,2)+0.5*dx*(i-1)]);

         
    end
end
%plot the square bottom left edges:
% cube_edges = cube_edges(2:end,2:end);
% cube_centres = cube_centres(2:end,2:end);

figure()
plot([cube_edges(:,1)],[cube_edges(:,2)],'o')




%place the points in 1 of the n cubes
coords_in_square = cell(1);
ids_in_square = cell(1);
id_sum =0;
for i=1:nr_cubes
    

    comparison1 = ones(nr_points,1);
    comparison2 = ones(nr_points,1);

    comparison1 = comparison1*cube_edges(i,1);  % bottom left x coord 
    comparison2 = comparison2*cube_edges(i,2);  % bottom right y coord

    %make an array of the top right coordinate, repeated for the length equiv
    %to nr_cubes
    comparison3 = comparison1 +dx; %top right x coord
    comparison4 = comparison2 +dx;  % top right y coord

    %see if all conditions are satisfied (need a score of 4 for a given point to lie within the square) 
    score = double(xs > comparison1) + double(ys > comparison2) + double(xs  < comparison3) + double(ys < comparison4);
    
    ids = find(score ==4);
    ids_in_square{i,1} = ids;
    coords_in_square{i,1} =xs(ids);
    coords_in_square{i,2} =ys(ids); %need to extend this in the z 
    
    id_sum = id_sum + size(ids,1);
 
end

%consistency check: The number of points should be conserved 
if id_sum ~= nr_points
    fprintf(strcat(strcat('Boundary points are ignored! The Error is ',num2str((nr_points - id_sum)*100/nr_points)),'percent'));
end
 
% find the point in each cube closest to the middle of the square 
for i=1:nr_cubes
    if size(ids_in_square{i,1},1) >1
        distances_from_centre{i,1} = coords_in_square{i,1} - cube_centres(i,1);
        distances_from_centre{i,2} = coords_in_square{i,2} - cube_centres(i,2);

        distances_squared{i,1} = power(distances_from_centre{i,1},2);
        distances_squared{i,2} = power(distances_from_centre{i,2},2);
        a =[distances_squared{i,1}(:)];
        b =[distances_squared{i,2}(:)];

        geometrical_distance_sq{i,1} = a+b;
        [mx,idx] = max(geometrical_distance_sq{1,1});
        ids_in_square{i,2} = ids_in_square{i,1}(idx);
    end
    if size(ids_in_square{i,1},1)  == 1
                ids_in_square{i,2} = ids_in_square{i,1};

    end
end

xs_selected = xs([ids_in_square{:,2}]);
ys_selected = ys([ids_in_square{:,2}]);


% plot the original points and the sampled points

figure()
hold all
plot(xs, ys, 'o')
plot(xs_selected,ys_selected,'o','Color','r')


