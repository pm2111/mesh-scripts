%load the heart squared file that you wish to plot!! 
addpath(genpath('D:/shared/'));

patients = {'02','05','09','18'};
distances = cell(1);
for z = 1:4
path = strcat('D:/ARVC meshing automatic/patients/patient',patients{z},'/mpp/');
cd (path);
load('HQ.mat');
load('HM.mat');


col = ['r','g','bl','m','b','c'];
figure()
hold all
for i=1:size(HQ,1)
    for j=1:6
            try
                scatter3(HQ{i,j}(:,1),HQ{i,j}(:,2),HQ{i,j}(:,3),50,'filled',col(j));
            end
    end
end
grid off
axis on
hplotMESH( RVms);



%HQ = HF;

% mesh = read_VTK('D:/ARVC meshing automatic/patients/patient09/HEART.vtk');
% 
% HM =mesh;
figure()
hold all
hplotMESH( HM);
hplotMESH( EPIms , 'b' ,'FaceAlpha',0.3,'ne');

figure()
hold all
hplotMESH( RVms , 'b' ,'FaceAlpha',0.3,'ne');


%algo to calculate the discrepancy between a mesh and a contour point

%  (1) get rid of contour lines above the mesh base
HQ_trimmed = cell(1);

%create plane by selecting 3 points on base of heart
normal = cross([p1.Position(1) - p2.Position(1),p1.Position(2)-p2.Position(2),p1.Position(3)-p2.Position(3)] , [p1.Position(1) - p3.Position(1),p1.Position(2)-p3.Position(2),p1.Position(3)-p3.Position(3)]);

plane_eq = @(x,y,z)normal(1)*(x  - p1.Position(1)) + normal(2)*(y - p1.Position(2)) + normal(3)*(z - p1.Position(3));
for i=1:size(HQ,1)
    for j=2:size(HQ,2)
%check if equation of plane will yield a +ve or -ve value for a known
%point:
%for each contour, check the value of the plane_eq func
        try
            plane_res = plane_eq(HQ{i,j}(:,1),HQ{i,j}(:,2),HQ{i,j}(:,3));

            if plane_eq(p4.Position(1),p4.Position(2),p4.Position(3)) > 0
                smaller_than_zero_idx = find(plane_res > 0);

                HQ_trimmed{i,j-1} = HQ{i,j}(smaller_than_zero_idx,:);
            else
                smaller_than_zero_idx = find(plane_res < 0);

                HQ_trimmed{i,j-1} = HQ{i,j}(smaller_than_zero_idx,:);
            end
        end
    end
end

%see if cutting out contour points works

col = ['r','g','bl','m','b','c'];
figure()
hold all
for i=1:size(HQ_trimmed,1)
    for j=1:6
            try
                scatter3(HQ_trimmed{i,j}(:,1),HQ_trimmed{i,j}(:,2),HQ_trimmed{i,j}(:,3),50,'filled',col(j));
            end
    end
end
grid off
axis on
hplotMESH( RVms);

%find the closest element of the mesh to a contour point
euclidian_dist = @(p1,p2)((p1(:,1)-p2(:,1)).^2 + (p1(:,2)-p2(:,2)).^2 +(p1(:,3)-p2(:,3)).^2).^0.5;
for i=1:size(HQ_trimmed,1)
    for j=1:size(HQ_trimmed,2)
        if j==1
            try
            closest_ind_mesh = knnsearch(EPIms.xyz, HQ_trimmed{i,j});

            distances{z,1} = euclidian_dist(EPIms.xyz(closest_ind_mesh,:),HQ_trimmed{i,j});
            end
        end
        if j==2
            try
                closest_ind_mesh = knnsearch(LVms.xyz, HQ_trimmed{i,j});

                distances{z,2} = euclidian_dist(LVms.xyz(closest_ind_mesh,:),HQ_trimmed{i,j});
            end
        end
        if j==3
            try
                closest_ind_mesh = knnsearch(RVms.xyz, HQ_trimmed{i,j});

                distances{z,3} = euclidian_dist(RVms.xyz(closest_ind_mesh,:),HQ_trimmed{i,j});
            end
        end
         if j==5
            try
                closest_ind_mesh = knnsearch(RVms.xyz, HQ_trimmed{i,j});

                distances{z,4} = euclidian_dist(RVms.xyz(closest_ind_mesh,:),HQ_trimmed{i,j});
            end
        end
        
        end
end

end
