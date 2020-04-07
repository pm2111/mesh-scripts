function [principal_cv,inter_fiber_cv,inter_sheet_cv]=cv_calc(path,mesh)

permutation = load(strcat(path,'permutation.txt'));
permutation = permutation+1;

fileinfo =h5info(strcat(path,'results.h5'));
upstroke = h5read(strcat(path,'results.h5'),strcat('/',fileinfo.Datasets(3).Name));
upstroke_heart = upstroke(1,permutation,1);

[x,y] = max(mesh.xyz);

idx=find(mesh.xyz(:,1)==1.2 & mesh.xyz(:,2)==0 & mesh.xyz(:,3)==0);

idx2=find(mesh.xyz(:,1)==0.8 & mesh.xyz(:,2)==0 & mesh.xyz(:,3)==0);

principal_cv = (mesh.xyz(idx,1)-mesh.xyz(idx2,1))/(0.001*(upstroke_heart(idx)-upstroke_heart(idx2)));


idx=find(mesh.xyz(:,1)==0 & mesh.xyz(:,2)==0.2 & mesh.xyz(:,3)==0);

idx2=find(mesh.xyz(:,1)==0 & mesh.xyz(:,2)==0.6 & mesh.xyz(:,3)==0);

inter_fiber_cv = (mesh.xyz(idx,2)-mesh.xyz(idx2,2))/(0.001*(upstroke_heart(idx)-upstroke_heart(idx2)));


idx=find(mesh.xyz(:,1)==0 & mesh.xyz(:,2)==0 & mesh.xyz(:,3)==0.2);

idx2=find(mesh.xyz(:,1)==0 & mesh.xyz(:,2)==0 & mesh.xyz(:,3)==0.6);


inter_sheet_cv =(mesh.xyz(idx,3)-mesh.xyz(idx2,3))/(0.001*(upstroke_heart(idx)-upstroke_heart(idx2)));
