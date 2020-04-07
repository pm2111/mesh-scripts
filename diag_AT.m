function [at]=diag_AT(mesh_path,sim_path)


permutation = load(strcat(sim_path,'permutation.txt'));
permutation = permutation+1;

fileinfo =h5info(strcat(sim_path,'results.h5'));
upstroke = h5read(strcat(sim_path,'results.h5'),strcat('/',fileinfo.Datasets(3).Name));
upstroke_heart = upstroke(1,permutation,1);

mesh=read_CHASTE(mesh_path);
[x]=max(mesh.xyz);

line = [linspace(0,1,100)*2;linspace(0,1,100)*0.8;linspace(0,1,100)*0.8];

idx=knnsearch(mesh.xyz,line');

at =upstroke_heart(idx);
