

vest=load('D:\ARVC meshing automatic\patients\patient09\mpp\fittedVEST.mat');
electrodes = load('D:\ARVC meshing automatic\patients\patient09\mpp\ECG_ELECTRODES.mat');

els = struct2cell(electrodes.ELECTRODES);
figure()
plotMESH(vest.VEST,'ne')
headlight
hold on
for i=1:12
    x=els{i};
    scatter3(x(1),x(2),x(3),40,'filled')
end