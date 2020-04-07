%ecg eikonal
addpath(genpath('C:/Users/petnov/Dropbox/'));
enableVTK;


leads=csvread('D:\ARVC meshing automatic\patients\patient06\results\Eikonal\PseudoECG__EndoSpeed0.15_EpiSpeed1_pred_PseudoECG_pred.csv');
for i=1:10        
    if i<10
    fid = fopen(strcat('D:\ARVC meshing automatic\patients\patient06\results\Eikonal\','reordered\0',int2str(i),'.txt'),'w');
        fprintf(fid,'%d \n',leads(i,:)');
        fclose(fid)
    else
            fid = fopen(strcat('D:\ARVC meshing automatic\patients\patient06\results\Eikonal\','reordered\',int2str(i),'.txt'),'w');
        fprintf(fid,'%d \n',leads(i,:)');
        fclose(fid)
    end
end
[I,II,III,aVR,aVL,aVF,V1w,V2w,V3w,V4w,V5w,V6w,LA]=my_ecg_from_sim_eikonal('06','D:\ARVC meshing automatic\patients\patient06\results\Eikonal\',62);

figure()
subplot(4,3,1)
plot(I)
subplot(4,3,2)
plot(II)
subplot(4,3,3)
plot(III)
subplot(4,3,4)
plot(aVR)
subplot(4,3,5)
plot(aVL)
subplot(4,3,6)
plot(aVF)