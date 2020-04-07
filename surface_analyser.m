%this script was written by Peter Marinov (peter.marinov@cs.ox.ac.uk) and
%takes a matlab ventricular surface mesh file as input
%OUTPUTS:
% - ventricular endocardial volume 


filename = '/home/scratch/meshes ready/patient05/HM.mat';

data = open(filename);

%visualise the data
figure()
title('Full surface mesh')
patch('vertices', data.HM.xyz, 'faces', data.HM.tri) ;
camlight('headlight')
alpha(0.3)  

figure()
title('Epicardial mesh')
patch('vertices', data.EPIms.xyz, 'faces', data.EPIms.tri) ;
camlight('headlight')
alpha(0.3)  


figure()
title('Endo RV mesh')
patch('vertices', data.RVms.xyz, 'faces', data.RVms.tri) ;
camlight('headlight')
alpha(0.3) 


figure()
title('Endo LV mesh')
patch('vertices', data.LVms.xyz, 'faces', data.LVms.tri) ;
camlight('headlight')
alpha(0.3) 