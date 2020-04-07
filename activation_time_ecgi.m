
%read ECGi input and compute the time until activation of the site


data = load('/data/ChrisReconstructions/ARVC18/ARVC18_Baseline.mat');

potentials = data.arvc18_baseline.epipots.potential;
times = data.arvc18_baseline.epipots.t;


for i= 1:size(potentials,1)

    [mx,ind] = max(potentials(i,:));

    time_to_peak(i) = times(ind);

end

h = figure();
patch('vertices', data.arvc18_baseline.nodes, 'faces', data.arvc18_baseline.mesh, 'facecolor', 'interp', 'cdata', time_to_peak, 'edgecolor', 'none') ;
colorbar()

data = load('/data/ChrisReconstructions/ARVC14/ARVC14_Baseline.mat');

potentials = data.arvc14_baseline.epipots.potential;
times = data.arvc14_baseline.epipots.t;


for i= 1:size(potentials,1)

    [mx,ind] = max(potentials(i,:));

    time_to_peak(i) = times(ind);

end

h = figure();
patch('vertices', data.arvc2_baseline.nodes, 'faces', data.arvc2_baseline.mesh, 'facecolor', 'interp', 'cdata', time_to_peak, 'edgecolor', 'none') ;
colorbar()

data = load('/data/ChrisReconstructions/ARVC5/ARVC5_Baseline.mat');

potentials = data.arvc5_baseline.epipots.potential;
times = data.arvc5_baseline.epipots.t;


for i= 1:size(potentials,1)

    [mx,ind] = max(potentials(i,:));

    time_to_peak(i) = times(ind);

end

h = figure();
patch('vertices', data.arvc5_baseline.nodes, 'faces', data.arvc5_baseline.mesh, 'facecolor', 'interp', 'cdata', time_to_peak, 'edgecolor', 'none') ;
colorbar()







data = load('/data/ChrisReconstructions/ARVC9/ARVC9_Baseline.mat');

potentials = data.arvc9_baseline.epipots.potential;
times = data.arvc9_baseline.epipots.t;


for i= 1:size(potentials,1)

    [mx,ind] = max(potentials(i,:));

    time_to_peak(i) = times(ind);

end

h = figure();
patch('vertices', data.arvc9_baseline.nodes, 'faces', data.arvc9_baseline.mesh, 'facecolor', 'interp', 'cdata', time_to_peak, 'edgecolor', 'none') ;
colorbar()






data = load('/data/ChrisReconstructions/ARVC2/ARVC2_Baseline.mat');

potentials = data.arvc2_baseline.epipots.potential;
times = data.arvc2_baseline.epipots.t;


for i= 1:size(potentials,1)

    [mx,ind] = max(potentials(i,:));

    time_to_peak(i) = times(ind);

end

h = figure();
patch('vertices', data.arvc2_baseline.nodes, 'faces', data.arvc2_baseline.mesh, 'facecolor', 'interp', 'cdata', time_to_peak, 'edgecolor', 'none') ;
colorbar()

%look at the potentials of the affected regions
[temp, originalpos] = sort(time_to_peak,'descend');


[mx, idx] = max(time_to_peak);

figure()
hold on
plot(potentials(idx,:))
plot(potentials(originalpos(2),:))
plot(potentials(originalpos(3),:))
plot(potentials(originalpos(50),:))




data = load('/data/ChrisReconstructions/ARVC2/ARVC2_Post_Exercise.mat');

potentials = data.arvc2_post_exercise.epipots.potential;
times = data.arvc2_post_exercise.epipots.t;



h = figure();
patch('vertices', data.arvc2_post_exercise.nodes, 'faces', data.arvc2_post_exercise.mesh, 'facecolor', 'interp', 'cdata', time_to_peak, 'edgecolor', 'none') ;
colorbar()


data = load('/data/ChrisReconstructions/ARVC2/ARVC2_Baseline.mat');


potentials = data.arvc2_baseline.epipots.potential;
times = data.arvc2_baseline.epipots.t;

data = load('/data/ChrisReconstructions/ARVC5/ARVC5_Baseline.mat');


potentials = data.arvc5_baseline.epipots.potential;
times = data.arvc5_baseline.epipots.t;

data = load('/data/ChrisReconstructions/ARVC18/ARVC18_Baseline.mat');


potentials = data.arvc18_baseline.epipots.potential;
times = data.arvc18_baseline.epipots.t;

data = load('/data/ChrisReconstructions/ARVC9/ARVC9_Baseline.mat');


potentials = data.arvc9_baseline.epipots.potential;
times = data.arvc9_baseline.epipots.t;



data = load('/data/ChrisReconstructions/ARVC2/ARVC2_Baseline.mat');


potentials = data.arvc2_baseline.epipots.potential;
times = data.arvc2_baseline.epipots.t;

%%%find plot the voltage regions of > 1mV (fibrotic)
for i= 1:size(potentials,1)

    [mx,ind] = max(potentials(i,:)); %find the max voltage of each heart electrode over time

    time_to_peak(i) = times(ind); 
    
    potential_amplitudes(i) = mx;

end


fibrotic_amplitudes = potential_amplitudes <0.3*max(potential_amplitudes); %constant is of arbitrary value
figure()
hist(fibrotic_amplitudes)

%find number of large scale deflections in electrogram
for j=1:size(potentials,1)
    y = potentials(j,:);
    dy = diff(y);
    sign = dy >0;
    counter =0;
    for i=1:size(sign,2)-1
       if isequal(   sign(i+1)  ,sign(i)) ==0  && potentials(j,i) > 0.2*max(potential_amplitudes(j)) %0.2 is an empirical constant
          counter = counter+1;
       end
    
    end
    fractionation_nr(j) =counter;
end

figure()
hist(fractionation_nr)

fractionation = fractionation_nr >15; %empirical number

fibrotic_patch = fibrotic_amplitudes + fractionation;
fibrotic = fibrotic_patch ==2;


%apical_anterior_rv = data.arvc9_baseline.labels.apical_anterior_rv;
% 
% lvaa = transpose(zeros(size(potentials,1),1));
% lvaa(apical_anterior_rv) =1; %useful for orientation in the meshes

figure()
hold all

patch('vertices', data.arvc2_baseline.nodes, 'faces', data.arvc2_baseline.mesh, 'facecolor', 'interp', 'cdata', int32(fibrotic), 'edgecolor', 'none') ;

patch('vertices', data.arvc2_baseline.nodes, 'faces', data.arvc2_baseline.mesh, 'facecolor', 'interp', 'cdata', lvaa, 'edgecolor', 'none') ;



colorbar()

