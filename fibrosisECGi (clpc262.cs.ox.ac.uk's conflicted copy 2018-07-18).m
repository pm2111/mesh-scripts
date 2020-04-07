%fibrosis from ECGi

filenames = dir('/data/ChrisReconstructions/baselines/');

data = load('/data/ChrisReconstructions/ARVC18/ARVC18_Baseline.mat');

for k=1:size(filenames,1)-2
    data = load(strcat(filenames(k+2).folder,'/',filenames(k+2).name));
    fields = fieldnames(data);

    activation = data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart;

    
    
    potentials = data.(fields{1}).epipots.potential;
    times = data.(fields{1}).epipots.t;

    %%%find plot the voltage regions of > 1mV (fibrotic)
    for i= 1:size(potentials,1)

        [mx,ind] = max(potentials(i,:)); %find the max voltage of each heart electrode over time

        time_to_peak(i) = times(ind); 

        potential_amplitudes(i) = mx;

    end


    fibrotic_amplitudes = potential_amplitudes <0.3*max(potential_amplitudes); %constant is of arbitrary value
    % figure()
    % hist(fibrotic_amplitudes)
    % title('fibrotic amplitudes =1')

    %find number of large scale deflections in electrogram
    for j=1:size(potentials,1)
        y = potentials(j,:);
        dy = diff(y);
        sign = dy >0;
        counter =0;
        for i=1:size(sign,2)-1
           if isequal(   sign(i+1)  ,sign(i)) ==0  && potentials(j,i) > 0.2*max(potential_amplitudes(j))%0.2 is an empirical constant
              counter = counter+1;
           end

        end
        fractionation_nr(j) =counter;
    end

    % figure()
    % hist(fractionation_nr)
    % title('fractionation nr')
    fractionation = fractionation_nr > mean(fractionation_nr)+std(fractionation_nr); %empirical number

    fibrotic_patch = fibrotic_amplitudes + fractionation;
    fibrotic = fibrotic_patch ==2;


    % 
    % figure()
    % hold all
    % 
    % patch('vertices', dafibrotic_patchta.(fields{1}).nodes, 'faces', data.(fields{1}).mesh, 'facecolor', 'interp', 'cdata', int32(fibrotic), 'edgecolor', 'none') ;
    % title('patient 2')

    fibrotic_cohort(k,:) = fibrotic;
    activations_cohort(k,:) = activation;
end

%now find correlation coeff between the fibrotic regions and activation
%times
for i=1:20
    [dummy,p] =corrcoef([transpose(fibrotic_cohort(i,:)) transpose(activations_cohort(i,:))]);
    correlations(1,i) = dummy(1,2); %this shows very low correlation with the current formulation
    correlations(2,i) = p(1,2);
end
