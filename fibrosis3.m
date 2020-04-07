%fibrosis from ECGi

filenames = dir('/data/ChrisReconstructions/baselines/');

data = load('/data/ChrisReconstructions/ARVC18/ARVC18_Baseline.mat');

for k=1:size(filenames,1)-2
    data = load(strcat(filenames(k+2).folder,'/',filenames(k+2).name));
    fields = fieldnames(data);

    activation = data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart;

    
    
    potentials = data.(fields{1}).epipots.potential;
    times = data.(fields{1}).epipots.t;
    dt = times(2) -times(1) ;%ms
    
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
        dy_dt = diff(y)/dt;
        qrs_start = abs(data.(fields{1}).epipots.t - data.(fields{1}).nodeparam.timestamps.activationStart(1)) < 0.5;
        
        qrs_end = abs(data.(fields{1}).epipots.t   - data.(fields{1}).nodeparam.timestamps.recoveryStart(1)) <0.5;
        qrs_start_idx = find(qrs_start,1);
        qrs_end_idx = find(qrs_end,1);

        qrs = y(qrs_start_idx:qrs_end_idx);
        qrs_time = times(qrs_start_idx:qrs_end_idx);
        qrs_dy_dx = dy_dt(qrs_start_idx:qrs_end_idx);
        local_maxg_idxs = islocalmin(qrs_dy_dx);
        
        steepest_qrs_deflections = qrs(local_maxg_idxs);
        
        figure()
        hold all
        plot(qrs_time,qrs)
        plot(qrs_time(local_maxg_idxs),steepest_qrs_deflections,'r*')
        
        dummy = linspace(1,size(potentials,2),size(potentials,2));
        local_maxg_idxs = dummy(local_maxg_idxs);
        local_max_idxs = dummy(local_max_idxs);
        
        baseline = mean(y);
        baseline = abs(baseline);
        
        counter =0;
        min_idx_selected = [];
        max_idx_selected = [];
        for i=2:size(steepest_qrs_deflections,2)-1
%             min_qualified = 0;
%             max_qualified =0;
%             if abs(min_amplitudes(i)-max_amplitudes(i)) > 0.3*potential_amplitudes(j)  %check if minima are large enough
%               counter = counter+1;
%               min_idx_selected = [min_idx_selected,local_min_idxs(i)];
%               qualified =1;
%             end
%             if qualified == 0
%                 if abs(max_amplitudes(i)-min_amplitudes(i-1)) > 0.3*potential_amplitudes(j)  %check if maxima are large enough
%                   counter = counter+1;
%                   max_idx_selected = [max_idx_selected,local_max_idxs(i)];
%                 end
%             
%             end
%             if qualified ==1
%                  if abs(max_amplitudes(i)) > 0.3*potential_amplitudes(j)  %check if maxima are large enough
%                   counter = counter+1;
%                   max_idx_selected = [max_idx_selected,local_max_idxs(i)];
%                  end
%             end
% 
            if abs(steepest_qrs_deflections(i))-baseline > 0.3*potential_amplitudes(j)  %check if minima are large enough
              counter = counter+1;
              min_idx_selected = [min_idx_selected,local_maxg_idxs(i)];
            end

            if abs(max_amplitudes(i)-baseline) > 0.3*potential_amplitudes(j)  %check if maxima are large enough
              counter = counter+1;
              max_idx_selected = [max_idx_selected,local_max_idxs(i)];
            end
                
        end
        fractionation_nr(j) =counter;
    end

%     
%      figure()
%     hold on
%      plot(times,y)
%     plot(times(local_min_idxs),y(local_min_idxs),'r*')
%     plot(times(local_max_idxs),y(local_max_idxs),'r*')
%     plot(times(max_idx_selected),y(max_idx_selected),'g*')
%     plot(times(min_idx_selected),y(min_idx_selected),'g*')
    
    % figure()
    % hist(fractionation_nr)
    % title('fractionation nr')
    fractionation = fractionation_nr > mean(fractionation_nr)+std(fractionation_nr); %empirical number

    fibrotic_patch = fibrotic_amplitudes + fractionation;
    fibrotic = fibrotic_patch ==2;


    % 
     figure()
     hold all
    % 
     patch('vertices', data.(fields{1}).nodes, 'faces', data.(fields{1}).mesh, 'facecolor', 'interp', 'cdata', int32(fibrotic), 'edgecolor', 'none') ;
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
