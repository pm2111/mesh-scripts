function [fractionation]= fractionation_counter(filename)
data = load(filename);
fields =fieldnames(data);
%QRS fractionation
  [value,Q_ind] = min(abs(data.(fields{1}).epipots.t - data.(fields{1}).surfparam.timestamps.qrsStart));
   [value2,T_ind] = min(abs(data.(fields{1}).epipots.t - data.(fields{1}).surfparam.timestamps.tWaveStart));

%find max - deflection and max amplitude for every node of the ecgi mesh    

 for k=1:size(data.(fields{1}).nodes,1)

     QT_int_phi = data.(fields{1}).epipots.potential(k,Q_ind:T_ind);
     QT_int_t = data.(fields{1}).epipots.t(Q_ind:T_ind);
     max_phi = max(QT_int_phi);
     min_phi = min(QT_int_phi);
     max_peak2peak = max_phi-min_phi;
     dx = QT_int_phi(2:end) -QT_int_phi(1:end-1);
     dt = QT_int_t(2:end) -QT_int_t(1:end-1);
     dQint_dt= dx./dt;
     max_dQ_int_dt_neg = max(-dQint_dt);
     max_neg_gradient(k) = max_dQ_int_dt_neg;
          positive = dQint_dt>0;

          for i=1:size(QT_int_phi,2)-2
        neg_deflection_start(i) = positive(i+1)==0 && positive(i) ==1;
        neg_deflection_end(i) = positive(i+1)==1 && positive(i) ==0;
          end
%      deflection_end = neg_deflection_start-neg_deflection_end; %end is if value is -1
%      single_deflections = neg_deflection_start==1 & neg_deflection_end==1;

     neg_deflection_start_ind = find(neg_deflection_start ==1);
     neg_deflection_end_ind =  find(neg_deflection_end ==1);
     %compute the deflection peak to peak amplitude
     amplitudes = [];
     max_dphi_dt_int = [];
     frac_dummy =0;

     if neg_deflection_end_ind(1) <= neg_deflection_start_ind(1) %check that the 1st deflection start is before 1st deflection end
         neg_deflection_start(1) =1;
     end
     neg_deflection_start_ind = find(neg_deflection_start ==1);
     neg_deflection_end_ind =  find(neg_deflection_end ==1);
     j=1;
     for i = neg_deflection_start_ind(1:size(neg_deflection_end_ind,2))
             amplitudes(j) = ((QT_int_phi(i))-(QT_int_phi(neg_deflection_end_ind(j))));
             max_dphi_dt_int(j) = max(-dQint_dt(i:neg_deflection_end_ind(j)));
          
             j =j+1;
     end
     peak2peak(k) = max(amplitudes);
 end
    
% perform analysis for the regions iin the ECGi dataset
labels = fieldnames(data.(fields{1}).labels);
fractionation =cell(1,size(labels,1));
for z =1:size(labels,1)
    indices = data.(fields{1}).labels.(labels{z,1});
    for l=1:size(indices,1)
        k = indices(l);
        
        QT_int_phi = data.(fields{1}).epipots.potential(k,Q_ind:T_ind);
        QT_int_t = data.(fields{1}).epipots.t(Q_ind:T_ind);
        
        dx = QT_int_phi(2:end) -QT_int_phi(1:end-1);
        dt = QT_int_t(2:end) -QT_int_t(1:end-1);
        dQint_dt= dx./dt;
        max_dQ_int_dt_neg = max(-dQint_dt);
        positive = dQint_dt>0;
        for i=1:size(QT_int_phi,2)-2
            neg_deflection_start(i) = positive(i+1)==0 && positive(i) ==1;
            neg_deflection_end(i+1) = positive(i+1)==1 && positive(i) ==0;
        end
        
        %single_deflections = neg_deflection_start==1 & neg_deflection_end==1;
        
        neg_deflection_start_ind = find(neg_deflection_start ==1);
        neg_deflection_end_ind =  find(neg_deflection_end ==1);
        % single_neg_deflections_end_ind = find(single_deflections ==1);
        %long_deflection_end_ind = find(deflection_end ==-1);
        %compute the deflection peak to peak amplitude
        amplitudes = [];
        max_dphi_dt_int = [];
        frac_dummy =0;
        if neg_deflection_end_ind(1) <= neg_deflection_start_ind(1) %check that the 1st deflection start is before 1st deflection end
            neg_deflection_start(1) =1;
        end
        neg_deflection_start_ind = find(neg_deflection_start ==1);
        neg_deflection_end_ind =  find(neg_deflection_end ==1);
        j=1;
        for i = neg_deflection_start_ind(1:size(neg_deflection_end_ind,2))
            amplitudes(j) = ((QT_int_phi(i))-(QT_int_phi(neg_deflection_end_ind(j))));
            max_dphi_dt_int(j) = max(-dQint_dt(i:neg_deflection_end_ind(j)));
%             stats(j,1) = amplitudes(j) >= 0.1*peak2peak(k);
%             stats(j,2)= max_dphi_dt_int(j) >= 0.7*max_dQ_int_dt_neg;
%             stats(j,3) =  amplitudes(j)>= 0.05*median(peak2peak);
%             stats(j,4) = max_dphi_dt_int(j) >= 0.05*median(max_neg_gradient);
            if amplitudes(j) >= 0.1*peak2peak(k) && max_dphi_dt_int(j) >= 0.7*max_dQ_int_dt_neg && amplitudes(j)>= 0.05*median(peak2peak) && max_dphi_dt_int(j) >= 0.05*median(max_neg_gradient)
                frac_dummy = frac_dummy +1;
            end
            j =j+1;
        end
        fractionation{1,z}(end+1) = frac_dummy;
    end
end
