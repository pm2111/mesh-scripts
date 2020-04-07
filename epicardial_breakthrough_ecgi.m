%epicardial breakthrough times

filenames = dir('/data/ChrisReconstructions/baselines/');

for i=1:size(filenames,1)-2
    data = load(strcat(filenames(i+2).folder,'/',filenames(i+2).name));
    fields = fieldnames(data);
   
    a =  data.(fields{1}).nodeparam.timestamps.activationTimes-data.(fields{1}).nodeparam.timestamps.activationStart;
    while min(a)<5
        [mn,idx ] = min(a);
        a(idx) =[];
    end
        
        
        activation_times(i) =  min(a);

    
        
end