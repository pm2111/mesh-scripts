cd /data/ecgi/
filenames = dir('/data/ChrisReconstructions/baselines/');
for i=1:size(filenames,1)-2
    
    plot_ecgi_rest(strcat(filenames(i+2).folder,'/',filenames(i+2).name))
    saveas(gcf,strcat(filenames(i+2).name,' activation time map lv view'),'png')
end