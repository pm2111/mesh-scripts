%bar plot fcn

y = [   4.3992,  2.7603,3.5362,3.5134];
figure()
bar(y)
set(gca,'Fontsize',35)
%xticks([1 2 3 4])
set(gca,'XTick',[1 2 3 4])
set(gca,'XTickLabel',{ 'Patient 1' ,'Patient 2' ,'Patient 3' ,'Patient 4' })
ylabel('Mean distance [mm]')

for i=1:4
    means(i) = mean(vertcat(distances{:,i}));
end


figure()
hold on
set(gca,'Fontsize',35)
%xticks([1 2 3 4])
bar(means)
set(gca,'XTick',[1 2 3 4])
set(gca,'XTickLabel',{ 'LV epi' ,'LV endo' ,'RV endo','septum'})
ylabel('Mean distance [mm]')