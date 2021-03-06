function ECGs = plot_ecgs(signal,leg,individual_figures)
%PLOT_ECGS_20141124 Summary of this function goes here
%   More efficient use of vertical space

    maxs = zeros(1,12);
    mins = zeros(1,12);
    row_max = zeros(1,3);
    row_min = zeros(1,3);
    row_height = zeros(1,3);

    % All on one figure, or one figure/lead
    % individual_figures = false;
    if(~individual_figures)
        figure
%         figure(1);
    end

    key_out = {'I'; 'II'; 'III'; 'aVR'; 'aVL'; 'aVF'; 'V1'; 'V2'; 'V3'; 'V4'; 'V5'; 'V6'};
    layout = fliplr(reshape(1:12,3,4)');
    figure('Position',[0 0 500 1000]); %left bottom width height
    colours = {'k--','r','g','b','r--','g--','b--','r:'};

    handles = [];

    % Grid: Horizontal: large unit = 200 ms, small unit = 40 ms
    %       Vertical:   large unit = 0.5 mV, small unit = 0.1 mV
    x_large_div = 200;            % 200 ms
    y_large_div = 0.5;            % 0.5 mV
    x_small_div = x_large_div;  % 40 ms
    y_small_div = y_large_div;  % 0.1 mV

    % First go through all the electrodes and generate the ECG leads
    % Also record min's and max's for each row
    ECGs={};
    for simu = 1:12
% 
%         for n=1:10
%             %disp(sprintf('%s=signal{simu}(:,%i);',key_in(n,:),n));
%            eval(sprintf('%s=signal{simu}(:,%i);',key_in(n,:),n)); % Nasty!
%         end
% 
% %         ECG = zeros(size(signal{simu},1), 12);
% %         ECG(:,1) = LA - RA; % I
% %         ECG(:,2) = LL - RA; % II
% %         ECG(:,3) = LL - LA; % III
% %         ECG(:,4) = RA - 0.5*(LA+LL); % aVR
% %         ECG(:,5) = LA - 0.5*(RA+LL); % aVL
% %         ECG(:,6) = LL - 0.5*(RA+LA); % aVF
% %         % Reference
% %         Vr = (RA + LA + LL)/3;
% %         ECG(:,7:12) = signal{simu}(:,5:10) - repmat(Vr, 1, 6);
% 
%         % Output
         ECGs = signal;
        ECG=signal{simu}';
%         % Update max and min values
      these_maxs = max(ECG);
     C = bsxfun(@gt,these_maxs,maxs);
       maxs(C) = these_maxs(1);
       these_mins = min(ECG);
       C = bsxfun(@lt,these_mins,mins);
         mins(C) = these_mins(1);
     end
    % For each row in the plot, find the largest (abs) max/min, round up, and
    % save for subplot later. (We know that 1:4 are bottom row, 5:8 middle
    % row...)
    for r=1:3
        inds = layout(4*(r-1)+1:4*r); % The indices of the frames of this row in maxs and mins
        row_max(r) = ceil(max(maxs(inds))/y_large_div)*y_large_div;
        row_min(r) = floor(min(mins(inds))/y_large_div)*y_large_div;
        row_height(r) = row_max(r)-row_min(r);
    end

    % Refactor row_lim into proportions
    row_height_rel = row_height ./ sum(row_height);
    % Cumulative sum
    row_cum_sum = [0 cumsum(row_height_rel)];

    % Done finding limits, now we can plot
    for simu = 1

        for n=1:12
            m = n;
            if individual_figures
                if simu == 1
                    handles(n)=figure('Position',[0 0 500 200]);
                else
                    figure(handles(n))
                end
                subaxis(1,1,1,'Margin',0.001)
            else
                x = mod(n-1,4); % 0 based
                y = ceil(n/4);  % 1 based!
                subplot('Position',[x/4 row_cum_sum(y) 1/4 row_height_rel(y)]);

                % Using the same figure, so grab the handle now, just once
                if isempty(handles)
                    handles = gcf;
                end
            end
            hold on;

            if ischar(colours{simu})
                plot(ECGs{n},colours{simu},'LineWidth',2);
            else
                plot(ECGs{n},'LineWidth',2,'color',colours{simu});
            end
            if simu==1 % First time through
                grid on
                set(gca,'GridLineStyle','-','LineWidth',0.1);
                set(gca, 'XColor', 'r')
                set(gca, 'YColor', 'r')
                grid minor
                set(gca,'MinorGridLineStyle','-');
                if(~individual_figures)
                    font_size = 0.3*y_large_div/row_height(y)*4;
                    text(40, row_max(y)-0.25, key_out(m,:), 'FontName', 'Times', 'FontUnits', 'normalized', 'FontSize', font_size);
                else
                    font_size = 0.1*4;
                    text(60, 0.1, key_out(m,:), 'FontName', 'Times', 'FontUnits', 'normalized', 'FontSize', font_size);
                end                
            end

            x_lim = xlim;
            if x_lim(2)<400
                x_lim(2)=400;
            end

            xlimceil = ceil(x_lim(2)/x_large_div)*x_large_div; % Round up in x
    %             xlimceil = 400; % override
            xlim([0 xlimceil]);
            set(gca,'XTick',0:x_small_div:xlimceil);
            set(gca,'XTickLabel',[]);

            if(~individual_figures)
                ylim([row_min(y) row_max(y)]);
                set(gca,'YTick',row_min(y):y_small_div:row_max(y));
            end
            set(gca,'YTickLabel',[]);

    %         pause
        end
    %     pause
    end

    if simu>1 % no point having a legend for one trace
        for h=handles
            figure(h)
            legend(leg,'Location','Best', 'FontSize', 12)
        end
    end

    % To get the aspect ratio square...
    max_dim = 800; % px
    % Total number of divisions in the figure
    y_fig_divs_total = sum(row_height)/y_small_div;
    x_fig_divs_total = 4*xlimceil/x_small_div;

    if x_fig_divs_total>y_fig_divs_total
        xfig = max_dim;
        yfig = xfig*y_fig_divs_total/x_fig_divs_total;
    else
        yfig = max_dim;
        xfig = yfig*x_fig_divs_total/y_fig_divs_total;
    end

    % Now resize the figure(s)
    if ~individual_figures % just for the big one for now
        set(handles,'Position',[0 0 xfig yfig]);
        refresh(handles)
    end

end
