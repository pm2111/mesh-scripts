function ECGs = plot_ecgs(signal,leg,individual_figures,col,ts)
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
       % figure
%         figure(1);
    end

    key_out = {'I'; 'II'; 'III'; 'aVR'; 'aVL'; 'aVF'; 'V1'; 'V2'; 'V3'; 'V4'; 'V5'; 'V6'};
    layout = fliplr(reshape(1:12,6,2))';
    %figure('Position',[0 0 500 1000]); %left bottom width height
  % colours = {'b-.','r--','g--','k:','r--','g--','b:','r:'};
  colours={'b','k--', 'g--','r--','r','m','k','y','c'};

    handles = [];

    % Grid: Horizontal: large unit = 200 ms, small unit = 40 ms
    %       Vertical:   large unit = 0.5 mV, small unit = 0.1 mV
    x_large_div = 100;            % 100 ms
    y_large_div = 1;            % 0.5 mV
    x_small_div = x_large_div/2;  % 50 ms
    y_small_div = y_large_div/2;  % 0.1 mV

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
       maxs =10*ones(1,12);
       these_mins = min(ECG);
       C = bsxfun(@lt,these_mins,mins);
         mins = -13*ones(1,12);
     end
    % For each row in the plot, find the largest (abs) max/min, round up, and
    % save for subplot later. (We know that 1:4 are bottom row, 5:8 middle
    % row...)
    for r=1
        inds = layout(6*(r-1)+1:6*r); % The indices of the frames of this row in maxs and mins
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
                   % handles(n)=figure('Position',[0 0 500 200]);
                else
                   % figure(handles(n))
                end
                subaxis(1,1,1,'Margin',0.001)
            else
                x = mod(n,13)-1; % 0 based
                y = ceil(n/13);  % 1 based!
                subplot('Position',[-0.5+x/12+n/500 row_cum_sum(y) 1/12 row_height_rel(y)]);

                % Using the same figure, so grab the handle now, just once
                if isempty(handles)
                    handles = gcf;
                end
            end
            hold on;

            if ischar(colours{simu})
%                 if n>0
%                 plot(ts,ECGs{n},'color',colours{col},'LineWidth',3);
%                 
%                 elseif n==8
%                    plot(ts(1:end-125+10*n),ECGs{n}(1:end-125+10*n),'color',colours{col},'LineWidth',3);
% 
%                 else
%                 plot(ts(1:end-135+10*n),ECGs{n}(1:end-135+10*n),'color',colours{col},'LineWidth',3);
%                 end
%             else
                plot(ts,ECGs{n},colours{col},'LineWidth',3);
            end
                grid on
                set(gca,'GridLineStyle','-','LineWidth',0.1);
                set(gca, 'XColor', 'r')
                set(gca, 'YColor', 'r')
                grid minor
                set(gca,'MinorGridLineStyle','-');
                if(~individual_figures)
                    font_size = 25;
                    text(40, maxs(1)-0.75, key_out(m,:), 'FontSize', font_size);
                else
                    font_size = 0.1*7;
                    text(60, 0.01, key_out(m,:), 'FontName', 'Times', 'FontUnits', 'normalized', 'FontSize', font_size);
                end                

            x_lim = xlim;
            if x_lim(2)<400
                x_lim(2)=200;
            end

            xlimceil = ceil(x_lim(2)/x_large_div)*x_large_div; % Round up in x
    %             xlimceil = 400; % override
            xlim([0 xlimceil]);
            set(gca,'XTick',0:x_small_div:xlimceil);
            set(gca,'XTickLabel',[]);

            if(~individual_figures)
                ylim([mins(1) maxs(1)]);
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
           % legend(leg,'Location','Best', 'FontSize', 12)
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
