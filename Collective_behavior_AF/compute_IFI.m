function [IFI_med, IFI_uppq, IFI_lowq] = compute_IFI(btakeoff,bat_nms,bat_clr,Fs,figopt,ss,scale)

n_tags = size(btakeoff,2);

IFI = cell(1,n_tags);
IFI_med = zeros(1,n_tags);
IFI_uppq = zeros(1,n_tags);
IFI_lowq = zeros(1,n_tags);

%=== FIGURE: IFI plot
if strcmp(scale,'linear')
    for i=1:n_tags
        ifi = diff(find(btakeoff(:,i))) / Fs; % intervals (sec)
        IFI{1,i} = ifi;
        IFI_med(1,i) = median(ifi);
        IFI_uppq(1,i) = quantile(ifi,0.75);
        IFI_lowq(1,i) = quantile(ifi,0.25);
    end
    
    figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.05 0.55 0.9 0.35]);
    tl = tiledlayout(1, n_tags, 'TileSpacing', 'tight');
    title(tl, sprintf('Inter-flight interval day%d',ss));
    xlabel(tl,'Time (sec)'); ylabel(tl,'Count')
    hedge = figopt.hedges; ifi_xlim = figopt.ifixlim;
    
    for i=1:n_tags
        nexttile(); histogram(IFI{1,i},hedge);
        set(gca,'yscale','log')
        title(bat_nms(i,:),'FontSize',12,'Color',bat_clr(i,:))
        xlim([0,ifi_xlim])
        hold on;
        plot(repmat(median(IFI{1,i}),[1,2]),ylim,'k--');
        hold off
        legend('',sprintf('Median: %.0fsec',median(IFI{1,i})));
    end
    
elseif strcmp(scale,'log')
    for i=1:n_tags
        ifi = diff(find(btakeoff(:,i))) / Fs; % intervals (sec)
        
        IFI_med(1,i) = median(ifi);
        IFI_uppq(1,i) = quantile(ifi,0.75);
        IFI_lowq(1,i) = quantile(ifi,0.25);
        
        IFI{1,i} = log10(ifi);
    end
    
    figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.05 0.55 0.9 0.35]);
    tl = tiledlayout(1, n_tags, 'TileSpacing', 'tight');
    title(tl, sprintf('Inter-flight interval day%d',ss));
    xlabel(tl,'Log10 time (sec)'); ylabel(tl,'Count')
    hedge = figopt.loghedges; ifi_xlim = figopt.logifixlim;
    
    for i=1:n_tags
        nexttile(); histogram(IFI{1,i},hedge);
        set(gca,'yscale','log')
        title(bat_nms(i,:),'FontSize',12,'Color',bat_clr(i,:))
        xlim([0,ifi_xlim])
        hold on;
        plot(repmat(median(IFI{1,i}),[1,2]),ylim,'k--');
        hold off
        legend('',sprintf('Median: %.0fsec',10^median(IFI{1,i})));
    end
else
   error('Unknown scale specified')
end
end