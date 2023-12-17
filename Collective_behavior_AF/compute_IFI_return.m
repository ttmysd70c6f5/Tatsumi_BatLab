function [IFI_ret_xmean, IFI_ret_ymean, IFI_ret_xstd, IFI_ret_ystd] = compute_IFI_return(btakeoff,bat_nms, bat_clr,Fs,figopt,ss,scale)

n_tags = size(btakeoff,2);    

IFI_return = cell(1,n_tags);
IFI_ret_xmean = zeros(1,n_tags);
IFI_ret_ymean = zeros(1,n_tags);
IFI_ret_xstd = zeros(1,n_tags);
IFI_ret_ystd = zeros(1,n_tags);


for i=1:n_tags
    if strcmp(scale,'linear')
        bin_len = figopt.ifibin; % bin duration (timepoints)
        bin_max = figopt.ifimax; % maximum ifi to count for IFI return plot
        IFIedges = 0:bin_len:bin_max;
        
        ifi = diff(find(btakeoff(:,i))) / Fs;
        ifi = ifi(ifi<bin_max);

        ifi_now = ifi(1:end-1);
        ifi_next = ifi(2:end);
        [ifi_cnt,~,~] = histcounts2(ifi_now,ifi_next,IFIedges,IFIedges);
        IFI_return{1,i} = ifi_cnt;

    elseif strcmp(scale,'log')
        bin_len = figopt.logifibin; % bin duration (timepoints)
        bin_max = figopt.logifimax; % maximum ifi to count for IFI return plot
        IFIedges = 0:bin_len:bin_max;
        
        ifi = log10(diff(find(btakeoff(:,i))) / Fs);
        ifi = ifi(ifi<bin_max);

        ifi_now = ifi(1:end-1);
        ifi_next = ifi(2:end);
        [ifi_cnt,~,~] = histcounts2(ifi_now,ifi_next,IFIedges,IFIedges);
        IFI_return{1,i} = log10(ifi_cnt);
    else
        error('Unknown scale specified')
    end
    

    IFI_ret_xmean(1,i) = mean(ifi_now);
    IFI_ret_ymean(1,i) = mean(ifi_next);
    IFI_ret_xstd(1,i) = std(ifi_now);
    IFI_ret_ystd(1,i) = std(ifi_next);
end


%=== FIGURE: IFI return plot
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.05 0.2 0.9 0.35]);
tl = tiledlayout(1, n_tags, 'TileSpacing', 'tight');
title(tl,sprintf('IFI return plot (%s): day%d',scale,ss)); xlabel(tl,'IFI(n) (sec)'); ylabel(tl,'IFI(n+1) (sec)')
for i=1:n_tags
   nexttile(); imagesc(IFI_return{1,i}); colormap jet; set(gca,'YDir','normal')
   title(bat_nms(i,:),'FontSize',12,'Color',bat_clr(i,:))
%    xticks(figopt.ifiedges); 
%    yticks(figopt.ifiedges); 
   colorbar
end


end