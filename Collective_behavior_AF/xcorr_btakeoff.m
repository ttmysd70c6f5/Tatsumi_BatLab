function xcorr_btakeoff(btakeoff,n_tags,t,bat_nms,bat_clr,figopt)


% options
maxlag = figopt.cc_maxlag; % maximum time to compute cc (sec)
mask = figopt.cc_mask;
bin_dur = figopt.cc_bin; % bin duration (sec)
bin_edges = 0:bin_dur:t(end);

% Compute xcorrs
xcorr_raw = zeros(2*maxlag/bin_dur+1,n_tags^2);

for i=1:n_tags^2
    id_col = i-n_tags*floor((i-1)/n_tags); % index of column of correlogram
    id_row = ceil(i/n_tags); % index of row of correlogram
    if id_col > id_row
        xcorr_raw(:,i) = zeros(size(xcorr_raw,1),1);
    
    else
        X = btakeoff(:,id_row); Y = btakeoff(:,id_col); % takeoffs of pairwise bats to compute xcorr
        N = zeros(length(bin_edges)-1,2);
        X_takeoff = t(X); Y_takeoff = t(Y);
        [N(:,1),~] = histcounts(X_takeoff,bin_edges);
        [N(:,2),~] = histcounts(Y_takeoff,bin_edges);
        [xcorr_raw(:,i), lags] = xcorr(N(:,1),N(:,2),maxlag/bin_dur);
    end
end

%=== FIGURE: Trajectories individual bats
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
tl=tiledlayout(n_tags, n_tags, 'TileSpacing', 'tight');
ylabel(tl,'Count'); xlabel(tl,'Time (sec)')
for i=1:n_tags^2
    nexttile(); 
    plot(bin_dur * lags,conv(xcorr_raw(:,i),mask,'same'),'k-'); hold on
    plot([0 0],ylim,'k--'); hold off
    if mod(i,n_tags)==1
        ylabel(bat_nms((i-1)/n_tags+1,:),'FontSize',12,'Color',bat_clr((i-1)/n_tags+1,:))
    end
    if i<=n_tags
        title(bat_nms(i,:),'FontSize',12,'Color',bat_clr(i,:))
    end
    if i<=n_tags*(n_tags-1)
       xticks([]);
    end
end

end