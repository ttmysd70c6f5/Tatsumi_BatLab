%% Set the directory of lfp data
lfp_dir = 'C:\Users\Tatsumi\Documents\KQTY_NP\FlyingPixels\29968\230531\ephys\20230531_115208_merged.LFP';

%% Export the median of lfp (Run once)
% Extract the median separately to avoid loading whole electrodes at once.
% My computer crashed because of the shortage of RAM storage when I load
% the data for the whole electrodes
ExportLFPMedian(lfp_dir)

%% LFP spectrum
nChan = 384; % Number of channels
lfp_dir = 'C:\Users\Tatsumi\Documents\KQTY_NP\FlyingPixels\29968\230531\ephys\20230531_115208_merged.LFP';
% lfp_dir = 'C:\Users\YartsevPixels1\Desktop\Tatsumi\Data\230531\20230531_115208_merged.LFP';

nChan_page = 12;
nPages = ceil(nChan/nChan_page);
svdir = 'C:\Users\Tatsumi\Documents\KQTY_NP\FlyingPixels\29968\230531\ephys\20230531_115208_merged.LFP\ripples';

% for page = 1:nPages
ch = 1;
for page = 1:nPages
    f = figure;
    f.PaperOrientation = 'landscape';
    f.PaperPosition = [0.05 0.05 0.9 0.9];
    f.PaperUnits = 'normalized';
    set(gcf,'Units','normalized','OuterPosition',[0.05 0.05 0.9 0.9])
    tl = tiledlayout(4,6,'TileSpacing','Compact');

    ch_now = 0;
    while ch_now < nChan_page
        
        RippleSpectrum_waveform_LowMemory(ch,lfp_dir) % Plot Waveform and spectrum

        pause(0.3)

        ch = ch + 1;
        ch_now = ch_now + 1;
    end

%     saveas(f,[svdir,'\ripples_',num2str(page+1000),'.pdf'])
    exportgraphics(gcf,[svdir,'\ripples_',num2str(page+1000),'.pdf'])
    pause(0.3)
    close all
end

