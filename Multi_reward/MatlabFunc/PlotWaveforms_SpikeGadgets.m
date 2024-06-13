function PlotWaveforms_SpikeGadgets(path_to_recording_dir,probeNum,FigDir)
% path_to_recording_dir: the directory containing the subdirectories for
%     kilosort outputs and spikeband data
% probeNum: probe ID
% FigDir: directory to save waveform plots

% waveData = getSpikeWaveforms(path_to_recording_dir,probeNum,'kilosort_outdir_name',sprintf('kilosort_outdir_probe%d',probeNum));
waveData = getSpikeWaveforms(path_to_recording_dir,probeNum);
if isempty(dir(FigDir))
    mkdir(FigDir)

    for i = 1:length(waveData)
    % for i = 1:1
        unit = waveData.waveforms(i).unit; % unit id
        t = waveData.t; % time range
        wv = waveData.waveforms(i).raw_mean; % waveform
        wv_std = waveData.waveforms(i).raw_std; % std
    
        y_max = max(mean(wv,1)) * 2;
        y_min = min(mean(wv,1)) * 2;
    
        % t = wv_window(1):1/30:wv_window(2);
    
        fig = figure('Visible','off');
        boundedline(t,wv,wv_std,'color',[200 200 200]/255,'alpha','transparency', 0.1)
        
        % hold on
        % plot(t,wv,'Color',[150,150,150,10]/255)
        % plot(t,mean(wv,1),'k','LineWidth',1)
        % hold off
        % title(sprintf('Unit %d',unit))
        xlabel('Time (ms)')
        ylabel('Amplitude (uv)')
        ylim([y_min, y_max])
    
        saveas(fig,fullfile(FigDir,sprintf("unit%d.png",unit)))
    end
else
    error('Saving directory is not empty.')
end
