%% Function for the analysis of Extracellular Field Recordings in the Hippocampus and the detection of Ripples

confirmed_by_inspection = zeros(4,1);

for CSC_ch = 0:4:15
    
    if ~exist('CSC_ch'),CSC_ch = input('CSC channel:');end
    disp(['Analyzing CSC',num2str(CSC_ch), ' ...']);
    
    %=== Load Files
    CSC_file = dir(fullfile(cd, ['*CSC' num2str(CSC_ch),'.mat']));load(CSC_file.name);      % CSC File
    load('TTL_timestamps.mat');                                                             % Reference Timestamps for alignment
    load('imp_bat.mat');                                                                    % Tag of the implanted bat
    BHV_file = dir(fullfile(fileparts(CSC_file.folder),'Ext_Beh*','Extracted_Behavior*'));  % Behavioral variables (velocity, acceleration, etc...)
    
    if ~isempty(BHV_file)
        load([BHV_file.folder,'\',BHV_file.name],'r','v_abs','a_abs','t','bflying','batdate','n_tags','T');
    else
        warning('-----No behavioral file detected-------');
        a_abs = randn(1e3,5);
        v_abs = randn(1e3,5);
        batdate = CSC_file.name(7:14);
        t = 1:1e3;
        bflying = zeros(1000,5);
    end
    load('bat_SN');
    
    %=== Options
    CSC_options.lookATpre = 0;
    CSC_options.lookATfly = 0;
    CSC_options.lookATpst = 0;
    CSC_options.lookATsorted = 0;
    CSC_options.checkAlgn = 1;
    CSC_options.extractRipples = 1;
    CSC_options.loadBehavior = 1;
    CSC_options.lookAccel = 0;
    CSC_options.savedata = 1;
    CSC_options.savefigures = 0;  fig_count = 1;
    
    %=== Create analysis folder for storing the results
    if CSC_options.savedata
        analysis_directory=fullfile(pwd,['CSC_Analysis_',batdate,]);
        if ~exist(analysis_directory,'dir')
            mkdir(analysis_directory);
        end
    end
    
    %% Assign relevant variables and define filters
    
    %=== Variables and Params
    TT = ceil((CSC_ch+1)/4);                                                    % Tetrode
    Fs = Estimated_channelFS_Transceiver(1);                                    % Voltage Sampling Frequency (Hz)
    raw_V = double(AD_count_int16*AD_count_to_uV_factor);                       % Raw Voltage Trace
    tstamps = Timestamps_of_first_samples_usec(1)+[0:length(raw_V)-1]/Fs*1e6;   % Timestamps (us)
    
    RP_F_band = [100 200];                  % Ripples frequency band (From MY: [80 160])
    RP_th = 5;                              % Threshold on zscored Ripple power (3.5)
    dec = 10;                               % Downsampling factor for Ripple detection
    
    %=== Ad HOC corrections
    if strcmp(CSC_file.name(1:18),'00000_20230322_CSC') && bat_SN==14407
        raw_V(1,214082000:end) = randn(size(raw_V(1,214082000:end)))*1e2;   % Remove large artifacts toward the end
    end
    
    %% Load single units
    
    if CSC_options.lookATsorted
        TTsFile = dir('SingleUnits*/*SingleUnits*');  load([TTsFile.folder '/' TTsFile.name]);
        n_cells = size(TT_unit,2);
        s = cell(n_cells,1);
        for nc = 1:n_cells
            s{nc} = TT_unit(nc).Timestamps;
            s{nc} = s{nc}(s{nc}>TTL_timestamps.fly(1) & s{nc}<TTL_timestamps.fly(end)); s{nc} = (s{nc}-TTL_timestamps.fly(1))'./1e6;
        end
    else
        s = [];
    end
    
    %% Look at pre-sleep session
    
    if CSC_options.lookATpre
        [~,sample_strt] = min(abs(tstamps-TTL_timestamps.pre(1)));
        [~,sample_stop] = min(abs(tstamps-TTL_timestamps.pre(end)));
        raw_V_pre = raw_V(:,sample_strt:sample_stop);
        tstamps_pre = tstamps(:,sample_strt:sample_stop);   tstamps_pre = (tstamps_pre-tstamps_pre(1))/1e6;
        Look_at_Ripples(raw_V_pre,tstamps_pre,Fs,[],[],[],s);
        close gcf
    end
    
    %% Look at fly session
    
    if CSC_options.lookATfly
        [~,sample_strt] = min(abs(tstamps-TTL_timestamps.fly(1)));
        [~,sample_stop] = min(abs(tstamps-TTL_timestamps.fly(end)));
        raw_V_fly = raw_V(:,sample_strt:sample_stop);
        tstamps_fly = tstamps(:,sample_strt:sample_stop);   tstamps_fly = (tstamps_fly-tstamps_fly(1))/1e6;
        Look_at_Ripples(raw_V_fly,tstamps_fly,Fs,t,v_abs(:,imp_bat),a_abs(:,imp_bat),s);
        close gcf
    end
    
    %% Look at post-sleep session
    
    if CSC_options.lookATpst
        [~,sample_strt] = min(abs(tstamps-TTL_timestamps.pst(1)));
        [~,sample_stop] = min(abs(tstamps-TTL_timestamps.pst(end)));
        raw_V_pst = raw_V(:,sample_strt:sample_stop);
        tstamps_pst = tstamps(:,sample_strt:sample_stop);   tstamps_pst = (tstamps_pst-tstamps_pst(1))/1e6;
        Look_at_Ripples(raw_V_pst,tstamps_pst,Fs,[],[],[],s);
        close gcf
    end
    
    %% Check alignment with behavior by looking at wingbeats artifacts
    
    if CSC_options.checkAlgn
        %=== Find samples of start - stop fly session
        [~,sample_strt] = min(abs(tstamps-TTL_timestamps.fly(1)));
        [~,sample_stop] = min(abs(tstamps-TTL_timestamps.fly(end)));
        
        %=== Cut data and reference timestamps (in s) to first TTL
        raw_V_fly = raw_V(:,sample_strt:sample_stop);
        tstamps_fly = tstamps(:,sample_strt:sample_stop);
        tstamps_fly = (tstamps_fly-tstamps_fly(1))/1e6;
        
        %=== Interpolate voltage signal at behavior time points, filter in the wingbeat band
        dws_V = interp1(tstamps_fly,raw_V_fly,t,'linear','extrap');
        V = zscore(bandpass(dws_V,[5 15],100));
        F = zscore(a_abs(:,imp_bat));
        [c,lags] = xcorr(V,F,10*100,'normalized');
        
        %=== Plot traces and cross-correlation
        figure('units','normalized','outerposition',[0.2 0.2 0.4 0.5]);
        tiledlayout(2,1,'TileSpacing','tight');
        nexttile;   plot(t,zscore(dws_V),t,zscore(v_abs(:,imp_bat)));       xlabel('Time (s)'); ylabel('z score');  legend('Voltage','Velocity');
        nexttile;   plot(lags/100,c);   xlabel('Time lag (s)'); ylabel('Wing beats VS voltage corr');
        hold on;    rectangle('Position',[-0.1 min(c) 0.2 range(c)],'FaceColor',[0 0 0 0.1],'EdgeColor','none');    sgtitle(['TT ',num2str(TT),', CSC ',num2str(CSC_ch)]);
        fig_count = saveFig(analysis_directory,[batdate, '_CSC',  num2str(CSC_ch)],fig_count,CSC_options.savefigures);
    end
    
    %% Ripple detection
    
    %=== Downsample raw voltage and timestamps
    ds_Fs = Fs/dec;                                                     % Downsampled Frequency
    ds_V = downsample(raw_V,dec);                                       % Downsampled Voltage
    ds_tstamps = downsample(tstamps,dec);                               % Downsampled TimeStamps (s)
    
    %=== Define filter for MU
    SP_F_band = [600 6000];                 % Spikes frequency band
    passband_norm = SP_F_band./(0.5*Fs);    % Normalize over sampling frequency
    [b,a] = ellip(4,0.1,70,passband_norm);  % Define filter transfer function coefficients
    
    %=== Ripple Detection (amplitude,time,duration)
    ripple_signal = bandpass(ds_V,RP_F_band,Fs/dec);                    % Bandpass signal at the Ripple Band
    ripple_power  = zscore(abs(hilbert(ripple_signal)));                % Calculate z-scored power as the magnitude of the Hilbert transform
    ripple_power  = smoothdata(ripple_power,'gaussian',0.05*Fs/dec);    % Smooth with a 50 ms Gaussian kernel
    [RPL_a,RPL_t,RPL_d] = findpeaks(ripple_power,Fs/dec,...             % Detect peaks in the ripple power...
        'MinPeakHeight',RP_th,...           % larger than RP_th...
        'MinPeakDistance',0.1);             % with minimal separation of 100 ms
    
    %=== Remove all the ripples happening during flight (+-0.5s)
    extended_flight = movmax(bflying(:,imp_bat),[50 50]);
    r2k = true(size(RPL_t));
    parfor i=1:numel(RPL_t)
        if RPL_t(i)>(TTL_timestamps.fly(1)-Timestamps_of_first_samples_usec(1))/1e6 && RPL_t(i)<(TTL_timestamps.fly(end)-Timestamps_of_first_samples_usec(1))/1e6
            [~,idx] = min(abs((RPL_t(i)-(TTL_timestamps.fly(1)-Timestamps_of_first_samples_usec(1))/1e6)-t));
            if extended_flight(idx,1),r2k(i) = 0;end
        end
    end
    RPL_a = RPL_a(r2k); RPL_t = RPL_t(r2k); RPL_d = RPL_d(r2k);
    
    %=== Calculate range of the raw voltage (+-50ms around ripple peak) or maximum derivative
    RPL_r = zeros(size(RPL_t));
    for i=1:numel(RPL_t)
        
        sample_ripple = round(RPL_t(i)*Fs);
        if sample_ripple-0.05*Fs>0 && sample_ripple+0.05*Fs<numel(raw_V)
            %RPL_r(1,i) = range(raw_V(round(sample_ripple-0.05*Fs):round(sample_ripple+0.05*Fs)));
            RPL_r(1,i) = max(abs(diff(raw_V(round(sample_ripple-0.05*Fs):round(sample_ripple+0.05*Fs)))))*Fs;
            RPL_r(1,i) = RPL_r(1,i)/1e6;    %uV/us
        end
    end
    
    %=== Remove some artifacts
    RPL_t=RPL_t(RPL_r<5);   RPL_a=RPL_a(RPL_r<5); RPL_d=RPL_d(RPL_r<5); RPL_r=RPL_r(RPL_r<5);

    %=== Plot a few examples (+-100ms around ripple peak)
    [rs_RPL_t,rs_idx] = datasample(RPL_t,50);              % Select random subset
    rs_RPL_a = RPL_a(rs_idx);
    rs_RPL_d = RPL_d(rs_idx);
    rs_RPL_r = RPL_r(rs_idx);
    [~,st_idx] = sort(rs_RPL_r);                            % Sort according to...
    rs_RPL_t = rs_RPL_t(st_idx);
    rs_RPL_a = rs_RPL_a(st_idx);
    rs_RPL_d = rs_RPL_d(st_idx);
    rs_RPL_r = rs_RPL_r(st_idx);
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(5,20,'TileSpacing','tight');
    for i =1:50
        sample_ripple = round(rs_RPL_t(i)*Fs);
        if round(sample_ripple-0.1*Fs)<0 || round(sample_ripple+0.1*Fs)>numel(raw_V),continue;end
        nexttile; 
        plot((round(-0.1*Fs):round(+0.1*Fs))/Fs,normalize(raw_V(round(sample_ripple-0.1*Fs):round(sample_ripple+0.1*Fs)),'range'),'k');
        hold on;    plot((round(-0.1*Fs):round(+0.1*Fs))/Fs,-1+normalize(filtfilt(b,a,raw_V(round(sample_ripple-0.1*Fs):round(sample_ripple+0.1*Fs))),'range')); hold off;
        xlim([-0.1 0.1]);   xlabel('s');    ylabel('uV');   %xticks([]); yticks([]);
        title(['A = ', num2str(rs_RPL_a(i),2), ' D = ' num2str(rs_RPL_d(i)*1e3,2),' ms, A = ', num2str(rs_RPL_r(i),2), ' uV/us']);
        nexttile;
        spectrogram_AF_v0(raw_V(round(sample_ripple-0.5*Fs):round(sample_ripple+0.5*Fs)),Fs,0.2,0.05,[70 200]);
        %plot(abs(diff(raw_V(round(sample_ripple-0.1*Fs):round(sample_ripple+0.1*Fs)))),'k');
        %     [cfs,f] = cwt(raw_V(round(sample_ripple-0.5*Fs):round(sample_ripple+0.5*Fs)),Fs);
        %     imagesc((round(-0.1*Fs):round(+0.1*Fs))/Fs,f(52:80),imgaussfilt(abs(cfs(52:80,:)),[2 500]));
        %     xlabel('Time (s)');ylabel('Frequency (Hz)');axis xy;
    end
    sgtitle(['TT ',num2str(TT),', CSC ',num2str(CSC_ch)]);
    fig_count = saveFig(analysis_directory,[batdate, '_CSC',  num2str(CSC_ch)],fig_count,CSC_options.savefigures);
    
%     %=== Plot whole session trace
%     figure('units','normalized','outerposition',[0 0.2 1 0.4]);
%     tiledlayout(2,1,'TileSpacing','tight');
%     ax(1) = nexttile;   plot((ds_tstamps-ds_tstamps(1))/60e6,ds_V); ylabel('Voltage (uV)');
%     text((TTL_timestamps.fly(1)-ds_tstamps(1))/60e6,2e3,'Start Fly Session');
%     text((TTL_timestamps.fly(end)-ds_tstamps(1))/60e6,2e3,'Stop Fly Session');
%     ax(2) = nexttile;   plot((ds_tstamps-ds_tstamps(1))/60e6,ripple_power); ylabel('Ripple Power (sd)');
%     xlabel('Minutes');  linkaxes(ax,'x');   sgtitle(['TT ',num2str(TT),', CSC ',num2str(CSC_ch)]);
%     fig_count = saveFig(analysis_directory,[batdate, '_CSC',  num2str(CSC_ch)],fig_count,CSC_options.savefigures);
    
    %% Plot spectrum associated with detected Ripples
    
    n=512; P=zeros(numel(RPL_t),n);
    for i = 1:numel(RPL_t)
        sample_ripple = round(RPL_t(i)*Fs/dec);
        if round(sample_ripple-0.5*Fs/dec)<0 || round(sample_ripple+0.5*Fs/dec)>numel(ds_V) || RPL_r(i)>5,continue;end  % Do not include artifacts and ripples at the edges of the session
        ripple_chunk = ds_V(sample_ripple-n/2:sample_ripple+n/2-1);
        Y = fft(ripple_chunk);
        P(i,:) = abs(Y/n).^2;
    end
    P_ave = mean(P,1);  f_ave = Fs/dec*(0:n/2)/n;
    % Y_all = fft(ds_V);
    % n_all = length(ds_V);
    % P_all = abs(Y_all/n_all).^2;
    % f_all = Fs/dec*(0:n_all/2)/n_all;
    
    figure;
    plot(f_ave,P_ave(1:n/2+1));
    xlim([50 200]);    set(gca, 'YScale', 'log');
    xlabel('Frequency (Hz)');   h = gca;    h.YAxis.Visible = 'off';
    sgtitle(['TT ',num2str(TT),', CSC ',num2str(CSC_ch)]);
    fig_count = saveFig(analysis_directory,[batdate, '_CSC',  num2str(CSC_ch)],fig_count,CSC_options.savefigures);
    
    %% Ask for Visual Inspection
    
    answer = questdlg('Are SWRs present on this Tetrode?');
    if strcmp(answer,'Yes')
        confirmed_by_inspection(TT) = 1;
    end
    
    %% Save detected Ripples, referenced to first TTL of the corresponding session
    
    idx_pre = find(RPL_t>(TTL_timestamps.pre(1)-Timestamps_of_first_samples_usec(1))/1e6 & RPL_t<(TTL_timestamps.pre(end)-Timestamps_of_first_samples_usec(1))/1e6);
    idx_fly = find(RPL_t>(TTL_timestamps.fly(1)-Timestamps_of_first_samples_usec(1))/1e6 & RPL_t<(TTL_timestamps.fly(end)-Timestamps_of_first_samples_usec(1))/1e6);
    idx_pst = find(RPL_t>(TTL_timestamps.pst(1)-Timestamps_of_first_samples_usec(1))/1e6 & RPL_t<(TTL_timestamps.pst(end)-Timestamps_of_first_samples_usec(1))/1e6);
    RPL_pre.t = RPL_t(idx_pre)'-(TTL_timestamps.pre(1)-Timestamps_of_first_samples_usec(1))/1e6;    RPL_pre.a = RPL_a(idx_pre)';    RPL_pre.d = RPL_d(idx_pre)';    RPL_pre.r = RPL_r(idx_pre)';
    RPL_fly.t = RPL_t(idx_fly)'-(TTL_timestamps.fly(1)-Timestamps_of_first_samples_usec(1))/1e6;    RPL_fly.a = RPL_a(idx_fly)';    RPL_fly.d = RPL_d(idx_fly)';    RPL_fly.r = RPL_r(idx_fly)';
    RPL_pst.t = RPL_t(idx_pst)'-(TTL_timestamps.pst(1)-Timestamps_of_first_samples_usec(1))/1e6;    RPL_pst.a = RPL_a(idx_pst)';    RPL_pst.d = RPL_d(idx_pst)';    RPL_pst.r = RPL_r(idx_pst)';
    save([analysis_directory,'\Ripples',CSC_file.name(6:end-4),'_AF.mat'],'RPL_pre','RPL_fly','RPL_pst','P_ave');
    
    %% Analysis of LFP and accelerometer signal
    
    if CSC_options.lookAccel
        
        %=== Extract Fly session
        [~,sample_strt] = min(abs(tstamps-TTL_timestamps.fly(1)));
        [~,sample_stop] = min(abs(tstamps-TTL_timestamps.fly(end)));
        raw_V_fly = raw_V(:,sample_strt:sample_stop);
        tstamps_fly = tstamps(:,sample_strt:sample_stop);   tstamps_fly = (tstamps_fly-tstamps_fly(1))/1e6;
        
        %=== Resample Voltage, Flight, Acceleration and Time at Fs/dec Hz
        us_V = downsample(raw_V_fly,dec);
        us_t = downsample(tstamps_fly,dec);
        us_b = interp1(t,bflying(:,imp_bat),us_t,'nearest','extrap');
        us_a = interp1(t,a_abs(:,imp_bat),us_t,'nearest','extrap');
        
        %=== Ripple Epochs
        us_r = zeros(size(us_a));
        us_r(round(RPL_fly.t*Fs/dec))=1;
        us_r = movmax(us_r,[0.1*Fs/dec 0.1*Fs/dec]);
        
        %=== Define Extended flight epochs (+-1s around flight)
        us_e = movmax(us_b,[1*Fs/dec 1*Fs/dec]);
        
        %=== Define distances from implanted bat
        imp_bat_dist = zeros(T,n_tags-1);
        oth_bat = [1:imp_bat-1,imp_bat+1:n_tags];
        c=1;
        for i = oth_bat
            imp_bat_dist(:,c) = vecnorm(r(:,:,imp_bat)-r(:,:,i),2,2);
            c=c+1;
        end
        us_d = interp1(t,min(imp_bat_dist,[],2),us_t,'nearest','extrap');
        
        %=== Define General Movement Epochs and plot
        us_m = (us_a-smoothdata(us_a,'movmean',1*Fs/dec)>5e-3 | us_a-smoothdata(us_a,'movmean',1*Fs/dec)<-5e-3) |  us_e;
        figure; plot(us_t,zscore(us_a));    hold on;    area(us_t,-7*us_e,'FaceAlpha',0.3,'LineStyle','none'); area(us_t,5*us_m,'FaceAlpha',0.3,'LineStyle','none');    hold off;
        
        %=== Look at FFT of the voltage signal during flight VS non-flight
        % (forcing to NaN and then interpolating gives similar spectra)
        figure('units','normalized','outerposition',[0.7 0 0.15 1]);
        tiledlayout(4,1,'TileSpacing','none');
        trace = us_V(find(us_m));       % Movement epochs
        n = 2^nextpow2(numel(trace));   Y = fft(trace,n);   P = abs(Y/n).^2;    f = Fs/dec*(0:n/2)/n;
        PSD = P(1:n/2+1);   sm_PSD = smoothdata(PSD,'movmedian',n/(Fs/dec));
        ax(1) = nexttile;   plot(f,mag2db(sm_PSD),'.g','LineWidth',2);     grid on;    xlabel('Hz');   legend('Movement');   set(gca, 'fontsize',16);
        trace = us_V(find(us_e));       % Flight epochs
        n = 2^nextpow2(numel(trace));   Y = fft(trace,n);   P = abs(Y/n).^2;    f = Fs/dec*(0:n/2)/n;
        PSD = P(1:n/2+1);   sm_PSD = smoothdata(PSD,'movmedian',n/(Fs/dec));
        ax(2) = nexttile;   plot(f,mag2db(sm_PSD),'.k','LineWidth',2);     grid on;    xlabel('Hz');   legend('Flight');     set(gca, 'fontsize',16);
        trace = us_V(find(~us_e));      % Non Flight epochs
        n = 2^nextpow2(numel(trace));   Y = fft(trace,n);   P = abs(Y/n).^2;    f = Fs/dec*(0:n/2)/n;
        PSD = P(1:n/2+1);   sm_PSD = smoothdata(PSD,'movmedian',n/(Fs/dec));
        ax(3) = nexttile;   plot(f,mag2db(sm_PSD),'.r','LineWidth',2);     grid on;    xlabel('Hz');   legend('Non-Flight');     set(gca, 'fontsize',16);
        trace = us_V(find(~us_m));      % Still epochs
        n = 2^nextpow2(numel(trace));   Y = fft(trace,n);   P = abs(Y/n).^2;    f = Fs/dec*(0:n/2)/n;
        PSD = P(1:n/2+1);   sm_PSD = smoothdata(PSD,'movmedian',n/(Fs/dec));
        ax(4) = nexttile;   plot(f,mag2db(sm_PSD),'.b','LineWidth',2);     grid on;    xlabel('Hz');   legend('Still');  set(gca, 'fontsize',16);
        linkaxes(ax,'x');   xlim([1 200]);
        % trace = us_V(find(us_r));      % Ripple epochs
        % n = 2^nextpow2(numel(trace));   Y = fft(trace,n);   P = abs(Y/n).^2;    f = Fs/dec*(0:n/2)/n;
        % PSD = P(1:n/2+1);   sm_PSD = smoothdata(PSD,'movmedian',n/(Fs/dec));
        % ax(5) = nexttile;   plot(f,mag2db(sm_PSD),'.c','LineWidth',2);     grid on;    xlabel('Hz');   legend('Ripple');
        % trace = us_V(find(us_d>1 & ~us_e));       % Close2bat epochs
        % n = 2^nextpow2(numel(trace));   Y = fft(trace,n);   P = abs(Y/n).^2;    f = Fs/dec*(0:n/2)/n;
        % PSD = P(1:n/2+1);   sm_PSD = smoothdata(PSD,'movmedian',n/(Fs/dec));
        % ax(1) = nexttile;   plot(f,mag2db(sm_PSD),'.m','LineWidth',2);     grid on;    xlabel('Hz');   legend('Far From Bat');
        
        %plot(us_t,zscore(us_V),us_t,us_r);
        
        % %=== Look at FFT depending on acceleration signal
        % plot(us_a); grid on;    hold on;    plot(smoothdata(us_a,'movmean',1*Fs/dec)); hold off;
        % plot(us_a-smoothdata(us_a,'movmean',1*Fs/dec));
    end
    
end

save([analysis_directory,'\SWR_Visual_inpection_AF.mat'],'confirmed_by_inspection');

figHandles = findall(0,'Type','figure');
for i = 1:numel(figHandles)
    saveas(figHandles(i),[analysis_directory, '/', [batdate, '_CSC',  num2str(CSC_ch)] '_figure' num2str(numel(figHandles)+1-i) '.png']);
end
close all;
 
%% Function for looking at voltage, velocity, acceleration, etc...

function Look_at_Ripples(signal,t_vector,Fs,t,v_abs,a_abs,s)

%=== Look at 60s chunks
figure('units','normalized','outerposition',[0 0.1 1 0.7]);
tiledlayout(6,10,'TileSpacing','none');
chunk_start_smp = round(1:60*Fs:t_vector(end)*Fs-1);

for i=1:numel(chunk_start_smp)-1
    
    %=== Filter definitions and params
    RP_F_band = [100 200];                  % Ripples frequency band (From MY: [80 160])
    SP_F_band = [600 6000];                 % Spikes frequency band
    passband_norm = SP_F_band./(0.5*Fs);    % Normalize over sampling frequency
    [b,a] = ellip(4,0.1,70,passband_norm);  % Define filter transfer function coefficients
    RP_th = 3;                              % Threshold on zscored Ripple power
    
    %=== Assign zero velocity and acceleration if not provided
    if isempty(t),t=t_vector;end
    if isempty(v_abs),v_abs=zeros(size(t_vector));end
    if isempty(a_abs),a_abs=zeros(size(t_vector));end
    
    %=== Get raw signal, spikes and ripple power
    raw_signal = signal(:,chunk_start_smp(i):chunk_start_smp(i+1)-1);
    hps_signal = filtfilt(b,a,raw_signal);  hps_signal(abs(hps_signal)>500)=0;
    rpl_signal = zscore(abs(hilbert(bandpass(raw_signal,RP_F_band,Fs))));
    time = t_vector(chunk_start_smp(i):chunk_start_smp(i+1)-1);  time = time-time(1);
    spt_signal = raw_signal;
    
    %=== Interpolate the velocity/acceleration signals at the corresponding timestamps
    [~,strt] = min(abs(t-t_vector(chunk_start_smp(i))));
    [~,stop] = min(abs(t-t_vector(chunk_start_smp(i+1)-1)));
    vel_signal = interp1(t(strt:stop),v_abs(strt:stop),t_vector(chunk_start_smp(i):chunk_start_smp(i+1)-1));
    acc_signal = interp1(t(strt:stop),a_abs(strt:stop),t_vector(chunk_start_smp(i):chunk_start_smp(i+1)-1));
    
    %=== Get spike times
    coord_x = [];   coord_y = [];
    if ~isempty(s)
        for nc = 1:size(s,1)
            s{nc} = s{nc}-t_vector(chunk_start_smp(i));
            spikes = s{nc}(s{nc}>0 & s{nc}<time(end));
            coord_x = cat(2,coord_x,[spikes';spikes']);
            coord_y = cat(2,coord_y,[ones(size(spikes'))-nc;zeros(size(spikes'))-nc]);
        end
    end
    
    %=== Plot
    ax(1) = nexttile(1,[1 8]);     plot(time,raw_signal,'k');                                                   xticks([]); ylabel('Raw Voltage (uV)');
    ax(2) = nexttile(11,[1 8]);    plot(coord_x,coord_y,'m-','LineWidth',1);                                    xticks([]); ylabel('Sorted units');
    ax(3) = nexttile(21,[1 8]);    plot(time,hps_signal,'r');                                                   xticks([]); ylabel('High Pass (uV)');
    ax(4) = nexttile(31,[1 8]);    plot(time,rpl_signal,'b');                                                   xticks([]); ylabel('Ripple Power (zscore)');    hold on;    refline(0,RP_th);   hold off;
    ax(5) = nexttile(41,[1 8]);
    yyaxis left;    plot(time,vel_signal,'g','LineWidth',3);   ylabel('velocity (m/s)');    yyaxis left;    ylim([0 10]);
    yyaxis right;   plot(time,acc_signal,'c','LineWidth',1);   ylabel('acceleration (g)');  yyaxis right;   ylim([0 4]);
    ax(6) = nexttile(51,[1 8]);    spectrogram_AF_v0(spt_signal,Fs,0.5,0.05,[70 180]);      xlabel('Time(s)');  ylabel('Frequency (Hz)');
    linkaxes(ax,'x');       xlim([time(1),time(end)]);
    
    %=== Calculate and plot the PSD
    x = raw_signal; n = 2^nextpow2(length(x));
    Y = fft(raw_signal,n);  f = Fs*(0:(n/2))/n; P = abs(Y/n).^2;
    nexttile(9,[6 2]);   plot(f,P(1:n/2+1),'-k'); hold on;   plot(f,smoothdata(P(1:n/2+1),'movmedian',3*n/Fs),'r','LineWidth',4);
    xlim([1 180]);    set(gca, 'YScale', 'log');    hold off;   xlabel('Frequency (Hz)');   h = gca;    h.YAxis.Visible = 'off';
    disp(['Chunk ',num2str(i)]);    pause;
    
end
end
