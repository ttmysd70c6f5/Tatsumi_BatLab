function ExtractCdp_AF_TY_v2(fname,outdir,use_sync,varargin)
%% Extract tracking data from Ciholas RTLS acquisition
% Data are recorded through a python script and saved in a ASCII file
% Data columns (comma delimited) are defines as follows:

% Column 1 = pos(1), sync(2) or acc(0) type index
% Column 2 = device Serial Number
% Column 3 = Network time (each tic is ~15.65 ps)
% Column 4,5,6 = x,y,z (position in mm or acceleration-sync)
% Column 7 = signal quality for pos or scale for acc
% Column 8 = number of receiving anchors for a tag or '0' for acc/sync

% The sync signal is driven by a TTL:
% 50 ms duration happening every 21, 13, 8, 5, 4s
%-----------------------------------------------------------------

%   *********USAGE EXAMPLES*****************

% Simply call the function by using a unique part of the filename to be processed 

% ExtractCdp_AF_v1('_cdp_fly_1');                                                                   Extract data from all detected serial numbers (excluding synch)
% ExtractCdp_AF_v1('_cdp_fly_1','Exclude',17106999);                                                Same as above, excluding specific serial numbers
% ExtractCdp_AF_v1('_cdp_fly_1','Include',[17106917 17106934 17107055 17106969 17106904 17106951]); Extract data ONLY from included serial numbers

% ATTENTION !! Please be aware that the time vector will be shifted by 1.05s, relative to the M-9 3s TTL

%Open a file to record a log
fid = fopen(fullfile(outdir,'extract_log.txt'),'w');

%Load data, keep non-zero entries, sort according to ascending network times and delete duplicates
file_name = dir(fullfile(cd, ['*' fname '*']));   file_name = file_name.name;
RTLS_data = load(file_name);
RTLS_data = RTLS_data(RTLS_data(:,2)~=0,:);
RTLS_data = sortrows(RTLS_data,3);
RTLS_data = unique(RTLS_data,'rows','stable');
detected_devices = unique(RTLS_data(:,2));

% % Ad hoc corrections
% if strcmp(file_name,'211218_cdp_fly_2.txt')
%     RTLS_data((RTLS_data(:,3)<6.3042e+14),:) = [];
% end

%User inputs overrides
tags_SN = detected_devices;
if nargin > 1
    eSN = [];   iSN = [];
    nparams=length(varargin);
    for i=1:2:nparams
        switch (varargin{i})
            case 'Exclude'
                eSN=varargin{i+1};
                tags_SN(any(tags_SN == eSN,2)) = [];
            case 'Include'
                iSN = varargin{i+1};
                tags_SN = [];   tags_SN = iSN';
        end
    end
end

%Parameters
% use_sync = 1;
rec_duration = 10800;                                                                           %approx rec duration (s), 10800 covers 3 hours
ntw_time = 15.65e-12;                                                                           %network tic interval (s)
sync_SN = 17040920; %17106963;                                                                  %sync tag serial number
if use_sync, tags_SN = setdiff(tags_SN,sync_SN,'stable');end
n_tags = length(tags_SN);
TTL_time_diff = [21; 13; 8; 5; 4];                                                              %TTL delays in s
TTL_abs_times = [0; cumsum(repmat(TTL_time_diff,round(rec_duration*2/sum(TTL_time_diff)),1))];
CDPmtdata.Fs = 100;                                                                             %Acquisition frequency CDP (Hz)
CDPmtdata.tag_SN = tags_SN;
CDPmtdata.sync_SN = sync_SN;


%Extract data and log CDP metadata, start from sync
sync_data = RTLS_data(RTLS_data(:,1)==2 & RTLS_data(:,2)==sync_SN,2:end);
CDPmtdata.sync = ~isempty(sync_data);
CDPmtdata.sync_duration = (sync_data(end,2)-sync_data(1,2))*ntw_time/60;
CDPmtdata.sync_samples = length(sync_data);
%Tag Data
j = 0;
for i = 1:n_tags
    if ~isempty(find(RTLS_data(:,2)==tags_SN(i)))
        j = j+1;
        tag_data{j} = RTLS_data(RTLS_data(:,1)==1 & RTLS_data(:,2)==tags_SN(i),2:end);
        tag_data{1,j}(:,[3:5]) = tag_data{1,j}(:,[3:5])/1000;
        CDPmtdata.tag_duration(j) = (tag_data{1,j}(end,2)-tag_data{1,j}(1,2))*ntw_time/60;
        CDPmtdata.tag_samples(j) = length(tag_data{1,j});
        
        if ~isempty(find((RTLS_data(:,1)==0 & RTLS_data(:,2)==tags_SN(i))))
            tag_ac_data{j} = RTLS_data(RTLS_data(:,1)==0 & RTLS_data(:,2)==tags_SN(i),2:end);
            %Acceleration values need to be converted
            indexes = [false(size(tag_ac_data{1,j}(:,1:2))), tag_ac_data{1,j}(:,[3:5])>double(intmax('int32'))];
            tag_ac_data{1,j}(indexes) = tag_ac_data{1,j}(indexes)-double(intmax('uint32'));
            tag_ac_data{1,j}(:,[3:5]) = tag_ac_data{1,j}(:,[3:5])*2/double(intmax('int32'));
            CDPmtdata.tag_ac_duration(j) = (tag_ac_data{1,j}(end,2)-tag_ac_data{1,j}(1,2))*ntw_time/60;
            CDPmtdata.tag_ac_samples(j) = length(tag_ac_data{1,j});
        end
    end
end

CDPmtdata.tags = n_tags;
disp(CDPmtdata);

%% Plot raw data (x position and sync + x acceleration and sync) referred to the first acquired sample

t0 = min(RTLS_data(:,3));

figure('units','normalized','outerposition',[0.5 0 0.5 1]);   %x_axis
for i=1:n_tags
    ax_0(i) = subplot(n_tags+1,1,i);    
    plot((tag_data{1, i}(:,2)-t0)*ntw_time/60, tag_data{1, i}(:,3));    ylabel('m');    title(['Tag',num2str(i)]);              
end
ax_0(n_tags+1) = subplot(n_tags+1,1,n_tags+1);
plot((sync_data(:,2)-t0)*ntw_time/60, sync_data(:,3));  ylabel('Synch');              
linkaxes(ax_0,'x');    xlabel('Time(s)');   sgtitle('x (POSITION)');

figure('units','normalized','outerposition',[0.5 0 0.5 1]);   %x_axis acceleration
for i=1:n_tags
    ax_a(i) = subplot(n_tags+1,1,i);
    plot((tag_ac_data{1, i}(:,2)-t0)*ntw_time/60, tag_ac_data{1, i}(:,3));  ylabel('g');    title(['Tag',num2str(i)]);             
end
ax_a(n_tags+1) = subplot(n_tags+1,1,n_tags+1);
plot((sync_data(:,2)-t0)*ntw_time/60, sync_data(:,3));  ylabel('Synch');              
linkaxes(ax_a,'x');    xlabel('Time(s)');   sgtitle('x (ACCELERATION)');


%% Interpolate time stamps (s) according to sync signal or network time
% [~, TTL_findpeak] = findpeaks(normalize(sync_data(:,3),'zscore'),sync_data(:,2),'MinPeakHeight',1,'MinPeakDistance',2/ntw_time,'MinPeakWidth',0.048/ntw_time);
% [~,~,TTL_risetime] =   risetime(normalize(sync_data(:,3),'zscore'),sync_data(:,2));
% [~,TTL_risetime,~] =   risetime(normalize(sync_data(:,3),'zscore'),sync_data(:,2));
% mean(TTL_findpeak-TTL_risetime)

if use_sync
    %Detect network times corresponding to peaks in the sync signal separated by at least 2s
    % [~, TTL_network_times] = findpeaks(normalize(sync_data(:,3),'zscore'),sync_data(:,2),'MinPeakHeight',1,'MinPeakDistance',2/ntw_time,'MinPeakWidth',0.048/ntw_time); % Avoid an artifact by turning the Master-9 off

    % Detect false-TTL (=artifact) based on the width of TTL signals
    [~, ~, TTL_width, ~] = findpeaks(normalize(sync_data(:,3),'zscore'),sync_data(:,2),'MinPeakHeight',1,'MinPeakDistance',2/ntw_time);
    med_cdp =   median(TTL_width*ntw_time);
    mad_cdp =   mad(TTL_width*ntw_time,1);
    true_TTL    =   ~(TTL_width*ntw_time<med_cdp - 5*mad_cdp);

    %Detect network times corresponding to peaks in the sync signal separated by at least 2s
    % [~, TTL_network_times2] = findpeaks(normalize(sync_data(:,3),'zscore'),sync_data(:,2),'MinPeakHeight',1,'MinPeakDistance',2/ntw_time,'MinPeakWidth',0.048/ntw_time); % Avoid an artifact by turning the Master-9 off

    % Detect network times corresponding to crossing 90% reference level in the sync signal
    [~,~,TTL_network_times] =   risetime(normalize(sync_data(:,3),'zscore'),sync_data(:,2));
    if length(TTL_network_times)==length(true_TTL)
        TTL_network_times       =   TTL_network_times(true_TTL);
    end
    
    %Check if all the TTL were correctly detected
    [~, control] = findpeaks([0; diff(TTL_network_times)*ntw_time],'MinPeakHeight',20);
    if all(unique(diff(control))==5)
        disp('All TTLs were detected, alignment OK');
        % fprintf(fid,'All TTLs were detected, alignment OK\n');
        sync_err = 0;
    %elseif any(unique(round(diff(TTL_network_times)*ntw_time,1),'stable') ~= TTL_time_diff) % This is a double check: 1st TTL interval should be 21s and all TTL intervals == TTL_time_diff
     %   warning('Undetected TTLs, check alignment process');
      %  sync_err = 1;
    else
        warning('Undetected TTLs, check alignment process');
        sync_err = 1;
    end
    figure('units','normalized','outerposition',[0.5 0 0.5 1]);
    subplot(211);   plot(diff(control),'.');    ylabel('TTL count');     xlabel('#TTL pair separated by 21s');
    subplot(212);   findpeaks(normalize(sync_data(:,3),'zscore'),sync_data(:,2),'MinPeakHeight',1,'MinPeakDistance',2/ntw_time,'MinPeakWidth',0.048/ntw_time);   xlabel('Packet#');     ylabel('Sync signal (a.u.)');
    sgtitle('---Synchronization Check---');
   
    %If no errors in TTL detection -> interpolate time from known TTL delays (M-9 is assumed to be the reference time!)
    if ~sync_err
        TTL_abs_times = TTL_abs_times(1:length(TTL_network_times));
        CDPmtdata.TTL_times = TTL_abs_times;
        cdp_t = interp1(TTL_network_times,TTL_abs_times,sync_data(:,2),'linear','extrap');
        sync_data(:,8) = cdp_t;
        for j = 1:n_tags
            tag_data{1,j}(:,8) = interp1(sync_data(:,2),cdp_t,tag_data{1,j}(:,2),'linear','extrap');
            if ~isempty(find((RTLS_data(:,1)==0 & RTLS_data(:,2)==tags_SN(j))))
                tag_ac_data{1,j}(:,8) = interp1(sync_data(:,2),cdp_t,tag_ac_data{1,j}(:,2),'linear','extrap');
            end
        end
        
    %If errors in TTL detection -> use network time between first and last TTL    
    else
        TTL_abs_times = [0 round(TTL_network_times(end)-TTL_network_times(1))*ntw_time];
        CDPmtdata.TTL_times = TTL_abs_times;
        cdp_t = interp1([TTL_network_times(1) TTL_network_times(end)],TTL_abs_times,sync_data(:,2),'linear','extrap');
        sync_data(:,8) = cdp_t;
        warning('First and last TTL will be used for synchronization');
        for j = 1:n_tags
            tag_data{1,j}(:,8) = interp1(sync_data(:,2),cdp_t,tag_data{1,j}(:,2),'linear','extrap');
            if ~isempty(find((RTLS_data(:,1)==0 & RTLS_data(:,2)==tags_SN(j))))
                tag_ac_data{1,j}(:,8) = interp1(sync_data(:,2),cdp_t,tag_ac_data{1,j}(:,2),'linear','extrap');
            end
        end
        
    end
    
%If not using sync signal -> use network time 
else
    sync_data(:,8) = (sync_data(:,2)-sync_data(1,2)).*ntw_time;
    CDPmtdata.TTL_times = TTL_abs_times(TTL_abs_times<sync_data(end,8));
    for j = 1:n_tags
        tag_data{1,j}(:,8) = (tag_data{1,j}(:,2)-sync_data(1,2)).*ntw_time;
        if ~isempty(find((RTLS_data(:,1)==0 & RTLS_data(:,2)==tags_SN(j))))
            tag_ac_data{1,j}(:,8) = (tag_ac_data{1,j}(:,2)-sync_data(1,2)).*ntw_time;
        end
    end
    true_TTL = [];
end

%% Filter position data with Quadratic regression and Save

tag_data_filt = tag_data;
for i = 1:n_tags
    tag_data_filt{1,i}(:,[3:5]) = smoothdata(tag_data_filt{1,i}(:,[3:5]),1,'loess',CDPmtdata.Fs*1);
end

%Save data
if ~isempty(find((RTLS_data(:,1)==0 & RTLS_data(:,2)==tags_SN(1))))
    save([outdir,'\extracted_', file_name(1:end-4), '.mat'],'tag_data','sync_data','CDPmtdata','tag_data_filt','tag_ac_data','true_TTL');
else
    save([outdir,'\extracted_', file_name(1:end-4), '.mat'],'tag_data','sync_data','CDPmtdata','tag_data_filt','true_TTL');
end

%% Plot Raw and Filtered position data

figure('units','normalized','outerposition',[0.5 0 0.5 1]);   %x_axis
for i=1:n_tags
    ax(i) = subplot(n_tags,1,i);
    plot(tag_data{1, i}(:,8), tag_data{1, i}(:,3));              hold on;     
    plot(tag_data_filt{1, i}(:,8), tag_data_filt{1, i}(:,3));    hold off;    ylim([-3 3]);  legend('raw','filt');
    ylabel('m');    title(['Tag',num2str(i)]);
end
linkaxes(ax,'x');    xlabel('Time(s)'); sgtitle('x (POSITION)');

figure('units','normalized','outerposition',[0.5 0 0.5 1]);   %y_axis
for i=1:n_tags
    ax(i) = subplot(n_tags,1,i);
    plot(tag_data{1, i}(:,8), tag_data{1, i}(:,4));              hold on;     
    plot(tag_data_filt{1, i}(:,8), tag_data_filt{1, i}(:,4));    hold off;    ylim([-3 3]);  legend('raw','filt');
    ylabel('m');    title(['Tag',num2str(i)]);
end
linkaxes(ax,'x');    xlabel('Time(s)'); sgtitle('y (POSITION)');

figure('units','normalized','outerposition',[0.5 0 0.5 1]);   %z_axis
for i=1:n_tags
    ax(i) = subplot(n_tags,1,i);
    plot(tag_data{1, i}(:,8), tag_data{1, i}(:,5));              hold on;     
    plot(tag_data_filt{1, i}(:,8), tag_data_filt{1, i}(:,5));    hold off;    ylim([0 2.5]); legend('raw','filt');
    ylabel('m');    title(['Tag',num2str(i)]);
end
linkaxes(ax,'x');    xlabel('Time(s)'); sgtitle('z (POSITION)');

%% Sanity check on aquisition intervals and position values
%Sync
figure('units','normalized','outerposition',[0.5 0.3 0.5 0.5]);   sgtitle('Sync data');
subplot(121);     plot(diff(sync_data(:,2))*ntw_time);         xlabel('Packet#');  ylabel('Network Time Difference(s)');
subplot(122);     histogram(diff(sync_data(:,2))*ntw_time);    set(gca,'YScale','log'); xlabel('Network Time Difference(s)');  ylabel('Counts'); 
disp(['Sync: Fs (network time VS interpolated): ' num2str(1/mean(diff(sync_data(:,2))*ntw_time),4) 'Hz VS ' num2str(1/mean(diff(sync_data(:,8))),4) ' Hz']);


if any(diff(sync_data(:,8))>0.5)
    warning('!! Long acquisition pause detected !!');
end

%Tags
figure('units','normalized','outerposition',[0.5 0 0.5 1]);   sgtitle('Tags data');
for i=1:n_tags
    ax(i) = subplot(n_tags,2,2*i-1);   plot(diff(tag_data{1, i}(:,2))*ntw_time);          ylabel('Network Time Difference(s)');          xlabel('Packet#');
    bx(i) = subplot(n_tags,2,2*i);     histogram(diff(tag_data{1, i}(:,2)*ntw_time));   set(gca,'YScale','log');     ylabel('Counts'); xlabel('Network Time Difference(s)');
    title(['Tag',num2str(i)]);
    disp(['Pos Tag' num2str(i) ': Fs (network time VS interpolated): ' num2str(1/mean(diff(tag_data{1, i}(:,2))*ntw_time),4) ' Hz VS ' num2str(1/mean(diff(tag_data{1, i}(:,8))),4) ' Hz']);
    fprintf(fid,['Pos Tag' num2str(i) ': Fs (network time VS interpolated): ' num2str(1/mean(diff(tag_data{1, i}(:,2))*ntw_time),4) ' Hz VS ' num2str(1/mean(diff(tag_data{1, i}(:,8))),4) ' Hz\n']);
    if any(diff(tag_data{1, i}(:,8))>0.5)
    warning('!! Long acquisition pause detected !!');
end
end
linkaxes(ax,'x');       linkaxes(bx,'x');


%% Plot acceleration
if ~isempty(find((RTLS_data(:,1)==0 & any(RTLS_data(:,2)==tags_SN'))))  
    figure('units','normalized','outerposition',[0.5 0 0.5 1]);   sgtitle('Tags data');
    for i=1:n_tags
        ax_1(i) = subplot(n_tags,2,2*i-1);   plot(diff(tag_ac_data{1, i}(:,2))*ntw_time);          ylabel('Network Time Difference(s)');          xlabel('Packet#');
        bx_1(i) = subplot(n_tags,2,2*i);     histogram(diff(tag_ac_data{1, i}(:,2)*ntw_time));   set(gca,'YScale','log');     ylabel('Counts'); xlabel('Network Time Difference(s)');
        disp(['Acc Tag' num2str(i) ': Fs (network time VS interpolated): ' num2str(1/mean(diff(tag_ac_data{1, i}(:,2))*ntw_time),4) ' Hz VS ' num2str(1/mean(diff(tag_ac_data{1, i}(:,8))),4) ' Hz']);
        fprintf(fid,['Acc Tag' num2str(i) ': Fs (network time VS interpolated): ' num2str(1/mean(diff(tag_ac_data{1, i}(:,2))*ntw_time),4) ' Hz VS ' num2str(1/mean(diff(tag_ac_data{1, i}(:,8))),4) ' Hz\n']);
        
    end
    linkaxes(ax_1,'x');       linkaxes(bx_1,'x');

  

    for j = 1:3
        figure('units','normalized','outerposition',[0.5 0 0.5 1]);   %x_axis
        for i=1:n_tags
            ax_2(i) = subplot(n_tags,1,i);
            plot(tag_ac_data{1, i}(:,8), tag_ac_data{1, i}(:,2+j));   sgtitle(['Dimesion' num2str(j)]);     ylabel('g units')
        end
        linkaxes(ax_2,'x');    xlabel('Time(s)');


    end 
    
end


fclose(fid);