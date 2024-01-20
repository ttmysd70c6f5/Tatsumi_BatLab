function Process_extracted_cdp_TY_v0()
%% Function for the pre-processing of collective behavior (March 2022)
% This script was modified from the script
% "Process_Collective_Behavior_AF_v0.m" written by Angelo Forli.
%
% Requires extracted data from Ciholas recordings

disp('...Processing Ciholas Data...');

%=== Load data and group name
extracted_CDPfile = dir(fullfile(cd, '*extracted_*'));          load(extracted_CDPfile.name);
batdate = extracted_CDPfile.name(11:16);

% %=== Ad-hoc corrections:
% if strcmp(extracted_CDPfile.name,'extracted_211229_cdp_fly_1.mat') || ...
%    strcmp(extracted_CDPfile.name,'extracted_210720_cdp_fly_1.mat') || ...
%    strcmp(extracted_CDPfile.name,'extracted_210720_cdp_fly_2.mat') || ...
%    strcmp(extracted_CDPfile.name,'extracted_210720_cdp_fly_3.mat')
%     CDPmtdata.tags = 6;
%     for i = 1:CDPmtdata.tags
%         tag_data{1,i} = tag_data{1,1};
%         tag_data_filt{1,i} = tag_data_filt{1,1};
%         tag_ac_data{1,i} = tag_ac_data{1,1};
%     end
%     for i = 1:CDPmtdata.tags
%         if i~=6
%             tag_data{1,i}(:,3:5) = tag_data{1,i}(:,3:5) + 1e-3*randn(size(tag_data{1,i}(:,3:5)));            % Add a bit of noise
%             tag_data_filt{1,i}(:,3:5) = tag_data_filt{1,i}(:,3:5) + 1e-3*randn(size(tag_data_filt{1,i}(:,3:5)));
%             tag_ac_data{1,i}(:,3:5) = tag_ac_data{1,i}(:,3:5) + 1e-3*randn(size(tag_ac_data{1,i}(:,3:5)));
%         end
%     end
% end

%=== Parameters and metadata
n_tags = CDPmtdata.tags;
r_lim = [-2.9 2.9; -2.6 2.6; 0 2.30];                                                                           %Room boundaries
edges_d = {r_lim(1,1):(r_lim(1,2)-r_lim(1,1))/10:r_lim(1,2) r_lim(2,1):(r_lim(2,2)-r_lim(2,1))/10:r_lim(2,2)};  %Edges for density histogram
Fs = 100;                                                                                                       %Sampling frequency (Hz) for common time
% if Group_name == 'D'
%     bat_nms = ['Dai'; 'Den'; 'Dia'; 'Dor'; 'Dum'; 'Im1'; 'Im2'];
% else
%     bat_nms = ['Far'; 'Fer'; 'Fia'; 'For'; 'Ful'; 'Im1'; 'Im2'];
% end
% bat_nms = bat_nms(1:n_tags,:);                                                                                  %Bat Names
bat_nms =   ['bt01';'bt02';'bt03';'bt04';'bt05';'bt06';'bt07';'bt08';'bt09';'bt10'];
bat_pairs = nchoosek(1:n_tags,2);                                                                               %Bat Pairs
bat_pair_nms = [bat_nms(bat_pairs(:,1),:), '-'.*ones(length(bat_pairs),1), bat_nms(bat_pairs(:,2),:)];          %Bat Pairs Names
bat_clr = lines(n_tags);                                                                                        %Bat Colors
v_th = 0.5;                                                                                                     %Velocity threshold (m/s) for flight segmentation

%=== Options
options.show_fig = 1;
options.save_data = 1;
options.use_r_corr = 1;
options.savemovie = 0;

fig_cnt = 1;
fig_dir = dir(fullfile(pwd,'Ext_Behavior*','*figure*.png'));
for i = 1:length(fig_dir)
    figstr = split(fig_dir(i).name,["figure",".png"]);
    if fig_cnt <= str2double(figstr{end-1})
        fig_cnt = str2double(figstr{end-1}) + 1;
    end
end

%=== Custom graded colormap(level,RGB,bat)
for i = 1:n_tags
    for j = 1:3
        custom_map(:,j,i) = linspace(1,bat_clr(i,j))';
    end
end

%=== Analysis folder for storing the results
if options.save_data
    analysis_directory=fullfile(pwd,['Ext_Behavior_',datestr(now, 'yymmdd_HHMM')]);
    if ~exist(analysis_directory,'dir')
        mkdir(analysis_directory);
    end
end

%% Remove duplicate samples and shift Ciholas Time by 1.05s (t = 0 corresponds to the first '3s' Master-9 TTL)

for i = 1:n_tags
    [~,ia,~] = unique(tag_data{1,i}(:,8),'stable');         tag_data{1,i} = tag_data{1,i}(ia,:);              tag_data{1,i}(:,8) = tag_data{1,i}(:,8)+1.05;
    [~,ia,~] = unique(tag_data_filt{1,i}(:,8),'stable');    tag_data_filt{1,i} = tag_data_filt{1,i}(ia,:);    tag_data_filt{1,i}(:,8) = tag_data_filt{1,i}(:,8)+1.05;
    [~,ia,~] = unique(tag_ac_data{1,i}(:,8),'stable');      tag_ac_data{1,i} = tag_ac_data{1,i}(ia,:);        tag_ac_data{1,i}(:,8) = tag_ac_data{1,i}(:,8)+1.05;
end

%% Ad-hoc corrections:
% 
% if strcmp(extracted_CDPfile.name,'extracted_211213_cdp_fly_1.mat')
%     for i = 1:n_tags
%         tag_data{1,i}(:,8) = tag_data{1,i}(:,8)-7.22;            % Determined by Cross-correlating the voltage artifacts and the acceleration signal from implanted bat
%         tag_data_filt{1,i}(:,8) = tag_data_filt{1,i}(:,8)-7.22;
%         tag_ac_data{1,i}(:,8) = tag_ac_data{1,i}(:,8)-7.22;
%     end
% elseif strcmp(extracted_CDPfile.name,'extracted_211209_cdp_fly_1.mat')
%     for i = 1:n_tags
%         tag_data{1,i}(:,3:5) = tag_data{1,i}(:,3:5)+1e-5*randn(size(tag_data{1,i}(:,3:5)));            % Add a bit of noise
%         tag_data_filt{1,i}(:,3:5) = tag_data_filt{1,i}(:,3:5)+1e-5*randn(size(tag_data_filt{1,i}(:,3:5)));
%         tag_ac_data{1,i}(:,3:5) = tag_ac_data{1,i}(:,3:5)+1e-5*randn(size(tag_ac_data{1,i}(:,3:5)));
%     end
% elseif strcmp(extracted_CDPfile.name,'extracted_211205_cdp_fly_1.mat')
%     %figure('units','normalized','outerposition',[0 0.3 1 0.5]);
%     for i = 1:n_tags
%         for j=3:5
%             %temp = tag_data{1,i}(:,j);
%             tag_data{1,i}(103090:105110,j) = NaN;
%             tag_data{1,i}(119100:123850,j) = NaN;
%             %plot(temp);    hold on; plot(tag_data{1,i}(:,j));    hold off;  title(num2str(i));
%             %xlim([103090-1e3 123850+1e3]);
%             %[x,y] = ginput;
%             %tag_data{1,i}(round(x),j) = y;
%             tag_data{1,i}(:,3:5) = fillmissing(tag_data{1,i}(:,3:5),'linear',1);
%         end
%         tag_data_filt{1,i}(103090:105110,3:5) = tag_data{1,i}(103090:105110,3:5);
%         tag_data_filt{1,i}(119100:123850,3:5) = tag_data{1,i}(119100:123850,3:5);
%     end
% end

%% Calculate evenly sampled kinematic variables (r,v,a) and wing-beats epochs

t = [CDPmtdata.TTL_times(1):1/Fs:CDPmtdata.TTL_times(end)]';       %Evenly sampled time vector from first to last TTL.
T = length(t);                                                     %Number of time samples
r = zeros(T,3,n_tags);                                             %3D position vector (sample,dim,id)

%=== Interpolate position at evenly spaced time points
for i = 1:n_tags
    r(:,:,i) =  interp1(tag_data_filt{1,i}(:,8), tag_data_filt{1,i}(:,[3:5]),t,'nearest','extrap'); %DO NOT use spline interpolation!
end

%=== Interpolate acceleration at evenly spaced time points
if exist('tag_ac_data')
    a = zeros(T,3,n_tags);
    for i = 1:n_tags
        a(:,:,i) =  interp1(tag_ac_data{1,i}(:,8), tag_ac_data{1,i}(:,[3:5]),t,'nearest','extrap'); %DO NOT use spline interpolation!
    end
    a_abs = squeeze(vecnorm(a,2,2));    %Modulus
    a_flt = bandpass(a_abs,[7 9],100);  %Filtered at the wing-beat frequency
end

%=== Calculate velocity and 2D-direction of motion
v = diff(r,1,1)./diff(t); v=cat(1,zeros(1,3,n_tags),v);   v_abs = squeeze(vecnorm(v,2,2));
% angle = squeeze(heading_angle(v(:,1,:),v(:,2,:)));

%=== Detect flight epochs based on wing-beat signal
if exist('tag_ac_data')
    for i = 1:n_tags
        [up,lo] = envelope(a_flt(:,i)+normrnd(0,1e-3,length(a_flt),1),10,'peak');   %Envelope of the acceleration signal (noise addedd to avoid problems with splining)
        env = normalize(up - lo,'range');                                           %Amplitude of the envelope
        env_th = otsuthresh(histcounts(env));                                       %Threshold (based on Otsu method). Can be set at 0.35
        wBeats(:,i) = movsum(env>env_th,2*Fs)>Fs/5;                                 %Euristic criterion for flight detection
%                         ax(i) = subplot(n_tags,1,i);  area(t,wBeats(:,i)*1,'FaceAlpha',0.3,'LineStyle','none');  hold on;
%                         plot(t,normalize(v_abs(:,i),'range'));
%                         plot(t,r(:,1,i),t,r(:,2,i));
%                         plot(t,normalize(a_flt(:,i),'range',[-1 1]));
%                         plot(t,normalize(movsum(env>env_th,2*Fs),'range'));
%                         hold off;
    end
%               linkaxes(ax,'x');
end

% % Uncomment to control interpolation quality
% ax(1) = subplot(311);  plot(tag_data{1,4}(:,8),tag_data{1,4}(:,3)','.',t,r(:,1,4),'-');
% ax(2) = subplot(312);  plot(tag_data{1,4}(:,8),tag_data{1,4}(:,4)','.',t,r(:,2,4),'-');
% ax(3) = subplot(313);  plot(tag_data{1,4}(:,8),tag_data{1,4}(:,5)','.',t,r(:,3,4),'-');
% linkaxes(ax,'x');
% ax(1) = subplot(311);  plot(tag_ac_data{1,1}(:,8),tag_ac_data{1,1}(:,3)','.',t,a(:,1,1),'-');
% ax(2) = subplot(312);  plot(tag_ac_data{1,1}(:,8),tag_ac_data{1,1}(:,4)','.',t,a(:,2,1),'-');
% ax(3) = subplot(313);  plot(tag_ac_data{1,1}(:,8),tag_ac_data{1,1}(:,5)','.',t,a(:,3,1),'-');
% linkaxes(ax,'x');

%% Correct position by median filtering when the bat is not flapping its wings

if exist('tag_ac_data')
    %=== Find stationary epochs
    stat_periods = repmat(~wBeats,1,1,3);    stat_periods = permute(stat_periods,[1 3 2]);
    r_stat = r;    r_stat(~stat_periods) = nan;
    
    %=== Filter position during stationary epochs
    for i = 1:n_tags
        r_stat(:,:,i) =  smoothdata(r_stat(:,:,i),1,'movmedian',Fs*5,'omitnan');
    end
    r_stat(~stat_periods) = nan;
    
    % % Uncomment to control filtering quality
    % ax(1) = subplot(311);  plot(tag_data{1,4}(:,8),tag_data{1,4}(:,3)',t,r_stat(:,1,4));
    % ax(2) = subplot(312);  plot(tag_data{1,4}(:,8),tag_data{1,4}(:,4)',t,r_stat(:,2,4));
    % ax(3) = subplot(313);  plot(tag_data{1,4}(:,8),tag_data{1,4}(:,5)',t,r_stat(:,3,4));
    % linkaxes(ax,'x');
    
    %=== Substitute median filtered data when bat is not flapping its wings
    r_stat = fillmissing(r_stat,'constant',0,1);
    r_corr = r_stat.*stat_periods+r.*(~stat_periods);
    r_corr = smoothdata(r_corr,1,'loess',Fs*1);     %DO NOT USE lowess!!
    
    % ax(1) = subplot(311);  plot(tag_data{1,4}(:,8),tag_data{1,4}(:,3)',t,r_corr(:,1,4));
    % ax(2) = subplot(312);  plot(tag_data{1,4}(:,8),tag_data{1,4}(:,4)',t,r_corr(:,2,4));
    % ax(3) = subplot(313);  plot(tag_data{1,4}(:,8),tag_data{1,4}(:,5)',t,r_corr(:,3,4));
    % linkaxes(ax,'x');
    
else
    %Aggressive median filtering on original data
    tag_data_stat = tag_data;
    r_stat = zeros(T,3,n_tags);
    for i = 1:n_tags
        tag_data_stat{1,i}(:,[3:5]) = smoothdata(tag_data_filt{1,i}(:,[3:5]),1,'movmedian',Fs*3);
        r_stat(:,:,i) =  csaps(tag_data_stat{1,i}(:,8), tag_data_stat{1,i}(:,[3:5])', 1, t)';
    end
    %Calculate xy-velocity and smooth
    v_filt = smoothdata(squeeze( vecnorm( v(:,1:2,:),2,2)), 1, 'movmedian', Fs*3);
    %Substitute median filtered data when velocity is less than threshold
    stat_periods = repmat((v_filt < v_th),1,1,3);    stat_periods = permute(stat_periods,[1 3 2]);
    r_corr = r_stat.*stat_periods+r.*(~stat_periods);
    r_corr = smoothdata(r_corr,1,'loess',Fs*1);
    r_corr_1 = smoothdata(r_corr,1,'loess',Fs*1);
end

%% Correct position, recalculate velocity and heading angle

if options.use_r_corr
    r_old = r;  v_old = v_abs;  %Store old versions, in case you need them
    r = r_corr; 
    v = diff(r,1,1)./diff(t); v=cat(1,zeros(1,3,n_tags),v);   v_abs = squeeze(vecnorm(v(:,1:3,:),2,2)); 
    % angle = squeeze(heading_angle(v(:,1,:),v(:,2,:)));
end

%% Flight segmentation

%=== Detect flight starts and stops by using risetime and falltime
%!!! Be careful with the levels of risetime and falltime, these
%    are fundamental in order to exclude stationary tails in the flight
bflying = zeros(T,n_tags); f_num = zeros(n_tags,1);   f_smp = cell(n_tags,2);
for i = 1:n_tags
    [bflying(:,i),f_num(i),f_smp{i,1},f_smp{i,2}] = FlightSegm_AF_v0(v_abs(:,i),v_th,Fs);
end

%Define staionary periods when the bat is not flying
stat_periods = repmat(~bflying,1,1,3);   stat_periods = permute(stat_periods,[1 3 2]);
r_qt = r;   r_qt(~stat_periods)=nan;     
% angle(squeeze(stat_periods(:,1,:))) = nan;

%% Make a few figures

if options.show_fig
    
    %=== FIGURE: Raw Position VS Corrected Position
    for j = 1:3
        figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        tiledlayout(n_tags,1,'TileSpacing','tight');
        for i=1:n_tags
            ax(i) = nexttile; 
            plot(tag_data{1, i}(:,8), tag_data{1, i}(:,2+j),t,r_corr(:,j,i));
            sgtitle(['D-' num2str(j)]);   ylim(r_lim(j,:));   legend('raw','corrected');    ylabel('m'); 
            
        end
        linkaxes(ax,'x');    xlabel('Time(s)'); xlim([0,t(T)]);

    end
     
    %=== FIGURE: Raw Velocity VS Corrected Velocity
    figure();   set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(n_tags,1,'TileSpacing','tight');
    for i=1:n_tags
        ax(i) = nexttile;   
        plot(t,v_old(:,i),'.');    hold on;     sgtitle('v');
        plot(t,v_abs(:,i),'.');    hold off;    legend('raw','corrected');      ylabel('m/s');
    end
    linkaxes(ax,'x');    xlabel('Time(s)'); xlim([0,t(T)]);  

    
    %=== FIGURE: Scatter plot all bats
    figure();       set(gcf, 'units','normalized','outerposition',[0 0.25 1 0.5]);
    for i=1:n_tags
        subplot(131);  plot3(r(:,1,i),r(:,2,i),r(:,3,i),'Color', bat_clr(i,:));  xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));  title('3D view');                 hold on;  axis equal;
        subplot(132);  plot3(r(:,1,i),r(:,2,i),r(:,3,i),'Color', bat_clr(i,:));  xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));  title('Top view');   view(0,90);  hold on;  axis equal;
        subplot(133);  plot3(r(:,1,i),r(:,2,i),r(:,3,i),'Color', bat_clr(i,:));  xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));  title('Door view');  view(-30,0); hold on;  axis equal;
    end
    hold off;
    

    %=== FIGURE: Trajectories individual bats
    figure();   set(gcf, 'units','normalized','outerposition',[0.2 0.25 0.7 0.35]);
    tiledlayout(1,n_tags,'TileSpacing','none');
    for i=1:n_tags
        nexttile;  plot3(r(:,1,i),r(:,2,i),r(:,3,i),'-','Color', bat_clr(i,:));  xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:));  title(bat_nms(i,:));  view(90,90);
        axis equal;
    end

    %=== FIGURE: Density Histograms heat-map
    figure();   set(gcf, 'units','normalized','outerposition',[0 0.25 1 0.35]);
    for i=1:n_tags
        sbpt = subplot(1,n_tags,i);
        % hist3(r(:,1:2,i),'edges',edges_d,'CdataMode','auto','EdgeColor','none','FaceColor','interp');
        hist3(r(:,1:2,i),'edges',edges_d,'CdataMode','auto','FaceColor','interp');
        xlabel('x');
        xlim(r_lim(1,:)); ylim(r_lim(2,:));   
        title(bat_nms(i,:));
        title(['bat' num2str(i)])
        view(90,90);  colormap(sbpt,custom_map(:,:,i)); % Change color scheme
        axis square;
    end
    

    %=== FIGURE: Velocity and flight segmentation
    figure();       set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(n_tags,1,'TileSpacing','tight');
    for i = 1:n_tags
        ax(i) = nexttile;   %ax(i) = subplot(n_tags,1,i);
        area(t,bflying(:,i)*5,'FaceAlpha',0.3,'LineStyle','none');  hold on;
        area(t,wBeats(:,i)*-1,'FaceAlpha',0.3,'LineStyle','none');  refline(0,-0.3);
        plot(t,v_abs(:,i),'.','Color', bat_clr(i,:));     plot(t,r(:,1,i),'k--');  ylabel('Velocity (m/s)');     hold off;
        legend('Fly','Wing-B','Vel','x(m)');
        title([num2str(f_num(i)) ' flights']);
    end
    linkaxes(ax,'x');   xlabel('Time(s)');
    

end

%% Save figures and data
if options.save_data
    figHandles = findall(0,'Type','figure');
    for i = 1:numel(figHandles)
        saveas(figHandles(i),[analysis_directory, '/', batdate '_figure' num2str(numel(figHandles)+1-i) '.png']);
    end
    close all;
    save([analysis_directory,'/Extracted_Behavior_', batdate, '.mat'],...
        'a','a_abs','a_flt','bat_clr','bat_nms','bat_pair_nms','bat_pairs',... % angle
        'batdate','bflying','CDPmtdata','edges_d','env','extracted_CDPfile',...
        'f_num','f_smp','Fs','n_tags','options',... % Group_names
        'r','r_lim','r_old','r_qt','stat_periods',...
        't','T','v','v_abs','v_th','wBeats');
end

%% Save movie at 10Hz with bat trajectories
% if options.savemovie
%     saveCollMovie_AF_v1(r,Fs,bat_nms,r_lim,analysis_directory,[batdate]);
% end
end

