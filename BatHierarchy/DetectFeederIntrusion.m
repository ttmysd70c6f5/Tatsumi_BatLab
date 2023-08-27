function DetectFeederIntrusion()

%% 1. Load data
% Reward signals
RewardFile = dir('Ext_Ephys*/*c3d*.mat');
load(fullfile(RewardFile(1).folder,RewardFile(1).name),'AnalogFrameRate','AnalogSignals') % reward signals from the feeder (1st col)
% Single Unit
% EphysFile = dir('Ext_Ephys*\SingleUnits*\*SingleUnits*.mat');
% load(fullfile(EphysFile.folder,EphysFile.name),'TT_unit')
% Behavior
BehFile = dir('Ext_Behavior*\Extracted_Behavior*.mat');
load(fullfile(BehFile.folder,BehFile.name))
% Position cluster
rClusFile = dir('Ext_Behavior*\clu_pos.mat');
load(fullfile(rClusFile.folder,rClusFile.name),'r_clus_id','clu_rwrd')
% Experiment info
% ExpFile1 = dir('Ext_Ephys*\imp_bat.mat');
% load(fullfile(ExpFile1.folder,ExpFile1.name),'imp_bat')
% Bat identification
idFile = dir('Ordered_IDs.txt');
bat_id = importdata(fullfile(idFile.folder,idFile.name));
% TTL
TTLFile = dir('Ext_Ephys*\TTL_timestamps.mat');
load(fullfile(TTLFile(1).folder,TTLFile(1).name),'TTL_timestamps')


fprintf('\n')
fprintf('***Experimental setup***\n')
fprintf('TTL interval: %d sec\n',diff(find(diff(AnalogSignals(:,2)) > 2,2)) / AnalogFrameRate)
% fprintf('Implanted bat: bat %d\n',imp_bat)
fprintf('Bat IDs: %d\n',bat_id)
fprintf('Length of TTL signals: %.2f min\n',(TTL_timestamps.fly(end) - TTL_timestamps.fly(1)) / 1e6 / 60)
fprintf('Length of behavior recording: %.2f min\n',T/Fs/60)

fprintf('TTL duration in TTL clock: %.2f sec\n',...
    (find(diff(AnalogSignals(:,2))>2,1,'last') - find(diff(AnalogSignals(:,2))>2,1,'first'))/AnalogFrameRate)
fprintf('TTL duration in internal clock: %.2f sec\n',...
    (TTL_timestamps.fly(end) - TTL_timestamps.fly(1)) / 1e6)

%% 2. Check the position of clusters
% clu_list = unique(r_clus_id(:,1:3));
% clu_list = clu_list(~isnan(clu_list));
% 
% figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0.05 0.9 0.85]);
% tiledlayout(2, 4, 'TileSpacing', 'tight');
% for cc = 1:length(clu_list)
%     nexttile
%     hold on
%     for bb = 1:3
%         x = r(r_clus_id(:,bb)==clu_list(cc),1,bb);
%         y = r(r_clus_id(:,bb)==clu_list(cc),2,bb);
%         z = r(r_clus_id(:,bb)==clu_list(cc),3,bb);
%         plot3(x,y,z,'.k')
%         xlim([-3 3])
%         ylim([-3 3])
%         zlim([0 3])
%         view([45 45])
%     end
%     hold off
%     title(sprintf('Cluser %d',clu_list(cc)))
%     xlim([-3 3])
%     ylim([-3 3])
% end
% 
% saveas(gcf,fullfile(BehFile.folder,'Analysis_Tatsumi','figure2.fig'))


%% 4. Extract the reward timings
t_rwrd = find(diff(AnalogSignals(:,1))>2);
t_rwrd = t_rwrd - find(diff(AnalogSignals(:,2))>2,1,'first');
t_rwrd = t_rwrd / AnalogFrameRate; % convert to sec
t_rwrd = round(t_rwrd*2)/2;

%% 5. Bat position when the rewad signals
% ver.1
%
% Rwr = cell(length(t_rwrd),5);
% for tt = 1:length(t_rwrd)
%     t_last = 0;
%     b_last = 0;
%     for bb = 1:3
%         if bflying(t == t_rwrd(tt),bb) == 1
%         else
%             b_now = bb;
%             t_now = find(bflying(t < t_rwrd(tt),bb),1,'last')+1; % who land last = who fly last, when land = when fly last time + 1
%             if t_last < t_now
%                 b_last = b_now;
%                 t_last = t_now;
%             end
%         end
%     end
%     Rwr{tt,1} = t_rwrd(tt);
%     Rwr{tt,2} = b_last;
%     Rwr{tt,3} = bat_id(b_last);
%     Rwr{tt,4} = t(t_last);
%     Rwr{tt,5} = t(t_last+find(bflying(t_last+1:end,b_last),1,'first')-1);
% end
% Rwr = cell2table(Rwr,'VariableNames',["t_r" "b_feed" "f_idf" "f_toff" "f_land"]);
% Rwr.stay = Rwr.f_land - Rwr.f_toff;
% 
% % ver.2
% cnt = 0;
% cnt2 = 0;
% Rwr = cell(length(t_rwrd),5);
% for tt = 1:length(t_rwrd)
%     t_feed = 0;
%     b_feed = 0;
%     dist_feed = 99;
%     for bb = 1:3
%         if bflying(t == t_rwrd(tt),bb) == 1
%         else
%             r_now = r(t == t_rwrd(tt),1:2,bb);
%             dist_now = pdist([r_now;3,1],'euclidean');
%             if dist_now < dist_feed
%                 b_feed = bb;
%                 dist_feed = dist_now;
%                 r_feed = r_now;
%                 t_feed = find(bflying(t < t_rwrd(tt),bb),1,'last')+1; % who land last = who fly last, when land = when fly last time + 1
%             end
%         end
%     end
%     if r_feed(1) < 2
%         cnt = cnt + 1;
%     else
%         cnt2 = cnt2 + 1;
%     end
% 
%     Rwr{tt,1} = t_rwrd(tt);
%     Rwr{tt,2} = b_feed;
%     Rwr{tt,3} = bat_id(b_feed);
%     Rwr{tt,4} = t(t_feed);
%     Rwr{tt,5} = t(t_feed+find(bflying(t_feed+1:end,b_feed),1,'first')-1);
%     Rwr{tt,6} = r_feed(1);
%     Rwr{tt,7} = r_feed(2);
%     Rwr{tt,8} = r_clus_id(t_feed,b_feed);
%     Rwr{tt,9} = r_clus_id(t==t_rwrd(tt),b_feed);
% end
% Rwr = cell2table(Rwr,'VariableNames',["t_r" "b_feed" "f_idf" "f_land" "f_toff" "x" "y" "clus_land" "clus_rwd"]);
% Rwr.stay = Rwr.f_toff - Rwr.f_land;

% ver.3
Rwr = cell(length(t_rwrd),5);
row_rmv = false(size(Rwr,1),1);
for tt = 1:length(t_rwrd)
    t_feed = 0;
    b_feed = 0;
    
    for bb = 1:3
        if ismember(r_clus_id(t==t_rwrd(tt),bb), clu_rwrd)
            t_feed_now = find(bflying(t <= t_rwrd(tt),bb),1,'last')+1; % who land last = who fly last, when land = when fly last time + 1
            if t_feed < t_feed_now
                t_feed = t_feed_now;
                b_feed = bb;
                r_feed = r(t_feed_now,1:2,bb);
            end
        end
    end

    if b_feed ~= 0 && ~isempty(find(bflying(t_feed+1:end,b_feed),1,'first'))
        Rwr{tt,1} = bat_id(b_feed);
        Rwr{tt,2} = b_feed;
        Rwr{tt,3} = t(t_feed);
        Rwr{tt,4} = t_rwrd(tt);
        Rwr{tt,5} = t(t_feed+find(bflying(t_feed+1:end,b_feed),1,'first')-1);
        Rwr{tt,6} = t_rwrd(tt) - t(t_feed);
        Rwr{tt,7} = r_feed(1);
        Rwr{tt,8} = r_feed(2);
        Rwr{tt,9} = r_clus_id(t_feed,b_feed);
        Rwr{tt,10} = r_clus_id(t==t_rwrd(tt),b_feed);
    else
        row_rmv(tt) = true;
    end
end

Rwr(row_rmv,:) = [];
Rwr = cell2table(Rwr,'VariableNames',["feeding_id" "b_feed" "t_land" "t_r" "t_toff" "wait" "x" "y" "clus_land" "clus_rwd"]);
Rwr.stay = Rwr.t_toff - Rwr.t_land;

%% 6. Remove overlapped extraction
row_rmv = false(size(Rwr,1),1);
stay_now = Rwr.stay(1);
for i = 2:size(Rwr,1)
    if stay_now == Rwr.stay(i)
        row_rmv(i) = true;
    else
        stay_now = Rwr.stay(i);
    end
end
fprintf('%d out of %d is overlapped\n', sum(row_rmv), size(Rwr,1))
Rwr(row_rmv,:) = [];

%% 7. Remove weired waiting time
Rwr(Rwr.wait > 5,:) = [];

%% 7. Evaluate the algorithm to extract feeding bats
% figure
% hold on
% for i=1:size(Rwr,1)
% %     plot3(r(t>=Rwr.f_toff(i)&t<=Rwr.f_land(i),1,Rwr.b_feed(i)),...
% %         r(t>=Rwr.f_toff(i)&t<=Rwr.f_land(i),2,Rwr.b_feed(i)),...
% %         r(t>=Rwr.f_toff(i)&t<=Rwr.f_land(i),3,Rwr.b_feed(i)),'ok')
% %     plot(r(t>=Rwr.f_toff(i)&t<=Rwr.f_land(i),1,Rwr.b_feed(i)),...
% %         r(t>=Rwr.f_toff(i)&t<=Rwr.f_land(i),2,Rwr.b_feed(i)),'.k')
%     plot(r(t==Rwr.t_r(i),1,Rwr.b_feed(i)),...
%         r(t==Rwr.t_r(i),2,Rwr.b_feed(i)),'.k')
% end
% hold off
% xlim([-3 3])
% ylim([-3 3])
% title('Position of feeding bat at reward timing')
% xlabel('x')
% ylabel('y')
% 
% groupcounts(Rwr,["clus_rwd"])

%% Identify intruders
% Notes: this script ignore the case that there are two intruders. Be
% careful when such cases are included in the further analysis.
Rwr.b_intruder = zeros(size(Rwr,1),1);
Rwr.intruder_id = zeros(size(Rwr,1),1);
Rwr.intruded = false(size(Rwr,1),1);
Rwr.mlt_intr = false(size(Rwr,1),1);
Rwr.feeder_win = false(size(Rwr,1),1);
Rwr.feeder_lose = false(size(Rwr,1),1);
Rwr.t_intr_come = zeros(size(Rwr,1),1);
Rwr.t_intr_leave = zeros(size(Rwr,1),1);

for tt = 1:size(Rwr,1)
    cnt_intr = 0;
    b_feeding = Rwr.b_feed(tt);
    b_intruder = 0;
    for bb = 1:3
        if bb ~= b_feeding
            is_intruded = ismember(r_clus_id(find(t==Rwr.t_r(tt))+1:find(t==Rwr.t_toff(tt)),bb),clu_rwrd);
            if sum(is_intruded) == 0 % no intruder
            elseif sum(diff([is_intruded(1);is_intruded]) == 1) == 0 % The case that the feeding bat is also an intruder
            else % The case the intruded bat comes back
                cnt_intr = cnt_intr + 1;
                b_intruder = bb;
                t_intr_first = find(is_intruded,1,'first');
                t_intr_last = find(is_intruded,1,'last');
            end
        end
    end
    if cnt_intr > 1 % if there are multiple intruders, ignore the case
        Rwr.intruded(tt) = true;
        Rwr.mlt_intr(tt) = true;
    elseif cnt_intr == 1 % Consider the cases when there is only one intruder
        Rwr.b_intruder(tt) = b_intruder;
        Rwr.intruder_id(tt) = bat_id(b_intruder);
        Rwr.intruded(tt) = true;
        if t_intr_last < length(is_intruded)
            Rwr.feeder_win(tt) = true;
        else
            Rwr.feeder_lose(tt) = true;
        end
        Rwr.t_intr_come(tt) = t(find(t==Rwr.t_r(tt))+t_intr_first);
        Rwr.t_intr_leave(tt) = t(find(t==Rwr.t_r(tt))+t_intr_last);
    end
end

fprintf('Total counts of rewards: %d\n', size(Rwr,1))
fprintf('Feeding bat was intruded: %d/%d\n',sum(Rwr.intruded),size(Rwr,1))
fprintf('Feeding bat was intruded by multiple intruders: %d/%d\n',sum(Rwr.mlt_intr),sum(Rwr.intruded))
fprintf('Feeding bat wins: %d/%d\n',sum(Rwr.feeder_win & ~Rwr.mlt_intr),sum(Rwr.intruded))
fprintf('Feeding bat loses: %d/%d\n',sum(Rwr.feeder_lose & ~Rwr.mlt_intr),sum(Rwr.intruded))

%% Social hierarchy
groupcounts(Rwr(Rwr.intruded==true & Rwr.mlt_intr==false,:),["feeding_id","intruder_id" "feeder_win","feeder_lose"])

%% Save file
save(fullfile(BehFile.folder,'Feeder_Intrusion.mat'),'Rwr')

%% Visualization of experimental features
% %===FIGURE: Reward depletion
% figure;
% % x = find(diff(AnalogSignals(:,1))>2)/AnalogFrameRate/60; % min
% y = pulsewidth(AnalogSignals(:,1))/AnalogFrameRate; % reward duration (
% plot(y,'k')
% title('Reward depletion')
% xlabel('Reward')
% ylabel('Reward duration (sec)')
% 
% %===FIGURE: feeding interval
% figure;
% edges = 0:5:ceil(max(diff(t_rwrd))/10)*10;
% histogram(diff(t_rwrd),edges)
% title('Feeding interval')
% xlabel('Interval (sec)')
% ylabel('Frequency')
