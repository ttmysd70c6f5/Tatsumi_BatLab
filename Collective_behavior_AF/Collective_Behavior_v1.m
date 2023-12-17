function Collective_Behavior()
% This is a script for analyzing the collective behavior data of bats by Angelo
% Dataset: Extracted_Behavior_(date).mat

%% Directory
currdir = split(cd,'\');
behdir = fullfile('C:\Users\Tatsumi\Documents\Behavior Analysis\Data_Tatsumi\Behavior',currdir{end});
ephydir = fullfile('C:\Users\Tatsumi\Documents\Behavior Analysis\Data_Tatsumi\Ephys\Ephys_Restrictive_Raw',currdir{end});

sessions = dir(fullfile(behdir,'**/Extracted_Behavior*.mat')); % extract file paths
TTLs = dir(fullfile(ephydir,'**/Ext_Ephys_1/*c3d_1-Bat_Cluster.mat')); % reward signals


%% Concatenate all sessions
% Concatenate all sessions
r_all = [];
bflying_all = [];
btakeoff_all = [];
t_all = [];
exp_d = []; % Experiment day

for i=1:length(sessions)
    load(fullfile(sessions(i).folder, sessions(i).name)) % load each session
    
    d_session = repmat(i,[T 1]); % the day of session
    [btakeoff,~] = findchanges(bflying,0,1); % Takeoff timing (logical)

    r_all = [r_all; r]; % concatenate all sessions
    bflying_all = [bflying_all; bflying]; % concatenate all sessions
    btakeoff_all = [btakeoff_all; btakeoff];
    t_all = [t_all; t];
    exp_d = [exp_d; d_session];
end

T_all = size(r_all,1); % Length of all sessions

save(fullfile(cd,'Extracted_All.mat'),...
    'bflying_all','btakeoff_all','t_all','T_all','exp_d',... % concatenated data
    'bat_clr','bat_nms','bat_pair_nms','bat_pairs','Fs','Group_name','n_tags') % common parameters
% save(fullfile(cd,'Extracted_All.mat'),'r','bflying_all','btakeoff_all','t_all','T_all','exp_d')

%% Load data
load(fullfile(cd,data.name))
load(fullfile(cd,'bpos.mat'),'bpos')
fprintf('Current session is %s \n', data.folder)

%% Check the spatial distribution of clusters
r_all = reshape(permute(r,[1 3 2]),[T*n_tags 3]); % concatenated positions for all bats
bpos_all = reshape(bpos,[T*n_tags 1]);

%=== FIGURE: 3D visualization of landing spot clusters
figure(); hold on
box on; ax = gca; ax.BoxStyle = 'full';
xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:)); xlabel("x"); ylabel("y"); zlabel("z"); view([135 50]);
for j=1:length(unique(bpos))-2
    scatter3(r_all(bpos==j,1),r_all(bpos_all==j,2),r_all(bpos_all==j,3),'.','MarkerEdgeAlpha',0.65)
end
scatter3(r_all(bpos_all==-1,1),r_all(bpos_all==-1,2),r_all(bpos_all==-1,3),'k.','MarkerEdgeAlpha',0.65);
hold off;

%% Ad-hoc correction of cluster identification
% Clutser identification was checked manually after DBSCAN-knnsearch.
% Reallocate cluster labels and index
% (-1=outliers,0=flying,1~100=home,101~200=feeder/banana)

if contains(extracted_CDPfile.folder,'Dataset_1_3')
elseif contains(extracted_CDPfile.folder,'Dataset_1_4')
elseif contains(extracted_CDPfile.folder,'Dataset_3_1')
elseif contains(extracted_CDPfile.folder,'Dataset_5_2')
    Avalnew = [-1 0 1 2 3 4 107 5 6];
    poslabel = ["Home 1" "Home 2" "Home 3" "Home 4" "Home 5" "Home 6" "Banana"];
    [~, ~, indAval] = unique(bpos);
    Anew = Avalnew(indAval);
    bposnew = reshape(Anew, size(bpos));
    n_pos = length(unique(bposnew)); % number of unique cluster index
    posval = unique(bposnew); % cluster index of positions
elseif contains(extracted_CDPfile.folder,'Dataset_5_3')
    Avalnew = [-1 0 1 2 3 104 105 106 107 1];
    poslabel = ["Home 1" "Home 2" "Home 3" "Feeder 1" "Feeder 2" "Feeder 3" "Feeder 4"];
    [~, ~, indAval] = unique(bpos);
    Anew = Avalnew(indAval);
    bposnew = reshape(Anew, size(bpos));
    n_pos = length(unique(bposnew)); % number of unique cluster index
    posval = unique(bposnew); % cluster index of positions
elseif contains(extracted_CDPfile.folder,'Dataset_5_4')
end

clear bpos
%% Visualizing transitions of crowdedness
tcrowd = zeros(T,length(unique(bposnew)));
for i=1:n_pos
   tcrowd(:,i) = sum(bposnew(:,:)==posval(i),2)/n_tags; % time course of crowdedness for each position (max = 1.0)
end

%=== FIGURE: timecourse of crowdedness
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.1 0.2 0.8 0.6]); hold on;
title('Crowdedness of each position'); xlabel('Time (sec)'); ylim([-0.5,n_pos+0.5]); grid
ynms = cell(1,n_pos);
for i=1:n_pos
    if posval(i) == -1
        ynms{1,i} = 'Outliers';
    elseif posval(i) == 0
        ynms{1,i} = 'Flying';
    else
        ynms{1,i} = poslabel(i-2);
    end
    plot(t,tcrowd(:,i)+i-1)
end
yticks(0:n_pos-1); yticklabels(ynms)
hold off

%% State transition
%% Count transitions of takeoff -> landing state
statenames = ["social-home" "isolated-home" "social-feeder" "isolated-feeder"];
num_state = length(statenames);

[P,state_cnt] = comp_transmat(bposnew,num_state,n_tags,[-1 0],'state',false);
P = P./sum(P,2,'omitnan');

%=== FIGURE: visualization of state transition without p-values
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.05 0.1 0.9 0.8]);
tl = tiledlayout(2, ceil(n_tags/2), 'TileSpacing', 'tight');
title(tl,'State transition graph')

xnodepos = [-0.65 0.95 -0.99 0.65];
ynodepos = [0.65 0.95 -0.99 -0.6];
mk_sz = state_cnt; mk_sz = 24.*mk_sz./max(mk_sz(:)); mk_sz = mk_sz+1;
for i=1:n_tags
%     mc = dtmc(P(any(P(:,:,i),2),any(P(:,:,i),1),i),'StateNames',statenames(any(P(:,:,i),2)));
    mc = dtmc(P(:,:,i),'StateNames',statenames);
    nexttile(); H = graphplot(mc,'ColorEdges',true);
%     H.XData = xnodepos(any(P(:,:,i),1)); H.YData = ynodepos(any(P(:,:,i),1));
    H.XData = xnodepos; H.YData = ynodepos; H.MarkerSize = mk_sz(:,i);
    title(bat_nms(i,:),'FontSize',12,'Color',bat_clr(i,:))
end

%% Randomization test
if ~exist(fullfile(cd,'P_randtest.mat'),'file')
    
    n_randtest = 10000;
    P_randtest = zeros(length(statenames),length(statenames),n_tags,n_randtest);
    statenames = ["social-home" "isolated-home" "social-feeder" "isolated-feeder"];
    num_state = length(statenames);

    for i=1:n_randtest
        P_randtest(:,:,:,i) = comp_transmat(bposnew,num_state,n_tags,[-1 0],'state',true);
    end
    P_randtest = P_randtest./sum(P_randtest,2,'omitnan');

    %=== SAVE: randomized transition matrice
    save(fullfile(cd,'P_randtest.mat'),'P_randtest')
else
    load(fullfile(cd,'P_randtest.mat'),'P_randtest')
    n_randtest = size(P_randtest,4);
end

% Compute p-value
p_val = zeros(num_state,num_state,n_tags);
for i=1:n_tags
    figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    tl = tiledlayout(num_state, num_state, 'TileSpacing', 'tight'); title(tl,sprintf('Randomization test: %s',bat_nms(i,:)))
    xlabel(tl,'Probability'); ylabel(tl,'Count')
    for j=1:num_state
        for k=1:num_state
            p_mean = mean(P_randtest(j,k,i,:),'omitnan');
            p_obs = P(j,k,i); % probabilities of observed matrix
            p_gen = P_randtest(j,k,i,:); % probabilities of generated transition matrix
            p_val(j,k,i) = sum(abs(p_mean-p_gen) > abs(p_mean - p_obs)) / n_randtest;

            %=== FIGURE: histogram of replicated probabilities
            nexttile();
            histogram(p_gen,(0:0.02:1)); hold on;
            plot(repmat(p_mean,[1,2]),ylim,'b--');
            plot(repmat(p_obs,[1,2]),ylim,'r--');
            plot(repmat(2*p_mean-p_obs,[1,2]),ylim,'k--');
            xlim([0,1]); xticks(0:0.2:1.0);
            title(sprintf('%s \\rightarrow %s',statenames{j},statenames{k}),'Interpreter','tex')
            legend('',sprintf('Mean: %.3f', p_mean),sprintf('Obs.: %.3f',p_obs)); hold off
        end
    end

    %=== Save: histogram of randomization test
    % saveas(gcf, fullfile(cd,'Figs_Tatsumi',sprintf('random_test_%d.png',i)))
end
    
%=== FIGURE: visualization of state transition with p-values
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.05 0.1 0.9 0.8]);
tl = tiledlayout(2, ceil(n_tags/2), 'TileSpacing', 'tight');
title(tl,'State transition graph')

xnodepos = [-0.65 0.95 -0.99 0.65];
ynodepos = [0.65 0.95 -0.99 -0.6];
mk_sz = state_cnt; mk_sz = 24.*mk_sz./max(mk_sz(:)); mk_sz = mk_sz+1;
for i=1:n_tags
%     mc = dtmc(P(any(P(:,:,i),2),any(P(:,:,i),1),i),'StateNames',statenames(any(P(:,:,i),2)));
    mc = dtmc(P(:,:,i),'StateNames',statenames);
    pvalnow = p_val(:,:,i);
    
    nexttile();
    H = graphplot(mc,'ColorEdges',true);
    H.XData = xnodepos; H.YData = ynodepos; H.MarkerSize = mk_sz(:,i);
%     H.XData = xnodepos(any(P(:,:,i),1)); H.YData = ynodepos(any(P(:,:,i),1));   
    H.EdgeCData(find(isnan(P(:,:,i)'))) = 0; % set a value for elements with 0 in P
    [prow,pcol] = find(isnan(P(:,:,i))); highlight(H,prow,pcol,'LineStyle',':')
    [prow,pcol] = find(pvalnow < 0.05 & pvalnow >= 0.01); highlight(H,prow,pcol,'LineWidth',3); labeledge(H,prow,pcol,'*')
    [prow,pcol] = find(pvalnow < 0.01); highlight(H,prow,pcol,'LineWidth',3); labeledge(H,prow,pcol,'**')
    
    H.EdgeFontSize = 12;
    title(bat_nms(i,:),'FontSize',12,'Color',bat_clr(i,:))
end

%% Spots transitions
%% Count transitions of takeoff -> landing spots
spotnames = ["Home 1","Home 2","Home 3", "Feeder 1", "Feeder 2", "Feeder 3", "Feeder 4"];
n_spot = length(unique(bposnew))-2;

[P,spot_cnt] = comp_transmat(bposnew,n_spot,n_tags,[-1 0],'spot',false);
P = P./sum(P,2,'omitnan');

%% Inter-flight interval (IFI)

[~,btakeoff] = findchanges(bposnew,[-1 0],1);
IFI = cell(1,n_tags);
for i=1:n_tags
    IFI{1,i} = diff(find(btakeoff(:,i)));
end

%=== FIGURE: IFI plot
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.05 0.2 0.9 0.35]);
tl = tiledlayout(1, n_tags, 'TileSpacing', 'tight');
title(tl,'Inter-flight interval'); xlabel(tl,'Time (min)'); ylabel(tl,'Count')
hedge = 0:1000:30000;
for i=1:n_tags
   nexttile(); histogram(IFI{1,i},hedge); title(bat_nms(i,:),'FontSize',12,'Color',bat_clr(i,:))
   xticks(0:6000:30000); xticklabels(xticks/6000); xlim([0,30000]); ylim([0,50])
   hold on; plot(repmat(median(IFI{1,i}),[1,2]),ylim,'k--'); hold off
   legend('',sprintf('Median: %.0fsec',median(IFI{1,i})/100));
end


%% Flight trajectory
% Extract flights by spline interpolation
[btakeoff,blanding] = findchanges(bflying,0,1);

n_smp = 10; % number of downsampled points
fl = [];

for i=1:n_tags
    if sum(btakeoff(:,i))==sum(blanding(:,i))
        fl_idx = zeros(sum(btakeoff(:,i)),2); % index of takeoffs and landings
        fl_idx(:,1) = find(btakeoff(:,i));
        fl_idx(:,2) = find(blanding(:,i));
        
        fl_now = zeros(size(fl_idx,1),n_smp*3); % downsampled trajectories(n_smp x 3D dimensionality)
        
        for j=1:size(fl_idx,1)
           r_flight = r(fl_idx(j,1):fl_idx(j,2),:,i); 
           r_intp = zeros(n_smp,3); % positions downsampled by spline interpolation
           for d=1:3
               r_intp(:,d) = spline(1:length(r_flight(:,d)),r_flight(:,d),linspace(1,length(r_flight(:,d)),10)); % spline interpolation
           end
           r_intp = reshape(r_intp,[1,n_smp*3]);
           fl_now(j,:) = r_intp;
        end
        
        fl = [fl;fl_now]; % sample x variable
        
    else
        error('Number of takeoffs and landings does not match')
    end
    
end

%% Hierarchical clustering
hc_tree = linkage(fl,'single','euclidean'); % hierarchical cluster tree
figure(); dendrogram(hc_tree,0)

hc_clus = cluster(hc_tree,'cutoff',1.2,'Criterion','distance');
[clus_cts,clus_id] = groupcounts(hc_clus);
clus_good = clus_id(clus_cts > 5);
length(unique(clus_good))

figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
tl = tiledlayout(3, 6, 'TileSpacing', 'tight');
for i=1:length(clus_good)
    fl_clus = fl(hc_clus==clus_good(i),:);
    nexttile(); hold on
    for j=1:size(fl_clus,1)
        xyz = reshape(fl_clus(j,:),[10,3]);
        
%         scatter3(spline(1:10,xyz(:,1),1:.1:10), spline(1:10,xyz(:,2),1:.1:10),spline(1:10,xyz(:,3),1:.1:10),1,repmat(clus_good(i),[91,1]));
        plot3(spline(1:10,xyz(:,1),1:.1:10), spline(1:10,xyz(:,2),1:.1:10),spline(1:10,xyz(:,3),1:.1:10),'-','Color','b');
        view([135 45])
    end
    hold off
end
% Segregate flights



%% Data visualization
%=== FIGURE: Trajectories individual bats
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.2 0.25 0.7 0.35]);
tiledlayout(1, n_tags, 'TileSpacing', 'none');
for i=1:n_tags
    nexttile;
    plot3(r(:,1,i),r(:,2,i),r(:,3,i),'-','Color', bat_clr(i,:)); 
    xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:)); title(bat_nms(i,:)); view(90,90);
    axis equal;
end

%% Data visualization
%=== FIGURE: Transition of trajectories individual bats
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
for i=1:n_tags
    subplot(5,n_tags,i);
    plot3(r(1:floor(T/5),1,i),r(1:floor(T/5),2,i),r(1:floor(T/5),3,i),'-','Color', bat_clr(i,:)); 
    xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:)); title(bat_nms(i,:)); view(90,90);
    
    subplot(5,n_tags,i+n_tags);
    plot3(r(floor(T/5)+1:floor(T/5)*2,1,i),r(floor(T/5)+1:floor(T/5)*2,2,i),r(floor(T/5)+1:floor(T/5)*2,3,i),'-','Color', bat_clr(i,:)); 
    xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:)); title(bat_nms(i,:)); view(90,90);
    
    subplot(5,n_tags,i+n_tags*2); 
    plot3(r(floor(T/5)*2+1:floor(T/5)*3,1,i),r(floor(T/5)*2+1:floor(T/5)*3,2,i),r(floor(T/5)*2+1:floor(T/5)*3,3,i),'-','Color', bat_clr(i,:)); 
    xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:)); title(bat_nms(i,:)); view(90,90);
    
    subplot(5,n_tags,i+n_tags*3); 
    plot3(r(floor(T/5)*3+1:floor(T/5)*4,1,i),r(floor(T/5)*3+1:floor(T/5)*4,2,i),r(floor(T/5)*3+1:floor(T/5)*4,3,i),'-','Color', bat_clr(i,:)); 
    xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:)); title(bat_nms(i,:)); view(90,90);
    
    subplot(5,n_tags,i+n_tags*4); 
    plot3(r(floor(T/5)*4+1:end,1,i),r(floor(T/5)*4+1:end,2,i),r(floor(T/5)*4+1:end,3,i),'-','Color', bat_clr(i,:)); 
    xlim(r_lim(1,:)); ylim(r_lim(2,:)); zlim(r_lim(3,:)); title(bat_nms(i,:)); view(90,90);
end

%% Data visualization
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

%% Data visualization
%=== FIGURE: Timecourse of bats pair distance
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.2 0.2 0.7 0.7]);
tiledlayout(2, size(bat_pairs,1)/2, 'TileSpacing', 'none');
for i=1:size(bat_pairs,1)
    nexttile; plot3(d(:,1,i),d(:,2,i),d(:,3,i),'-'); 
    xlim(r_lim(1,:)*2); ylim(r_lim(2,:)*2); title(bat_pair_nms(i,:));view(90,90);
    axis equal;
end

%%
figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.2 0.25 0.7 0.35]);
tiledlayout(1, n_tags, 'TileSpacing', 'none');
for i=1:n_tags
    nexttile;
    scatter3(r(logical(1-bflying(:,i)),1,i),r(logical(1-bflying(:,i)),2,i),r(logical(1-bflying(:,i)),3,i),'.');
    view(-75,15)
%     axis equal
end
