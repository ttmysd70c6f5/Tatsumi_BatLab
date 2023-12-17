function state_transitions_sessions()

sessions = dir('**/Extracted_Behavior*.mat'); % extract file paths

%% Make figs
load(fullfile(sessions(1).folder, sessions(1).name),'n_tags','bat_nms')
for i = 1:n_tags
    %=== FIGURE: visualization of state transition with p-values
    figs(i).f = figure();
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]); tl = tiledlayout(3, ceil(length(sessions)/3), 'TileSpacing', 'tight');
    title(tl,sprintf('State transition graph: bat %s',bat_nms(i,:)))
end

%%
for ss = 1:length(sessions)

%% Load data
load(fullfile(sessions(ss).folder, sessions(ss).name))
load(fullfile(sessions(ss).folder, 'bpos.mat'))
fprintf('Current session is %s \n', sessions(ss).folder)

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

%% Count transitions of takeoff -> landing state
statenames = ["social-home" "isolated-home" "social-feeder" "isolated-feeder"];
num_state = length(statenames);

% compute state transitions
[P,state_cnt] = comp_transmat(bposnew,num_state,n_tags,[-1 0],'state',false);
P = P./sum(P,2,'omitnan');

% Randomization test
if ~exist(fullfile(sessions(ss).folder,'P_randtest.mat'),'file')
    
    n_randtest = 10000;
    P_randtest = zeros(length(statenames),length(statenames),n_tags,n_randtest);
    statenames = ["social-home" "isolated-home" "social-feeder" "isolated-feeder"];
    num_state = length(statenames);

    for i=1:n_randtest
        P_randtest(:,:,:,i) = comp_transmat(bposnew,num_state,n_tags,[-1 0],'state',true);
    end
    P_randtest = P_randtest./sum(P_randtest,2,'omitnan');

    %=== SAVE: randomized transition matrice
    save(fullfile(sessions(ss).folder,'P_randtest.mat'),'P_randtest')
else
    load(fullfile(sessions(ss).folder,'P_randtest.mat'),'P_randtest')
    n_randtest = size(P_randtest,4);
end

% Compute p-values
p_val = zeros(num_state,num_state,n_tags);
for i=1:n_tags
    for j=1:num_state
        for k=1:num_state
            p_mean = mean(P_randtest(j,k,i,:),'omitnan');
            p_obs = P(j,k,i); % probabilities of observed matrix
            p_gen = P_randtest(j,k,i,:); % probabilities of generated transition matrix
            p_val(j,k,i) = sum(abs(p_mean-p_gen) > abs(p_mean - p_obs)) / n_randtest;
        end
    end
end
%%
%=== FIGURE
xnodepos = [-0.65 0.95 -0.99 0.65];
ynodepos = [0.65 0.95 -0.99 -0.6];
mk_sz = state_cnt; mk_sz = 24.*mk_sz./max(mk_sz(:)); mk_sz = mk_sz+1;
for i=1:n_tags
    figure(figs(i).f)
%     ax = subplot(3, ceil(length(sessions)/3), ss);
    nexttile()  
    mc = dtmc(P(:,:,i),'StateNames',statenames); % compute markov matrix
    pvalnow = p_val(:,:,i);
    
    % plot
    H = graphplot(mc,'ColorEdges',true);
    H.XData = xnodepos; H.YData = ynodepos; H.MarkerSize = mk_sz(:,i);
    H.EdgeCData(find(isnan(P(:,:,i)'))) = 0; % set a value for elements with 0 in P
    [prow,pcol] = find(isnan(P(:,:,i))); highlight(H,prow,pcol,'LineStyle',':')
    [prow,pcol] = find(pvalnow < 0.05 & pvalnow >= 0.01); highlight(H,prow,pcol,'LineWidth',3); labeledge(H,prow,pcol,'*')
    [prow,pcol] = find(pvalnow < 0.01); highlight(H,prow,pcol,'LineWidth',3); labeledge(H,prow,pcol,'**')
    H.EdgeFontSize = 12;
%     title(bat_nms(i,:),'FontSize',12,'Color',bat_clr(i,:))
end



end


end
