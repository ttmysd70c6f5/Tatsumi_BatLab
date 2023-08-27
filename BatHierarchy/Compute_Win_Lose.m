function p = Compute_Win_Lose(RwrFile)
% Input: the list of directories to concatenate and compute win-lose rate

%% 1. Concatenate file
% RwrFile = dir('**\Ext_Behavior*\Feeder_Intrusion.mat');
Rwr_cnct = [];
for i = 1:length(RwrFile)
    load(fullfile(RwrFile(i).folder,RwrFile(i).name),'Rwr');
    Rwr_cnct = [Rwr_cnct; Rwr];
end

id_list = unique(Rwr_cnct.feeding_id);
id_pair = id_list(nchoosek([1:length(id_list)],2));

%% Social hierarchy
groupcounts(Rwr_cnct(Rwr_cnct.intruded==true & Rwr_cnct.mlt_intr==false,:),["feeding_id","intruder_id" "feeder_win","feeder_lose"])

%% Visualize win-lose counts
win_lose_cnt = zeros(size(id_pair,1),2,3);
for i=1:size(id_pair,1)
    win_lose_cnt(i,1,2) = sum(Rwr_cnct.feeder_win(Rwr_cnct.feeding_id==id_pair(i,1) & Rwr_cnct.intruder_id==id_pair(i,2)));
    win_lose_cnt(i,2,2) = sum(Rwr_cnct.feeder_lose(Rwr_cnct.feeding_id==id_pair(i,1) & Rwr_cnct.intruder_id==id_pair(i,2)));

    win_lose_cnt(i,1,3) = sum(Rwr_cnct.feeder_lose(Rwr_cnct.feeding_id==id_pair(i,2) & Rwr_cnct.intruder_id==id_pair(i,1)));
    win_lose_cnt(i,2,3) = sum(Rwr_cnct.feeder_win(Rwr_cnct.feeding_id==id_pair(i,2) & Rwr_cnct.intruder_id==id_pair(i,1)));
end

win_lose_cnt(:,:,1) = win_lose_cnt(:,:,2) + win_lose_cnt(:,:,3); % total win-lose counts

figure; set(gcf, 'units', 'normalized', 'outerposition', [0.2 0.3 0.6 0.4]);
tiledlayout(1, 3, 'TileSpacing', 'tight');
for i = 1:3
    nexttile
    imagesc(win_lose_cnt(:,:,i)./sum(win_lose_cnt(:,:,i),2))
    set(gca,'FontSize',12)
    xticks([1,2])
    yticks(1:size(id_pair,1))
    xtickangle(0)
    colorbar
    colormap jet
    ylbls = cell(size(id_pair,1),1);
    if i == 1
        for j = 1:size(id_pair,1)
            ylbls{j,1} = [num2str(id_pair(j,1)) 'vs' num2str(id_pair(j,2))];
        end
        xticklabels({'Left win','Left lose'})
        yticklabels(ylbls)
        title('Total')
    elseif i == 2
        for j = 1:size(id_pair,1)
            ylbls{j,1} = [num2str(id_pair(j,1)) '->' num2str(id_pair(j,2))];
        end
        xticklabels({'Def. win','Def. lose'})
        yticklabels(ylbls)
        title('Defender <- Intruder')
    elseif i == 3
        for j = 1:size(id_pair,1)
            ylbls{j,1} = [num2str(id_pair(j,1)) '<-' num2str(id_pair(j,2))];
        end
        xticklabels({'Int. win','Int. lose'})
        yticklabels(ylbls)
        title('Intruder -> Defender')
    end
end

%% chi squared test
win_lose_cnt_2 = zeros(2,2,size(win_lose_cnt,1));
for i = 1:size(win_lose_cnt_2,3)
    win_lose_cnt_2(:,:,i) = [win_lose_cnt(i,:,2); win_lose_cnt(i,end:-1:1,3)];
end

chi_test_tbl_est = zeros(2,2,size(win_lose_cnt,1)); % estimated counts
for i = 1:size(win_lose_cnt_2,3)
    for j = 1:2
        for k = 1:2
            chi_test_tbl_est(j,k,i) = sum(win_lose_cnt_2(:,k,i)) * sum(win_lose_cnt(j,:,i) / sum(sum(win_lose_cnt(:,:,i))));
        end
    end
end

p = zeros(size(win_lose_cnt_2,3),1);
for i = 1:length(p)
    chi = sum(sum((win_lose_cnt_2(:,:,i) - chi_test_tbl_est(:,:,i)).^2./chi_test_tbl_est(:,:,i)));
    deg_free = (size(chi_test_tbl_est(:,:,i),1)-1)*(size(chi_test_tbl_est(:,:,i),2)-1);
    p(i) = 1-chi2cdf(chi,deg_free);
end

% %% chi squared test
% chi_test_tbl = zeros(size(id_pair,1),2,size(win_lose_cnt,3));
% for i = 1:size(chi_test_tbl,1)
%     for j = 1:size(chi_test_tbl,2)
%         for k = 1:size(chi_test_tbl,3)
%             chi_test_tbl(i,j,k) = sum(win_lose_cnt(:,j,k)) * sum(win_lose_cnt(i,:,k) / sum(sum(win_lose_cnt(:,:,k))));
%         end
%     end
% end
% 
% p = zeros(3,1);
% for k = 1:size(chi_test_tbl,3)
%     chi = sum(sum((win_lose_cnt(:,:,k) - chi_test_tbl(:,:,k)).^2./chi_test_tbl(:,:,k)));
%     deg_free = (size(chi_test_tbl(:,:,k),1)-1)*(size(chi_test_tbl(:,:,k),2)-1);
%     p(k) = 1-chi2cdf(chi,deg_free);
% end
