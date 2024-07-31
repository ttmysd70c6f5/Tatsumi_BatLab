function PlotOtherFlightModulation(BhvData,SocialFlights,spikeData,target,FigDir)
% Common parameters
t_margin = [-1 1];
% fr_th = 0.5; % firing rate exclusion criteria (Hz)
fs_bhv = BhvData.fs;
n_bats = BhvData.n_tags;

fr_all = {};
% bbb = 8;
for bbb = 1:n_bats
    bname = SocialFlights(bbb).bname;
    
    clu_idx = SocialFlights(bbb).clu_idx;
    clu_list = SocialFlights(bbb).clu_list;
    n_clu = length(clu_list);
    for clu = 1:n_clu
        % clu = 1;
        clu_id = clu_list(clu);
        
        
        is_social = SocialFlights(bbb).social(clu_idx==clu_id);
        a_obs = SocialFlights(bbb).a_abs_obs(clu_idx==clu_id);
        
        % 1. Firing rate curve
        % Parameters
        w_step = 0.1; % sliding window (sec)
        w_step_2 = 0.01;
        w_sz = 0.5; % bin size (sec)
        pts_margin = ceil(w_sz/w_step); % marginal points (remove later)
        t_bin = t_margin(1)-pts_margin*w_step:w_step:t_margin(2)+pts_margin*w_step;
        t_edge = [t_bin(1)-w_step/2, t_bin+w_step/2];
        % t_bin_2 = t_margin(1)-pts_margin*1/w_step_2:1/w_step_2:t_margin(2)+pts_margin*1/w_step_2;
        % t_edge_2 = [t_bin_2(1)-1/w_step_2/2, t_bin_2+1/w_step_2/2];
        t_bin_2 = t_margin(1)-pts_margin*w_step_2:w_step_2:t_margin(2)+pts_margin*w_step_2;
        t_edge_2 = [t_bin_2(1)-w_step_2/2, t_bin_2+w_step_2/2];
        
        % Extract event-aligned spike activity
        % target_units =  find(spikeData.good_units.fr >= fr_th); % Exclusion criteria
        % n_unit = length(target_units);
        n_unit = spikeData.numGoodUnits;
        
        n1 = sum(is_social); n2 = sum(~is_social);
        sp_mat = NaN(n1+n2,length(t_bin),n_unit); % event aligned spikes
        sp_mat_2 = NaN(n1+n2,length(t_bin_2),n_unit); % event aligned spikes
    
        if n1 >= 5 && n2 >= 5 % Check exclusion criteria
            t_flight = SocialFlights(bbb).t_sec(clu_idx==clu_id);
            for unit = 1:n_unit
                sptimes_unit_sec = spikeData.good_units.spikeTimes_usec{unit} / 1e6; % spike times of the current unit in sec
                for ff = 1:n1+n2
                    sp_mat(ff,:,unit) = histcounts(sptimes_unit_sec,t_flight(ff)+t_edge); % [flight,time,unit]
                end
            end
        
            % Poisson c-test
            n_sp_1 =  squeeze(sum(movsum(sp_mat(is_social,:,:),w_sz/w_step,2),1));
            n_sp_2 =  squeeze(sum(movsum(sp_mat(~is_social,:,:),w_sz/w_step,2),1)); % # of spikes for each bin [time,unit]
            ps_ctest = NaN(length(t_bin)-2*pts_margin,n_unit); % [time,unit]
            p_ctest = NaN(n_unit,1); % [unit]
            for unit = 1:n_unit
                for i = 1:length(t_bin)-2*pts_margin
                    ps_ctest(i,unit) = poisson_c_test(n_sp_1(pts_margin+i,unit),n_sp_2(pts_margin+i,unit),n1,n2) * n_unit; % p value with Bonferroni correction
                end
            end
        
            % Identify the optimal bin based on c-test results
            opt_bins = zeros(n_unit,1); % optimal bin [unit,1]
            for unit = 1:n_unit
                [~,opt_bins(unit)] = min(ps_ctest(:,unit));
                p_ctest(unit) = ps_ctest(opt_bins(unit),unit);
            end
            t_opt_bins = t_bin(opt_bins+pts_margin);
        
            % Firing rate
            fr1 = n_sp_1(1+pts_margin:end-pts_margin,:)/w_sz/n1;
            fr2 = n_sp_2(1+pts_margin:end-pts_margin,:)/w_sz/n2; % bin firing rate [time,unit]
            sem1 = squeeze(std(movsum(sp_mat(is_social,:,:),w_sz/w_step,2)/w_sz,0,1)) / sqrt(n1); % SEM for social flights [time,unit]
            sem2 = squeeze(std(movsum(sp_mat(~is_social,:,:),w_sz/w_step,2)/w_sz,0,1)) / sqrt(n2); % SEM for social flights [time,unit]
            sem1 = sem1(1+pts_margin:end-pts_margin,:); sem2 = sem2(1+pts_margin:end-pts_margin,:);
        
            % ANOVA test
            p_anova = NaN(n_unit,1); % [unit]
            for unit = 1:n_unit
                p_anova(unit) = anova1([fr1(ceil(size(fr1,1)):end,unit),fr2(ceil(size(fr2,1)):end,unit)],{'social','non-social'},'off') * n_unit; % p value with Bonferroni correction
            end
        
            % GLM for each neuron
            for unit = 1:n_unit
                n_sp_glm = movsum(sp_mat(:,:,unit),w_sz/w_step,2); n_sp_glm = n_sp_glm(:,opt_bins(unit)+pts_margin); % spikes [F,T]
                social_glm = double(is_social);
                a_glm = a_obs;
            
                T_glm = table(a_glm,social_glm,n_sp_glm);
                mdl = stepwiseglm(T_glm,'constant','Distribution','poisson','Link','log');
                out_glm(unit).coef = mdl.CoefficientNames;
                out_glm(unit).intercept = mdl.Coefficients.pValue(strcmp(mdl.CoefficientNames,'(Intercept)'));
                out_glm(unit).social = mdl.Coefficients.pValue(strcmp(mdl.CoefficientNames,'social_glm'));
                out_glm(unit).a_abs = mdl.Coefficients.pValue(strcmp(mdl.CoefficientNames,'a_glm'));
                out_glm(unit).p = mdl.Coefficients.pValue;
            end
        
            % Spike arrays for raster plots
            sp_mat_2 = NaN(n1+n2,length(t_bin_2),n_unit);
            for unit = 1:n_unit
                sptimes_unit_sec = spikeData.good_units.spikeTimes_usec{unit} / 1e6; % spike times of the current unit in sec
                for ff = 1:n1+n2
                    sp_mat_2(ff,:,unit) = histcounts(sptimes_unit_sec,t_flight(ff)+t_edge_2); % [flight,time,unit]
                end
            end
            
            % Peri-event firing rate (social vs non-social)
            pefr = NaN(n_unit,2); % firing rate around events (Hz) [unit, social/asocial]
            for unit = 1:n_unit
                pefr(unit,1) = sum(sp_mat_2(is_social,1+pts_margin:end-pts_margin,unit),'all')/(n1*diff(t_margin));
                pefr(unit,2) = sum(sp_mat_2(~is_social,1+pts_margin:end-pts_margin,unit),'all')/(n2*diff(t_margin));
            end
        
            % Remove margin
            sp_mat = sp_mat(:,1+pts_margin:end-pts_margin,:);
            sp_mat_2 = sp_mat_2(:,1+pts_margin:end-pts_margin,:);

            % Output structure
            for unit = 1:n_unit
                unit_id = spikeData.good_units.cluster_id(unit);
                fr = spikeData.good_units.fr(unit);
                fr_all = [fr_all; {bname,unit_id,fr,clu_id,p_ctest(unit),p_anova(unit),out_glm(unit).social,fr1(:,unit),fr2(:,unit)}];
            end
            
        
        
        
            %%% FIGURE: Firing rate curve and raster plot
            for unit = 1:n_unit
            % unit = 35;
                unit_id = spikeData.good_units.cluster_id(unit);
                fr = spikeData.good_units.fr(unit);
                depth = spikeData.good_units.depth(unit);
                region = spikeData.good_units.region{unit};
                %%% Initialize figure bundle
                fig = figure('visible','off');
                % fig = figure;
                set(fig,'units','normalize','Position',[0.2 0.1 0.3 0.8])
                % tl = tiledlayout(2,1,'TileSpacing','tight');
                title_str = sprintf(['%s_%d\n' ...
                    '%s,%dum, %.1fHz(social:%.1fHz, Non-social:%.1fHz)\n' ...
                    'C-test:%.3f,ANOVA:%.3f,GLM(social):%.3f,GLM(acc):%/3f'], ...
                    bname,unit_id, ...
                    region,depth,fr,pefr(unit,1),pefr(unit,2), ...
                    p_ctest(unit),p_anova(unit),out_glm(unit).social,out_glm(unit).a_abs);
                % title(title_str,'Interpreter','None','FontSize',14)
                % title_str = sprintf('Unit%d, %s, %dum, %.1fHz (Social: %.1fHz, Non-social: %.1fHz)',unit_id,region,depth,fr,pefr,pefr_asocial);
                % title(tl,title_str,'Interpreter','None','FontSize',14)
            
                % Plot 1: Firing rate curve
                % nexttile
                subplot('Position',[0.1 0.52 0.8 0.40])
                t_fig = t_margin(1):w_step:t_margin(2); % time vector for figure
               
                hold on
                % plot(t_fig,e_test_poisson(1,:)/w_sz/n1,'r','DisplayName','Social')
                % plot(t_fig,e_test_poisson(2,:)/w_sz/n2,'k','DisplayName','Non-social')
                boundedline(t_fig,fr1(:,unit),sem1(:,unit),'color',[250 50 50]/255,'alpha','transparency', 0.1)
                boundedline(t_fig,fr2(:,unit),sem2(:,unit),'color',[100 100 100]/255,'alpha','transparency', 0.1)
                hold off
                ylabel('Firing rate (Hz)')
                % xlabel('Time from takeoffs (sec)')
                xlim(t_margin)
                % title_str = sprintf('Others'' takeoffs: p=%.3f \nGLM:%s',min(stat_test_poisson(3,:)),coef_sig{:});
                % title_str = sprintf('Others'' landings: \nGLM:%s',coef_sig{:});
                % title(title_str,'Interpreter','None')
                title(title_str,'Interpreter','None','FontSize',14)
                set(gca,'FontSize',14,'XColor','none','TickDir','none')
            
                
                % Plot 2: raster plot
                % nexttile
                subplot('Position',[0.1 0.18 0.8 0.32])
                t_bins = t_margin(1):1/w_step_2:t_margin(2);
                % A = [spmat2_social; spmat2_asocial]; % array for raster plot
                % A = sp_mat_2(:,:,unit);
                A = [sp_mat_2(is_social,:,unit); sp_mat_2(~is_social,:,unit)];
                A(A==0) = nan; A(A>0) = 1; % convert to binary
                h = imagesc([t_margin(1) t_margin(2)],[0.5 n1+n2-0.5],A);
                set(gca,'YDir','normal','FontSize',16)
                set(h,'AlphaData',~isnan(A))
                colormap gray
            
                hold on
                yregion(0,n1,'FaceColor',[250 50 50]/255,'FaceAlpha',0.2)
                yregion(n1,n1+n2,'FaceColor',[100 100 100]/255,'FaceAlpha',0.2)
                hold off
                xlim(t_margin)
                ylim([0,n1+n2])
                xlabel(sprintf('Time from others'' %s (sec)',target))
                yticks([0.5, n1-0.5, n1+n2-0.5]);
                yticklabels({1,n1,n1+n2})
                set(gca,'TickDir','none','FontSize',14)
            
                % the arrows
                xO = -0.2;  
                yO = 0.5;
                x_fd1 = t_margin(1);
                x_fd2 = x_fd1;
                y_fd1 = round(n1/2);
                y_fd2 = round(n1+n2/2);
                % patch(...
                %     [x_fd1+xO x_fd2+xO; x_fd1+xO x_fd1+xO; x_fd1 x_fd2], ...
                %     [y_fd1-yO y_fd2-yO; y_fd1+yO y_fd2+yO; y_fd1 y_fd2], 'k', 'clipping', 'off')
            
                % Texts
                text(x_fd1-0.1, y_fd1, 'Social', 'fontsize', 14,'HorizontalAlignment','center','Rotation',90)
                text(x_fd2-0.1, y_fd2, 'Non-social', 'fontsize', 14,'HorizontalAlignment','center','Rotation',90)
            
            
                % Plot 3: flight structure
                subplot('Position',[0.1 0.05 0.8 0.05])
                hold on
                xregion([find(is_social) find(is_social)+1],'FaceColor',[250 50 50]/255,'FaceAlpha',0.2)
                xregion([find(~is_social) find(~is_social)+1],'FaceColor',[100 100 100]/255,'FaceAlpha',0.2)
                hold off
                xlim([1,n1+n2])
                xticks([1,n1+n2])
                set(gca,'FontSize',14,'TickDir','none','YColor','none')
            
                %%% Save figure
                saveDir = fullfile(FigDir,bname,sprintf('Clu%d',clu_id),region);
                if isempty(dir(saveDir)); mkdir(saveDir); end
                saveas(fig,fullfile(saveDir,sprintf('Unit%d.jpg',unit_id)))
                close(fig)
               
            end
        end
    end
end
fr_field = {'bname','unit_id','fr','clu_id','p_ctest','p_anova','glm_social','fr1','fr2'};
save(fullfile(FigDir,'pefr_all.mat'),'fr_all','fr_field')