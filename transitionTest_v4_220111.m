%% Calculate active-rest transition probabilities based on Lim et al 2011
% Written by Jerrald Rector on 11-Nov-2021

%  Quantify fragmentation of activity for comparison across groups
%  using a more intuitive interpretation than the alpha/GINI

clear; close all; clc;

cd 'H:\Onderzoek\Jerrald\Optimistic\Optimistic Full data'

% Load data created by "OptimPA_datasetcreation_210212.m"
% load('OptimData4.mat')

% Note: Data are Euclidian norm minus one (ENMO) sampled 1x/5 seconds

%% Begin transition analysis 
% tic
% for ID = 1:size(OptimPA,2) 
%     IDstr = num2str(ID);
%     ind = OptimPA(ID).grpind;
%     timept = OptimPA(ID).timept;
% 
%     %% Collect active/rest bout data & remove padding
% 
%     actfreq = OptimPA(ID).actfreq;
%     sedfreq = OptimPA(ID).sedfreq;
% 
%     % Remove zero padding 
%     a_pad = actfreq(:,2)==0; % Find zeros
%     actfreq(a_pad,:) = []; % Drop unobserved bout lengths
%     clear a_pad
% 
%     r_pad = sedfreq(:,2)==0; % Find zeros
%     sedfreq(r_pad,:) = [];
%     clear r_pad
% 
% 
%     % Visualize in histogram
% %     figure
% %     subplot(211)
% %         bar(actfreq(:,1),actfreq(:,2))
% %         ylabel('# Runs of Activity')
% %     subplot(212)
% %         bar(sedfreq(:,1),sedfreq(:,2))
% %         ylabel('# Runs of Rest')
% %         xlabel('Duration of run in minutes')
% 
%     %% Collect percent of runs of duration >= t
% 
%     d_a = diff(actfreq(:,1)); % Define duration step intervals
%     d_r = diff(sedfreq(:,1));
% 
%     for t_a = length(actfreq):-1:1
%         NtA(t_a,1) = sum(actfreq(t_a:end,2))/sum(actfreq(:,2))*100; 
%     end
% 
%     for t_r = length(sedfreq):-1:1
%         NtR(t_r,1) = sum(sedfreq(t_r:end,2))/sum(sedfreq(:,2))*100; 
%     end
% 
%     % plot Nt
%     % figure
%     % subplot(211)
%     %     bar(actfreq(:,1),NtA)
%     %     ylabel('% Runs of Activity of duration >= t')
%     % subplot(212)
%     %     bar(sedfreq(:,1),NtR)
%     %     ylabel('% Runs of Rest of duration >= t')
%     %     xlabel('t = time in minutes')
% 
% 
%     % Calculate pAR(t)
%     pAR = (NtA(1:end-1)-NtA(2:end))./(NtA(1:end-1).*d_a);
%     pRA = (NtR(1:end-1)-NtR(2:end))./(NtR(1:end-1).*d_r);
% 
%     % Plot pAR(t) & pRA(t)
% %     figure
% %     subplot(211)
% %         plot(actfreq(1:end-1,1),pAR,'o')
% %         ylabel('A->R transition probability')
% %         ylim([0 1])
% %     subplot(212)
% %         plot(sedfreq(1:end-1,1),pRA,'o')
% %         ylabel('R->A transition probability')
% %         xlabel('t = time in minutes')
% %         ylim([0 1])
% 
%     %% Fit LOWESS regression & identify constant region
%     tA = actfreq(1:end-1,1);
%     wA = 0.3*length(tA); % Set span at 30% of total # of data pts
%     [trendA,winA] = smoothdata(pAR,'lowess',wA,'SamplePoints', tA);
% 
%     tR = sedfreq(1:end-1,1);
%     wR = 0.3*length(tR); 
%     [trendR,winR] = smoothdata(pRA,'lowess',wR,'SamplePoints', tR);
% 
%     % figure
%     % subplot(211)
%     %     plot(actfreq(1:end-1,1),pAR,'o')
%     %     hold on
%     %     plot(tA,trendA)
%     %     hold off
%     %     ylabel('')
%     %     ylim([0 1])
%     % subplot(212)
%     %     plot(sedfreq(1:end-1,1),pRA,'o')
%     %     hold on
%     %     plot(tR,trendR)
%     %     hold off
%     %     ylim([0 1])
% 
%     % Find constant region
%     sdTrendA = std(trendA,'omitnan'); % standard deviation of pAR(t)
%     upperA = mean(trendA) + sdTrendA;
%     lowerA = mean(trendA) - sdTrendA;
%     inboundA = trendA >= lowerA & trendA <= upperA;
% 
%     sdTrendR = std(trendR,'omitnan'); % standard deviation of pRA(t)
%     upperR = mean(trendR) + sdTrendR;
%     lowerR = mean(trendR) - sdTrendR;
%     inboundR = trendR >= lowerR & trendR <= upperR;
% 
%     % Find longest stretch of inbound (code from 
%     % https://www.mathworks.com/matlabcentral/answers/
%     % 404502-find-the-longest-sequence-of-consecutive-non-zero-values-in-a-vector)
%     zposA = find(~[0; inboundA; 0]); % Find zero positions
%     [~, grpidxA] = max(diff(zposA)); % Find largest gap between zeros = longest stretch
%     cregA = [zposA(grpidxA),zposA(grpidxA+1)-2]; % Note start and stop index for constant region
%     regionA = pAR(zposA(grpidxA):zposA(grpidxA+1)-2); % Collect pAR values within longest stretch
% 
%     zposR = find(~[0; inboundR; 0]); % Find zero positions
%     [~, grpidxR] = max(diff(zposR)); % Find largest gap between zeros = longest stretch
%     cregR = [zposR(grpidxR),zposR(grpidxR+1)-2]; % Note start and stop index for constant region
%     regionR = pRA(zposR(grpidxR):zposR(grpidxR+1)-2);
% 
%     % Calculated weighted average of constant region
%     wgtsA = sqrt(actfreq(zposA(grpidxA):zposA(grpidxA+1)-2,2));
%     wgtsA = wgtsA./max(wgtsA);
%     kAR = mean(regionA.*wgtsA);
% 
%     wgtsR = sqrt(sedfreq(zposR(grpidxR):zposR(grpidxR+1)-2,2));
%     wgtsR = wgtsR./max(wgtsR);
%     kRA = mean(regionR.*wgtsR);
% 
%     % Plot kAR against pAR / kRA against pRA
%     figure('units','normalized','outerposition',[0 0 .75 0.75]);
%     subplot(121)
%         plot(actfreq(1:end-1,1),pAR,'o');
%         hold on
%         plot(tA,trendA);
%         plot(actfreq(1:end-1,1),ones(length(pAR))*kAR,'--');
%         plot(cregA(1):cregA(2),0.5*ones(1,cregA(2)-cregA(1)+1),'k');
%         hold off
%         ylim([0 1])
%         ylabel('A -> R Transition proability')
%         legend('pAR','lowess fit', sprintf('k_{AR} = %.3f', kAR),'Constant region')
%     subplot(122)
%         plot(sedfreq(1:end-1,1),pRA,'o')
%         hold on
%         plot(tR,trendR)
%         plot(sedfreq(1:end-1,1),ones(length(pRA))*kRA,'--')
%         plot(cregR(1):cregR(2),0.5*ones(1,cregR(2)-cregR(1)+1),'k')
%         hold off
%         ylim([0 1])
%         xlim([0 150])
%         ylabel('R -> A Transition proability')
%         legend('pRA','lowess fit', sprintf('k_{RA} = %.3f', kRA),'Constant region')
%     suptitle(['Transition probabilities for ' ind,' - ' char(timept)])
% 
%     % Save plot with calculated parameters
%     loc = 'H:\Onderzoek\Jerrald\Optimistic\Optimistic Full data\OptimPA_transition_plots\';
%     print([loc, IDstr, '_OptimPA - ', ind,' - ' char(timept) '.png'],'-dpng')  
%     
%     close
% 
%     % Add relevant variables to OptimPA struct
%     OptimPA(ID).kAR  = kAR;
%     OptimPA(ID).kRA  = kRA;
%     OptimPA(ID).t_cA = cregA(1);
%     OptimPA(ID).t_cR = cregR(1);
% 
%     clearvars -except OptimPA % clean up
% 
% end
% toc

%% Save dataset so far

% save('OptimData5.mat','OptimPA');

%% Begin here with healthy/DM1 comparison of kAR & kRA

load('OptimData5.mat')

%% Create dataset for analysis

trans = struct2table(OptimPA);
trans = trans(:,[20, 18, 19, 17, 11, 16, 21:24]);

%% Plot transition vars against alpha/GINI
figure
subplot(221); scatter(trans.alpha,trans.kAR); ylabel('kAR'); xlabel('alpha')
subplot(222); scatter(trans.alpha,trans.kRA); ylabel('kRA'); xlabel('alpha')
subplot(223); scatter(trans.alpha,trans.Gindex); ylabel('kAR'); xlabel('GINI')
subplot(224); scatter(trans.alpha,trans.Gindex); ylabel('kRA'); xlabel('GINI')

%% Plot difference between groups
BL = logical(trans.PrePost==0);
FU = logical(trans.PrePost==1); % create filter to exclude follow up

kAR = trans.kAR;
kRA = trans.kRA;
trans.logitkAR = log(kAR./(1 - kAR));
trans.logitkRA = log(kRA./(1 - kRA));

figure
subplot(211); histogram(kAR); xlabel('Active -> Rest')
subplot(212); histogram(kRA); xlabel('Rest -> Active')

figure
subplot(211); histogram(logitkAR); xlabel('Active -> Rest')
subplot(212); histogram(logitkRA); xlabel('Rest -> Active')

figure
subplot(121)
    boxplot(kRA(BL==1),trans.healthy(BL==1))
%     xlabel('Group')
    ylabel('Rest -> Active probability')
    xticklabels({'DM1','Healthy'})
    axis square
subplot(122)
    boxplot(kAR(BL==1),trans.healthy(BL==1))
%     xlabel('Group')
    ylabel('Active -> Rest probability')
    xticklabels({'DM1','Healthy'})
    axis square

%% Prelim statistical comparison

modelspec = 'healthy ~ ZkAR';
mdlkARz = fitglm(trans,modelspec,'link','logit','Distribution','binomial', 'Exclude', FU)

% Estimated Coefficients:
%                    Estimate      SE        tStat       pValue  
%                    ________    _______    _______    __________
% 
%     (Intercept)     2.6433     0.98816      2.675     0.0074737
%     kAR            -49.786      13.808    -3.6055    0.00031153
% 
% 
% 78 observations, 76 error degrees of freedom
% Dispersion: 1
% Chi^2-statistic vs. constant model: 21, p-value = 4.69e-06

modelspec = 'healthy ~ kRA';
mdlkRA = fitglm(trans,modelspec,'link','logit','Distribution','binomial', 'Exclude', FU)

% Estimated Coefficients:
%                    Estimate      SE       tStat       pValue  
%                    ________    ______    _______    __________
% 
%     (Intercept)    -4.2094     1.2485    -3.3715    0.00074762
%     kRA             160.27     62.316      2.572      0.010112
% 
% 
% 78 observations, 76 error degrees of freedom
% Dispersion: 1
% Chi^2-statistic vs. constant model: 8.34, p-value = 0.00388

%% Compare pre-post among patients
healthy = trans.healthy==1;
% modelspec = 'logitkRA ~ 1 + PrePost*grp + (1 + PrePost | subject)';
% % modelspec = 'logit(kRA) ~ 1 + PrePost + (1|subject)';
% lmekRA = fitlme(trans, modelspec, 'Exclude', healthy)
% 
% 
% modelspec = 'logitkAR ~ 1 + PrePost*grp + (1|subject)';
% % modelspec = 'logit(kAR) ~ 1 + PrePost + (1|subject)';
% % modelspec = 'logitkAR ~ 1 + grp + (1|subject)';
% lmekAR = fitlme(trans, modelspec, 'Exclude', healthy)
 
interven = double(logical(trans.grp == 3 | trans.grp ==4));
interven(healthy) = nan;

respond = double(logical(trans.grp == 2 | trans.grp == 4));
respond(healthy) = nan;

trans.interven = interven;
trans.respond = respond;

trans.subject  = categorical(trans.subject);
trans.healthy  = categorical(trans.healthy);
trans.PrePost  = categorical(trans.PrePost);
trans.interven = categorical(trans.interven);
trans.respond  = categorical(trans.respond);

% Fit models
modelspec = 'logitkRA ~ 1 + PrePost*interven*respond + (1|subject)';
lmekRA = fitlme(trans, modelspec, 'Exclude', healthy)

% Model fit statistics:
%     AIC      BIC       LogLikelihood    Deviance
%     15.45    43.157    2.2748           -4.5497 
% 
% Fixed effects coefficients (95% CIs):
%     Name                              Estimate     SE          tStat       DF     pValue        Lower        Upper   
%     '(Intercept)'                       -4.1036    0.064927     -63.203    110    2.7073e-88      -4.2322     -3.9749
%     'PrePost'                         -0.021466    0.067839    -0.31643    110       0.75228     -0.15591     0.11297
%     'interven'                         0.070558     0.09182     0.76844    110       0.44387     -0.11141     0.25252
%     'respond'                          0.086313    0.093446     0.92368    110       0.35768    -0.098874      0.2715
%     'PrePost:interven'                 -0.10263    0.095938     -1.0698    110       0.28706     -0.29276    0.087494
%     'PrePost:respond'                 -0.010656    0.097636    -0.10914    110       0.91329     -0.20415     0.18284
%     'interven:respond'                -0.071901     0.13101    -0.54883    110       0.58423     -0.33153     0.18773
%     'PrePost:interven:respond'          0.17319     0.13688      1.2652    110       0.20846    -0.098082     0.44446
% 
% Random effects covariance parameters (95% CIs):
% Group: subject (59 Levels)
%     Name1                Name2                Type         Estimate    Lower      Upper  
%     '(Intercept)'        '(Intercept)'        'std'        0.16946     0.12447    0.23071
% 
% Group: Error
%     Name             Estimate    Lower      Upper  
%     'Res Std'        0.18578     0.15511    0.22252

modelspec = 'logitkAR ~ 1 + PrePost*interven*respond + (1|subject)';
lmekAR = fitlme(trans, modelspec, 'Exclude', healthy)

% Model fit statistics:
%     AIC       BIC       LogLikelihood    Deviance
%     175.93    203.63    -77.963          155.93  
% 
% Fixed effects coefficients (95% CIs):
%     Name                              Estimate    SE         tStat      DF     pValue        Lower       Upper    
%     '(Intercept)'                      -2.0331    0.12682    -16.032    110    1.4971e-30     -2.2844      -1.7817
%     'PrePost'                          0.12257    0.13718    0.89347    110       0.37355    -0.14929      0.39443
%     'interven'                        -0.44549    0.17934     -2.484    110        0.0145    -0.80091    -0.090071
%     'respond'                         -0.39575    0.18252    -2.1683    110      0.032293    -0.75746    -0.034039
%     'PrePost:interven'                 0.16763      0.194    0.86406    110       0.38943    -0.21684       0.5521
%     'PrePost:respond'                  0.18772    0.19744    0.95078    110        0.3438    -0.20355      0.57899
%     'interven:respond'                  0.7221    0.25589      2.822    110     0.0056654       0.215       1.2292
%     'PrePost:interven:respond'        -0.60246     0.2768    -2.1765    110      0.031654      -1.151    -0.053912
% 
% Random effects covariance parameters (95% CIs):
% Group: subject (59 Levels)
%     Name1                Name2                Type         Estimate    Lower      Upper  
%     '(Intercept)'        '(Intercept)'        'std'        0.31638     0.22679    0.44135
% 
% Group: Error
%     Name             Estimate    Lower      Upper  
%     'Res Std'        0.37568     0.31366    0.44997

%% Plot interaction using fixed effects estimates

PrePost  = [0; 0; 0; 0; 1; 1; 1; 1];
interven = [0; 0; 1; 1; 0; 0; 1; 1];
respond  = [0; 1; 0; 1; 0; 1; 0; 1];

betaAR = fixedEffects(lmekAR);
for ii = 1:length(PrePost)
    yhatAR(ii,1) = betaAR(1) + betaAR(2)*PrePost(ii) + betaAR(3)*interven(ii)...
         + betaAR(4)*respond(ii) + betaAR(5)*PrePost(ii)*interven(ii)...
         + betaAR(6)*PrePost(ii)*respond(ii) + betaAR(7)*interven(ii)*respond(ii)...
         + betaAR(8)*PrePost(ii)*interven(ii)*respond(ii);
end

betaRA = fixedEffects(lmekRA);
for ii = 1:length(PrePost)
    yhatRA(ii,1) = betaRA(1) + betaRA(2)*PrePost(ii) + betaRA(3)*interven(ii)...
         + betaRA(4)*respond(ii) + betaRA(5)*PrePost(ii)*interven(ii)...
         + betaRA(6)*PrePost(ii)*respond(ii) + betaRA(7)*interven(ii)*respond(ii)...
         + betaRA(8)*PrePost(ii)*interven(ii)*respond(ii);
end

% Plot
figure
subplot(121)
    plot([0; 1],[yhatAR(1),yhatAR(5)],'-o','Linewidth',2)
    hold on
    plot([0; 1],[yhatAR(2),yhatAR(6)],'-o','Linewidth',2)
    plot([0; 1],[yhatAR(3),yhatAR(7)],'-o','Linewidth',2)
    plot([0; 1],[yhatAR(4),yhatAR(8)],'-o','Linewidth',2)
    hold off
    legend('con-noresp','con-resp','int-noresp','int-resp','Location','Southeast')
    xlabel('Time point'); ylabel('logit(kAR)')
    xlim([-0.1 1.1]); ylim([-2.5 -1.8])
    xticks([0 1])
    xticklabels({'Baseline','Follow up'})
    title('Active -> Rest Probability')
subplot(122)
    plot([0; 1],[yhatRA(1),yhatRA(5)],'-o','Linewidth',2)
    hold on
    plot([0; 1],[yhatRA(2),yhatRA(6)],'-o','Linewidth',2)
    plot([0; 1],[yhatRA(3),yhatRA(7)],'-o','Linewidth',2)
    plot([0; 1],[yhatRA(4),yhatRA(8)],'-o','Linewidth',2)
    hold off
%     legend('con-noresp','con-resp','int-noresp','int-resp','Location','Southwest')
    xlabel('Time point'); ylabel('logit(kRA)')
    xlim([-0.1 1.1]); ylim([-4.2 -3.9])
    xticks([0 1])
    xticklabels({'Baseline','Follow up'})
    title('Rest -> Active Probability')
    
%% Export kAR and kRA values in csv for Els

% Copy OptimPA struct, convert to table, & remove unnecessary fields

OptimPA2 = struct2table(OptimPA);
remove = fieldnames(OptimPA);
remove([1:3,11:13,16:24]) = []; % delete all except these fields
OptimPA2 = removevars(OptimPA2, remove);
% Rearrange columns 
OptimPA2 = OptimPA2(:,[11 2 1 8 3 9 10 12 13]);

% Save dataset
save('transition.mat','OptimPA2')
% Export table using writetable
writetable(OptimPA2,'OptimisticFrag2.csv')



