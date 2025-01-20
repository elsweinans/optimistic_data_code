%% OPTIMISTIC BASELINE ANALYSES | DM1 VS HEALTHY CONTROLS
% Written by Jerrald Rector & Els Weinans 13-Jan-2022

% Explore full dataset and run binary logistic regression models 

clear; close all; clc;

% Load full dataset w/age, sex, Cd, and fragmentation metrics 
load('dataBL.mat') % Latest version with new corrdims

%% Test correlation between variables of interest

vnames = {'avgs','stds','CV','ACs','corrdims','ZkAR','ZkRA'};
y = table2array(dataBL(:,vnames));
[testcoeff, testp] = corr(y);

A = tril(testcoeff,-1); B = tril(testp,-1); % Extract lower triangle only
% Concatenate coefficients and p-values (rounded to 3 places)
concattest = arrayfun(@(x,y) [num2str(round(x,3)) ' (' num2str(round(y,3)) ')'],A,B,'un',0);
emp = concattest{1,1}; % find template empty cell
concattest(contains(concattest, emp)) = {'-'};
% Collect in prelim table
corcoeftbl = array2table(concattest,'RowNames',vnames,'VariableNames',vnames);

% corcoeftbl =
% 
%   6�6 table
% 
%                      avgs               stds                CV               corrdims              ZkAR          ZkRA
%                 _______________    _______________    _______________    ________________    ________________    ____
% 
%     avgs        '-'                '-'                '-'                '-'                 '-'                 '-' 
%     stds        '0.938 (0)'        '-'                '-'                '-'                 '-'                 '-' 
%     CV          '0.211 (0.063)'    '0.502 (0)'        '-'                '-'                 '-'                 '-' 
%     corrdims    '-0.544 (0)'       '-0.652 (0)'       '-0.468 (0)'       '-'                 '-'                 '-' 
%     ZkAR        '-0.584 (0)'       '-0.629 (0)'       '-0.39 (0)'        '0.489 (0)'         '-'                 '-' 
%     ZkRA        '0.345 (0.002)'    '0.303 (0.007)'    '0.049 (0.671)'    '-0.342 (0.002)'    '-0.306 (0.006)'    '-' 

% Note: Correlations between predictors are all highly significant
%       Range |.21 to .94| - May pose multicolinearity issues...

% Save table to spreadsheet
% writetable(corcoeftbl,'OptimPAprelimStats_corcoef.xlsx','WriteRowNames',true,'Sheet',1) % Updated version

%% Boxplots of indicators by group

figure('position',[100,100,700,260])
subplot(1,3,1)
boxplot(dataBL.avgs,dataBL.healthy)
ylabel('Average activity')
xticklabels({'DM1','Healthy'})

subplot(1,3,2)
boxplot(dataBL.CV,dataBL.healthy)
ylabel('Coefficient of variation')
xticklabels({'DM1','Healthy'})

subplot(1,3,3)
boxplot(dataBL.ACs,dataBL.healthy)
ylabel('Activity lag-1 autocorrelation')
xticklabels({'DM1','Healthy'})

f = gcf;
exportgraphics(f,'trad_boxplots2.png','Resolution',600)

figure('position',[100,100,700,260])
subplot(1,3,1)
boxplot(dataBL.corrdims,dataBL.healthy)
ylabel('Correlation dimension')
xticklabels({'DM1','Healthy'})

subplot(1,3,2)
boxplot(dataBL.kAR,dataBL.healthy)
ylabel('Active-Rest (kAR)')
xticklabels({'DM1','Healthy'})

subplot(1,3,3)
boxplot(dataBL.kRA,dataBL.healthy)
ylabel('Rest-Active (kRA)')
xticklabels({'DM1','Healthy'})

f = gcf;
exportgraphics(f,'compl_boxplots.png','Resolution',600)

% Check distributions in histograms
figure
for n = 1:length(vnames)
    subplot(3,3,n)
        histogram(dataBL.(vnames{n}))
        title(vnames{n})
end

%% Mann whitney U tests (non-parametric)
h = logical(dataBL.healthy); % Healthy flag

[p1,~,stats1] = ranksum(dataBL.avgs(~h)    ,dataBL.avgs(h))
[p2,~,stats2] = ranksum(dataBL.CV(~h)      ,dataBL.CV(h))
[p3,~,stats3] = ranksum(dataBL.ACs(~h)     ,dataBL.ACs(h))
[p4,~,stats4] = ranksum(dataBL.corrdims(~h),dataBL.corrdims(h))
[p5,~,stats5] = ranksum(dataBL.kAR(~h)     ,dataBL.kAR(h))
[p6,~,stats6] = ranksum(dataBL.kRA(~h)     ,dataBL.kRA(h))



% All significantly different at < 0.01 level

% Collect stats
test(1,1) = stats1.ranksum;
test(2,1) = stats2.ranksum;
test(3,1) = stats3.ranksum;
test(4,1) = stats4.ranksum;
test(5,1) = stats5.ranksum;
test(6,1) = stats6.ranksum;

pval(1,1) = p1;
pval(2,1) = p2;
pval(3,1) = p3;
pval(4,1) = p4;
pval(5,1) = p5;
pval(6,1) = p6;

% Gather median and IQR for each
vnames2 = {'avgs','CV','corrdims','ACs','ZkAR','ZkRA'};
for v = 1:length(vnames2)
    mDM1{v,1} = [num2str(median(dataBL.(vnames2{v})(~h),'omitnan'),3), ...
        ' (' num2str(iqr(dataBL.(vnames2{v})(~h)),3),')'];
    mHLY{v,1} = [num2str(median(dataBL.(vnames2{v})(h),'omitnan'),3), ...
        ' (' num2str(iqr(dataBL.(vnames2{v})(h)),3),')'];
end

% Combine into table
mannU = table;
mannU.ind = vnames2';
mannU.mDM1 = mDM1;
mannU.mHLY = mHLY;
mannU.U = test;
mannU.pval = pval;

% mannU =
% 
%   5�5 table
% 
%        ind              mDM1                 mHLY            U         pval   
%     __________    _________________    _________________    ____    __________
% 
%     'avgs'        '0.0172 (0.0111)'    '0.0368 (0.0164)'    1910     1.013e-06
%     'CV'          '2.66 (0.555)'       '3.15 (0.99)'        2037    0.00064789
%     'corrdims'    '1.32 (0.607)'       '0.631 (0.188)'      2857     9.182e-10
%     'ZkAR'        '0.162 (1.04)'       '-1.24 (1.54)'       2688    3.2424e-05
%     'ZkRA'        '-0.273 (1.12)'      '0.5 (1.12)'         2101     0.0076823

% NOTE: Values are median(IQR)

% Save table to spreadsheet
%loc = 'H:\Onderzoek\Jerrald\Optimistic\Optimistic Full data\';
writetable(mannU,'OptimPAprelimStats_mannU.xlsx','Sheet',3) % Updated version

%% Compare transformed kAR & kRA with t-test (parametric)
% Label outliers
ol_kAR = logical(dataBL.ZkAR > 2 | dataBL.ZkAR < -2);
ol_kRA = logical(dataBL.ZkRA > 2 | dataBL.ZkRA < -2);

olkAR = dataBL.ZkAR; % First remove outliers for kAR/kRA
olkAR(ol_kAR) = nan;
dataBL.olkAR = olkAR;

olkRA = dataBL.ZkRA;
olkRA(ol_kRA) = nan;
dataBL.olkRA = olkRA;

% [h3a,p3a] = ttest2(dataBL.olkAR(dataBL.healthy==0)  ,dataBL.olkAR(dataBL.healthy==1))
% % p = 3.1303e-06
% [h4a,p4a] = ttest2(dataBL.olkRA(dataBL.healthy==0)  ,dataBL.olkRA(dataBL.healthy==1))
% % p = 0.0093

%% Basic (unadjusted) logistic regression plots

modelspec = 'healthy ~ avgs';
mdl1 = fitglm(dataBL,modelspec,'link','logit','Distribution','binomial')
modelspec = 'healthy ~ CV';
mdl2 = fitglm(dataBL,modelspec,'link','logit','Distribution','binomial')
modelspec = 'healthy ~ ACs';
mdl3 = fitglm(dataBL,modelspec,'link','logit','Distribution','binomial')

modelspec = 'healthy ~ corrdims';
mdl4 = fitglm(dataBL,modelspec,'link','logit','Distribution','binomial')
modelspec = 'healthy ~ ZkAR';
mdl5 = fitglm(dataBL,modelspec,'link','logit','Distribution','binomial','Exclude',ol_kAR)
modelspec = 'healthy ~ ZkRA';
mdl6 = fitglm(dataBL,modelspec,'link','logit','Distribution','binomial','Exclude',ol_kRA)

figure
subplot(3,1,1)
hold on
scatter(dataBL.avgs,dataBL.healthy,'k','filled')
x=linspace(min(dataBL.avgs),max(dataBL.avgs),100);
beta = mdl1.Coefficients.Estimate;
plot(x,1./(1+exp(-(beta(1)+beta(2)*x))),'-r','linewidth',1.5)
xlabel('Average activity')
yticks([0 1])
yticklabels({'DM1','Healthy'})

subplot(3,1,2)
hold on
scatter(dataBL.CV,dataBL.healthy,'k','filled')
x=linspace(min(dataBL.CV),max(dataBL.CV),100);
beta = mdl2.Coefficients.Estimate;
plot(x,1./(1+exp(-(beta(1)+beta(2)*x))),'-r','linewidth',1.5)
xlabel('CoV Activity')
yticks([0 1])
yticklabels({'DM1','Healthy'})

subplot(3,1,3)
hold on
scatter(dataBL.ACs,dataBL.healthy,'k','filled')
x=linspace(min(dataBL.ACs),max(dataBL.ACs),100);
beta = mdl3.Coefficients.Estimate;
plot(x,1./(1+exp(-(beta(1)+beta(2)*x))),'-r','linewidth',1.5)
xlabel('AC Activity')
yticks([0 1])
yticklabels({'DM1','Healthy'})

f = gcf;
exportgraphics(f,'logmodel_trad.png','Resolution',600)

figure
subplot(3,1,1)
hold on
scatter(dataBL.corrdims,dataBL.healthy,'k','filled')
x=linspace(min(dataBL.corrdims),max(dataBL.corrdims),100);
beta = mdl4.Coefficients.Estimate;
plot(x,1./(1+exp(-(beta(1)+beta(2)*x))),'-r','linewidth',1.5)
xlabel('Correlation dimension')
yticks([0 1])
yticklabels({'DM1','Healthy'})

subplot(3,1,2)
hold on
scatter(dataBL.ZkAR,dataBL.healthy,'k','filled')
x=linspace(min(dataBL.ZkAR),max(dataBL.ZkAR),100);
beta = mdl5.Coefficients.Estimate;
plot(x,1./(1+exp(-(beta(1)+beta(2)*x))),'-r','linewidth',1.5)
xlabel('k_{AR}')
yticks([0 1])
yticklabels({'DM1','Healthy'})

subplot(3,1,3)
hold on
scatter(dataBL.ZkRA,dataBL.healthy,'k','filled')
x=linspace(min(dataBL.ZkRA),max(dataBL.ZkRA),100);
beta = mdl6.Coefficients.Estimate;
plot(x,1./(1+exp(-(beta(1)+beta(2)*x))),'-r','linewidth',1.5)
xlabel('k_{RA}')
yticks([0 1])
yticklabels({'DM1','Healthy'})

f = gcf;
exportgraphics(f,'logmodel_compl.png','Resolution',600)

%% Full (adjusted) binary logistic analysis
% Adjust progressively for age & sex, then avgs & CV (not stds)

% Define p-value stars
p05 = {'*'};
p01 = {'**'};
p001 = {'***'};
ns ={''};

%% Center variables for analysis

% dataBL.corrdimz = zscore(dataBL.corrdim);
dataBL.corrdimz = zscore(dataBL.corrdims); 
dataBL.agec = dataBL.age - mean(dataBL.age,'omitnan');
dataBL.sexc = categorical(double(dataBL.sex) - 1,[0 1],{'Female' 'Male'});
dataBL.avgsz = zscore(dataBL.avgs);
dataBL.CVc = dataBL.CV - mean(dataBL.CV,'omitnan');
dataBL.ACz = zscore(dataBL.ACs);

%% Model 1 - Each indicator individually adjusted for age & sex
modelspec = 'healthy ~ agec + sexc + corrdimz';
mdl1Cd = fitglm(dataBL,modelspec,'link','logit','Distribution','binomial')

% Estimated Coefficients:
%                    Estimate        SE        tStat        pValue  
%                    _________    ________    ________    __________
% 
%     (Intercept)      -3.8454      1.2769     -3.0116     0.0025991
%     corrdimz         -4.7794      1.3763     -3.4727    0.00051534
%     agec           -0.044391    0.038907      -1.141       0.25389
%     sexc_Male       -0.49351     0.96059    -0.51376       0.60742

% Convert results to odds ratio (95% CI) p-value
est = mdl1Cd.Coefficients.Estimate;
se  = mdl1Cd.Coefficients.SE;
OR  = exp(est);
UB  = exp(est + 1.96*se);
LB  = exp(est - 1.96*se);
pVal = mdl1Cd.Coefficients.pValue;
p = num2cell(pVal);
p(pVal <= 0.05) = p05;
p(pVal <= 0.01) = p01;
p(pVal <= 0.001) = p001;
p(pVal > 0.05) = ns;

% Get row names
coef = mdl1Cd.Coefficients;
rNames = coef.Properties.RowNames;
% allNames = rNames; % collect unique variable names for use in table

% Concatenate OR (LB - UB) & p-value
for k = 1:height(coef)
    model1cd{k,1} = rNames{k,1};
    model1cd{k,2} = [num2str(OR(k),2),' (',num2str(LB(k),2),' - ',num2str(UB(k),3),')',p{k}];
end

clear est se OR UB LB pval coef rNames


modelspec = 'healthy ~ agec + sexc + ZkAR';
mdl1kARz = fitglm(dataBL,modelspec,'link','logit','Distribution','binomial','Exclude',ol_kAR)

% Estimated Coefficients:
%                    Estimate        SE        tStat        pValue  
%                    _________    ________    ________    __________
% 
%     (Intercept)      -1.7027     0.55763     -3.0534     0.0022624
%     ZkAR             -1.6559     0.46562     -3.5563    0.00037615
%     agec           -0.090124    0.035617     -2.5304      0.011394
%     sexc_Male       -0.49489     0.70054    -0.70644       0.47991

% Convert results to odds ratio (95% CI) p-value
est = mdl1kARz.Coefficients.Estimate;
se  = mdl1kARz.Coefficients.SE;
OR  = exp(est);
UB  = exp(est + 1.96*se);
LB  = exp(est - 1.96*se);
pVal = mdl1kARz.Coefficients.pValue;
p = num2cell(pVal);
p(pVal <= 0.05) = p05;
p(pVal <= 0.01) = p01;
p(pVal <= 0.001) = p001;
p(pVal > 0.05) = ns;

% Get row names
coef = mdl1kARz.Coefficients;
rNames = coef.Properties.RowNames;
% allNames = rNames; % collect unique variable names for use in table

% Concatenate OR (LB - UB) & p-value
for k = 1:height(coef)
    model1ar{k,1} = rNames{k,1};
    model1ar{k,2} = [num2str(OR(k),2),' (',num2str(LB(k),2),' - ',num2str(UB(k),3),')',p{k}];
end

% Fix variable name for creating table
model1ar{strcmpi(model1ar,'ZkAR')} = 'olkAR';

clear est se OR UB LB pval coef rNames


modelspec = 'healthy ~ agec + sexc + ZkRA';
mdl1kRAz = fitglm(dataBL,modelspec,'link','logit','Distribution','binomial','Exclude',ol_kRA)

% Estimated Coefficients:
%                    Estimate        SE        tStat       pValue  
%                    _________    ________    ________    _________
% 
%     (Intercept)      -1.2737     0.46783     -2.7225    0.0064782
%     ZkRA             0.78369     0.37698      2.0788     0.037631
%     agec           -0.081757    0.029593     -2.7627    0.0057328
%     sexc_Male       -0.45133      0.6297    -0.71674      0.47354

% Convert results to odds ratio (95% CI) p-value
est = mdl1kRAz.Coefficients.Estimate;
se  = mdl1kRAz.Coefficients.SE;
OR  = exp(est);
UB  = exp(est + 1.96*se);
LB  = exp(est - 1.96*se);
pVal = mdl1kRAz.Coefficients.pValue;
p = num2cell(pVal);
p(pVal <= 0.05) = p05;
p(pVal <= 0.01) = p01;
p(pVal <= 0.001) = p001;
p(pVal > 0.05) = ns;

% Get row names
coef = mdl1kRAz.Coefficients;
rNames = coef.Properties.RowNames;
% allNames = rNames; % collect unique variable names for use in table

% Concatenate OR (LB - UB) & p-value
for k = 1:height(coef)
    model1ra{k,1} = rNames{k,1};
    model1ra{k,2} = [num2str(OR(k),2),' (',num2str(LB(k),2),' - ',num2str(UB(k),3),')',p{k}];
end

% Fix variable name for creating table
model1ra{strcmpi(model1ra,'ZkRA')} = 'olkRA';

clear est se OR UB LB pval coef rNames

%% Model 2 - Each indicator individually additionally adjusted for avgs + CV +AC
modelspec = 'healthy ~ agec + sexc + avgsz + CVc + ACz + corrdimz';
mdl2Cd = fitglm(dataBL,modelspec,'link','logit','Distribution','binomial')

% Estimated Coefficients:
%                    Estimate        SE        tStat       pValue  
%                    _________    ________    ________    _________
% 
%     (Intercept)      -4.4633      1.6331      -2.733    0.0062759
%     corrdimz          -4.244      1.4292     -2.9696    0.0029823
%     agec           -0.030488    0.042469    -0.71789      0.47283
%     avgsc             82.827       56.22      1.4733      0.14068
%     CVc               1.2438      1.1518      1.0799      0.28021
%     sexc_Male       -0.56131     0.99896    -0.56189      0.57419

% Convert results to odds ratio (95% CI) p-value
est = mdl2Cd.Coefficients.Estimate;
se  = mdl2Cd.Coefficients.SE;
OR  = exp(est);
UB  = exp(est + 1.96*se);
LB  = exp(est - 1.96*se);
pVal = mdl2Cd.Coefficients.pValue;
p = num2cell(pVal);
p(pVal <= 0.05) = p05;
p(pVal <= 0.01) = p01;
p(pVal <= 0.001) = p001;
p(pVal > 0.05) = ns;

% Get row names
coef = mdl2Cd.Coefficients;
rNames = coef.Properties.RowNames;
% allNames = rNames; % collect unique variable names for use in table

% Concatenate OR (LB - UB) & p-value
for k = 1:height(coef)
    model2cd{k,1} = rNames{k,1};
    model2cd{k,2} = [num2str(OR(k),2),' (',num2str(LB(k),2),' - ',num2str(UB(k),3),')',p{k}];
end

clear est se OR UB LB pval coef rNames

% test2 = table(model2cd);
%                       model2cd                  
%     ____________________________________________
% 
%     '(Intercept)'    '0.012 (0.00047 - 0.283)**'
%     'corrdimz'       '0.014 (0.00087 - 0.236)**'
%     'agec'           '0.97 (0.89 - 1.05)'       
%     'CVc'            '3.5 (0.36 - 33.2)'        
%     'sexc_Male'      '0.57 (0.081 - 4.04)'      
%     'avgsz'          '3 (0.69 - 13.2)'  


modelspec = 'healthy ~ agec + sexc + avgsz + CVc + ACz + ZkAR';
mdl2kARz = fitglm(dataBL,modelspec,'link','logit','Distribution','binomial','Exclude',ol_kAR)

% Estimated Coefficients:
%                    Estimate        SE        tStat      pValue  
%                    _________    ________    _______    _________
% 
%     (Intercept)      -2.0617     0.68334    -3.0171    0.0025518
%     ZkAR             -1.1119     0.55967    -1.9868     0.046948
%     agec           -0.050752    0.036918    -1.3747      0.16922
%     CVc               2.0628     0.90229     2.2862     0.022245
%     sexc_Male        -0.8552     0.82991    -1.0305      0.30279
%     avgsz              1.563     0.56232     2.7796    0.0054429

% Convert results to odds ratio (95% CI) p-value
est = mdl2kARz.Coefficients.Estimate;
se  = mdl2kARz.Coefficients.SE;
OR  = exp(est);
UB  = exp(est + 1.96*se);
LB  = exp(est - 1.96*se);
pVal = mdl2kARz.Coefficients.pValue;
p = num2cell(pVal);
p(pVal <= 0.05) = p05;
p(pVal <= 0.01) = p01;
p(pVal <= 0.001) = p001;
p(pVal > 0.05) = ns;

% Get row names
coef = mdl2kARz.Coefficients;
rNames = coef.Properties.RowNames;
% allNames = rNames; % collect unique variable names for use in table

% Concatenate OR (LB - UB) & p-value
for k = 1:height(coef)
    model2ar{k,1} = rNames{k,1};
    model2ar{k,2} = [num2str(OR(k),2),' (',num2str(LB(k),2),' - ',num2str(UB(k),3),')',p{k}];
end

% Fix variable name for creating table
model2ar{strcmpi(model2ar,'ZkAR')} = 'olkAR';

clear est se OR UB LB pval coef rNames


modelspec = 'healthy ~ agec + sexc + avgsz + CVc + ACz + ZkRA';
mdl2kRAz = fitglm(dataBL,modelspec,'link','logit','Distribution','binomial','Exclude',ol_kRA)

% Estimated Coefficients:
%                    Estimate        SE        tStat      pValue  
%                    _________    ________    _______    _________
% 
%     (Intercept)      -9.3882      3.0125    -3.1164    0.0018305
%     CV                2.5492      0.8962     2.8444    0.0044496
%     ZkRA             0.97675     0.54826     1.7816     0.074821
%     agec           -0.037106    0.034712     -1.069      0.28508
%     sexc_Male       -0.93356     0.83839    -1.1135      0.26549
%     avgsz             1.6058     0.51848     3.0972    0.0019535

% Convert results to odds ratio (95% CI) p-value
est = mdl2kRAz.Coefficients.Estimate;
se  = mdl2kRAz.Coefficients.SE;
OR  = exp(est);
UB  = exp(est + 1.96*se);
LB  = exp(est - 1.96*se);
pVal = mdl2kRAz.Coefficients.pValue;
p = num2cell(pVal);
p(pVal <= 0.05) = p05;
p(pVal <= 0.01) = p01;
p(pVal <= 0.001) = p001;
p(pVal > 0.05) = ns;

% Get row names
coef = mdl2kRAz.Coefficients;
rNames = coef.Properties.RowNames;
% allNames = rNames; % collect unique variable names for use in table

% Concatenate OR (LB - UB) & p-value
for k = 1:height(coef)
    model2ra{k,1} = rNames{k,1};
    model2ra{k,2} = [num2str(OR(k),2),' (',num2str(LB(k),2),' - ',num2str(UB(k),3),')',p{k}];
end

% Fix variable name for creating table
model2ra{strcmpi(model2ra,'ZkRA')} = 'olkRA';

clear est se OR UB LB pval coef rNames



% Model 2a - Both fragmentation indicators entered in one model[?]
    % With outliers removed (see code above)
modelspec = 'healthy ~ agec + sexc + avgsz + CVc + ACz + olkAR + olkRA';
mdl2kXXz = fitglm(dataBL,modelspec,'link','logit','Distribution','binomial')

% Estimated Coefficients:
%                    Estimate        SE        tStat      pValue  
%                    _________    ________    _______    _________
% 
%     (Intercept)      -2.5142     0.86149    -2.9185    0.0035174
%     olkAR            -1.2072     0.58676    -2.0575     0.039641
%     olkRA            0.97749      0.5742     1.7024     0.088688
%     agec           -0.049101    0.037199      -1.32      0.18685
%     CVc               2.7111      1.1007     2.4631     0.013774
%     sexc_Male        -1.1246       0.925    -1.2158      0.22405
%     avgsz             1.4384     0.58467     2.4602     0.013885

% Convert results to odds ratio (95% CI) p-value
est = mdl2kXXz.Coefficients.Estimate;
se  = mdl2kXXz.Coefficients.SE;
OR  = exp(est);
UB  = exp(est + 1.96*se);
LB  = exp(est - 1.96*se);
pVal = mdl2kXXz.Coefficients.pValue;
p = num2cell(pVal);
p(pVal <= 0.05) = p05;
p(pVal <= 0.01) = p01;
p(pVal <= 0.001) = p001;
p(pVal > 0.05) = ns;

% Get row names
coef = mdl2kXXz.Coefficients;
rNames = coef.Properties.RowNames;
% allNames = rNames; % collect unique variable names for use in table

% Concatenate OR (LB - UB) & p-value
for k = 1:height(coef)
    model2kxx{k,1} = rNames{k,1};
    model2kxx{k,2} = [num2str(OR(k),2),' (',num2str(LB(k),2),' - ',num2str(UB(k),3),')',p{k}];
end

clear est se OR UB LB pval coef rNames

% Decide whether to use this model 2a (or the 2 above with individual kAR/kRA values)


% Model 3 - All metrics in one model together (maybe issues)
modelspec = 'healthy ~ agec + sexc + avgsz + CVc + olkAR + olkRA + corrdimz';
mdl3kXXCd = fitglm(dataBL,modelspec,'link','logit','Distribution','binomial')

% Estimated Coefficients:
%                    Estimate      SE        tStat       pValue 
%                    ________    _______    ________    ________
% 
%     (Intercept)     -10.298     5.3234     -1.9344    0.053065
%     olkAR           -5.2687     3.1412     -1.6773    0.093486
%     olkRA            2.3719     1.6285      1.4565     0.14526
%     corrdimz        -10.586     5.2529     -2.0152    0.043882
%     agec           -0.18639    0.13446     -1.3862     0.16569
%     CVc            0.097675     3.0122    0.032427     0.97413
%     sexc_Male       -2.6597      2.729     -0.9746     0.32976
%     avgsz          -0.84296     1.6399    -0.51404     0.60722

% Convert results to odds ratio (95% CI)
est = mdl3kXXCd.Coefficients.Estimate;
se  = mdl3kXXCd.Coefficients.SE;
OR  = exp(est);
UB  = exp(est + 1.96*se);
LB  = exp(est - 1.96*se);
pVal = mdl3kXXCd.Coefficients.pValue;
p = num2cell(pVal);
p(pVal <= 0.05) = p05;
p(pVal <= 0.01) = p01;
p(pVal <= 0.001) = p001;
p(pVal > 0.05) = ns;

% Get row names
coef = mdl3kXXCd.Coefficients;
rNames = coef.Properties.RowNames;
allNames = rNames; % collect unique variable names for use in table

% Concatenate OR (LB - UB) & p-value
for k = 1:height(coef)
    model3all{k,1} = rNames{k,1};
    model3all{k,2} = [num2str(OR(k),2),' (',num2str(LB(k),2),' - ',num2str(UB(k),3),')',p{k}];
end

clear est se OR UB LB pval coef %rNames

%% Combine results into table

Tbl1 = table; % store results
Tbl1.var = rNames;
for v = 1:length(rNames) % Loop thru variables (of largest nested model)
    for m1 = 1:length(model1cd)
       if strcmpi(rNames(v,1),model1cd(m1,1))
           Tbl1.Model_1(v) = model1cd(m1,2);
       end
    end
    for m2 = 1:length(model1ar)
       if strcmpi(rNames(v,1),model1ar(m2,1))
           Tbl1.Model_2(v) = model1ar(m2,2);
       end
    end
    for m3 = 1:length(model1ra)
       if strcmpi(rNames(v,1),model1ra(m3,1))
           Tbl1.Model_3(v) = model1ra(m3,2);
       end
    end
    for m4 = 1:length(model2cd)
       if strcmpi(rNames(v,1),model2cd(m4,1))
           Tbl1.Model_4(v) = model2cd(m4,2);
       end
    end
    for m5 = 1:length(model2ar)
       if strcmpi(rNames(v,1),model2ar(m5,1))
           Tbl1.Model_5(v) = model2ar(m5,2);
       end
    end
    for m6 = 1:length(model2ra)
       if strcmpi(rNames(v,1),model2ra(m6,1))
           Tbl1.Model_6(v) = model2ra(m6,2);
       end
    end
    for m7 = 1:length(model3all)
       if strcmpi(rNames(v,1),model3all(m7,1))
           Tbl1.Model_7(v) = model3all(m7,2);
       end
    end
end

% Rearrange rows 
Tbl1_1 = Tbl1([5 6 7 8 4 2 3],:);

%% Save table in Excel sheet
writetable(Tbl1_1,'OptimPAprelimStats2.xlsx','Sheet',2) % Updated version

