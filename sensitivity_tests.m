%% Basic test and sensitivity to parameters

%% basic analysis
% load('cont_noresp_group.mat')
% load('cont_resp_group.mat')
% load('int_noresp_group.mat')
% load('int_resp_group.mat')
% load('healthy_group.mat')

timing = "b";
smoothwindow = 1;
avgwindow = 12;
ranger = [0.03 0.3];
D = 6;
distancetype = "euclidean";
n=51840;

group = "cont_noresp";
% slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
% save('z_slopes_cont_noresp4_ranger003_03','slopes')

[avgs,stds,AC,avgsraw,stdsraw] = calc_means_stds(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('cont_noresp_avg_std','avgs','stds','AC','avgsraw','stdsraw')  

group = "cont_resp";
% slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
% save('z_slopes_cont_resp4_ranger003_03','slopes')

[avgs,stds,AC,avgsraw,stdsraw] = calc_means_stds(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('cont_resp_avg_std','avgs','stds','AC','avgsraw','stdsraw') 

group = "int_noresp";
% slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
% save('z_slopes_int_noresp4_ranger003_03','slopes')

[avgs,stds,AC,avgsraw,stdsraw] = calc_means_stds(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('int_noresp_avg_std','avgs','stds','AC','avgsraw','stdsraw') 

group = "int_resp";
% slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
% save('z_slopes_int_resp4_ranger003_03','slopes')

[avgs,stds,AC,avgsraw,stdsraw] = calc_means_stds(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('int_resp_avg_std','avgs','stds','AC','avgsraw','stdsraw') 

group = "healthy";
% slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
% save('z_slopes_healthy4_ranger003_03','slopes')

[avgs,stds,AC,avgsraw,stdsraw] = calc_means_stds(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('healthy_avg_std','avgs','stds','AC','avgsraw','stdsraw') 

load('dataBL')

corrdims=[];
load('z_slopes_cont_noresp4_ranger003_03');
corrdims=[corrdims; slopes];
load('z_slopes_cont_resp4_ranger003_03');
corrdims=[corrdims; slopes];
load('z_slopes_int_noresp4_ranger003_03');
corrdims=[corrdims; slopes];
load('z_slopes_int_resp4_ranger003_03');
corrdims=[corrdims; slopes];
load('z_slopes_healthy4_ranger003_03');
corrdims=[corrdims; slopes];

avgs2=[];
load('cont_noresp_avg_std');
avgs2=[avgs2; avgs];
load('cont_resp_avg_std');
avgs2=[avgs2; avgs];
load('int_noresp_avg_std');
avgs2=[avgs2; avgs];
load('int_resp_avg_std');
avgs2=[avgs2; avgs];
load('healthy_avg_std');
avgs2=[avgs2; avgs];

dataBL.avgs=avgs2;

stds2=[];
load('cont_noresp_avg_std');
stds2=[stds2; stds];
load('cont_resp_avg_std');
stds2=[stds2; stds];
load('int_noresp_avg_std');
stds2=[stds2; stds];
load('int_resp_avg_std');
stds2=[stds2; stds];
load('healthy_avg_std');
stds2=[stds2; stds];

dataBL.stds=stds2;

ACs=[];
load('cont_noresp_avg_std');
ACs=[ACs; AC];
load('cont_resp_avg_std');
ACs=[ACs; AC];
load('int_noresp_avg_std');
ACs=[ACs; AC];
load('int_resp_avg_std');
ACs=[ACs; AC];
load('healthy_avg_std');
ACs=[ACs; AC];

dataBL.ACs=ACs;

figure
boxplot(dataBL.corrdim,dataBL.healthy)
ylabel('correlation dimension')
xticklabels({'myotonic dystrophy','healthy'})

figure
boxplot(dataBL.avgs,dataBL.healthy)
ylabel('average activity')
xticklabels({'myotonic dystrophy','healthy'})

figure
boxplot(dataBL.stds,dataBL.healthy)
ylabel('standard deviation activity')
xticklabels({'myotonic dystrophy','healthy'})

dataBL.CV=dataBL.stds./dataBL.avgs;

figure
boxplot(dataBL.CV,dataBL.healthy)
ylabel('CV activity')
xticklabels({'myotonic dystrophy','healthy'})

figure
boxplot(dataBL.ACs,dataBL.healthy)
ylabel('AC activity')
xticklabels({'myotonic dystrophy','healthy'})

%% no smooting/smoothing with window size 3, smoothing with window size 10
load('cont_noresp_group.mat')
load('cont_resp_group.mat')
load('int_noresp_group.mat')
load('int_resp_group.mat')
load('healthy_group.mat')

timing = "b";
smoothwindow = 360; %half an hour
ranger = [0.05 0.5];
D = 6;
distancetype = "euclidean";
n=51840;

group = "cont_noresp";
slopes = calc_corrdims(group,timing,smoothwindow,ranger,D,distancetype,n);
save('z_slopes_cont_noresp_ws100','slopes')

group = "cont_resp";
slopes = calc_corrdims(group,timing,smoothwindow,ranger,D,distancetype,n);
save('z_slopes_cont_resp_ws100','slopes')

group = "int_noresp";
slopes = calc_corrdims(group,timing,smoothwindow,ranger,D,distancetype,n);
save('z_slopes_int_noresp_ws100','slopes')

group = "int_resp";
slopes = calc_corrdims(group,timing,smoothwindow,ranger,D,distancetype,n);
save('z_slopes_int_resp_ws100','slopes')

group = "healthy";
slopes = calc_corrdims(group,timing,smoothwindow,ranger,D,distancetype,n);
save('z_slopes_healthy_ws100','slopes')



%% analyse sensitivity to smoothing
load('dataBL')

corrdims=[];
load('z_slopes_cont_noresp4');
corrdims=[corrdims; slopes];
load('z_slopes_cont_resp4');
corrdims=[corrdims; slopes];
load('z_slopes_int_noresp4');
corrdims=[corrdims; slopes];
load('z_slopes_int_resp4');
corrdims=[corrdims; slopes];
load('z_slopes_healthy4');
corrdims=[corrdims; slopes];

dataBL.corrdim=corrdims;

corrdimsws1=[];
load('z_slopes_cont_noresp_ws1');
corrdimsws1=[corrdimsws1; slopes];
load('z_slopes_cont_resp_ws1');
corrdimsws1=[corrdimsws1; slopes];
load('z_slopes_int_noresp_ws1');
corrdimsws1=[corrdimsws1; slopes];
load('z_slopes_int_resp_ws1');
corrdimsws1=[corrdimsws1; slopes];
load('z_slopes_healthy_ws1');
corrdimsws1=[corrdimsws1; slopes];

dataBL.corrdimws1=corrdimsws1;

% boxplots

figure
subplot(2,1,1)
boxplot(dataBL.corrdim,dataBL.healthy)
ylabel('correlation dimension')
xticklabels({'myotonic dystrophy','healthy'})

subplot(2,1,2)
boxplot(dataBL.corrdimws1,dataBL.healthy)
ylabel('correlation dimension window size = 1')
xticklabels({'myotonic dystrophy','healthy'})


%% visualisation for number of points used (example with healthy ind 002)
group = "healthy"; 
load('healthy_group.mat');
groupind=healthy_group(2);

smoothwindow = 1;
ranger = [0.05 0.5];
D = 6;
distancetype = "euclidean";
ns=[4320:4320:51840]; % 6h, 12h, 18h, 24h, 1d6h, ..., 3d
nrns=length(ns);
nrsims = 100;

filename=group + "/" + groupind + ".csv";
data = readtable(filename);
data1=data.x; 
mvavg=smoothdata(data1,'gaussian',smoothwindow);
mvavg=mvavg(5500:5500+120960); %only part with data

nrsteps=4; %10
slopes=zeros(nrns,nrsims);
for j = 8:nrns  
    j
    n=ns(j);
    for k = 1:nrsims
        startval = randi([1 120960-n+1],1,1);
        mvavg2 = mvavg(startval:startval+n);
        rs=logspace(log(ranger(1)*std(mvavg2))/log(10),log(ranger(2)*std(mvavg2))/log(10),nrsteps); %mvavg2 ipv data1
        Cds=zeros(nrsteps,1);
        for i=1:nrsteps
           r=rs(i);
           Cds(i)=calc_Cd(mvavg2,D,distancetype,r); 
        end    

        slopefit=polyfit(log(rs'),log(Cds),1);
        slopes(j,k)=slopefit(1);
        save('sensitivity_datalength_cont_noresp_1','slopes')
    end
end

dataforCI=sort(slopes'); 

figure
hold on
%errorbar(mean(slopes'),std(slopes'),'linewidth',2)
errorbar([1:12],mean(slopes'),mean(slopes')-dataforCI(5,:),dataforCI(96,:)-mean(slopes'),'linewidth',2)
%errorbar([1:12],mean(slopes'),mean(slopes')-min(slopes'),max(slopes')-mean(slopes'),'linewidth',2)
xticklabels({'0','0.5 day','1 day','1.5 days','2 days','2.5 days','3 days'})
xlabel('length of time used for analysis')
ylabel('correlation dimension (mean and CI)')


%% ranger from 0.02-0.2
load('cont_noresp_group.mat')
load('cont_resp_group.mat')
load('int_noresp_group.mat')
load('int_resp_group.mat')
load('healthy_group.mat')

timing = "b";
smoothwindow = 1; 
avgwindow = 12;
ranger = [0.02 0.2];
D = 6;
distancetype = "euclidean";
n=51840;

group = "cont_noresp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_cont_noresp_ranger02','slopes')

group = "cont_resp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_cont_resp_ranger02','slopes')

group = "int_noresp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_int_noresp_ranger02','slopes')

group = "int_resp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_int_resp_ranger02','slopes')

group = "healthy";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_healthy_ranger02','slopes')

%% ranger 0.08-0.8
clear all
close all

load('cont_noresp_group.mat')
load('cont_resp_group.mat')
load('int_noresp_group.mat')
load('int_resp_group.mat')
load('healthy_group.mat')

timing = "b";
smoothwindow = 1; %half an hour
ranger = [0.08 0.8];
D = 6;
distancetype = "euclidean";
n=51840;

group = "cont_noresp";
slopes = calc_corrdims(group,timing,smoothwindow,ranger,D,distancetype,n);
save('z_slopes_cont_noresp_ranger08','slopes')

group = "cont_resp";
slopes = calc_corrdims(group,timing,smoothwindow,ranger,D,distancetype,n);
save('z_slopes_cont_resp_ranger08','slopes')

group = "int_noresp";
slopes = calc_corrdims(group,timing,smoothwindow,ranger,D,distancetype,n);
save('z_slopes_int_noresp_ranger08','slopes')

group = "int_resp";
slopes = calc_corrdims(group,timing,smoothwindow,ranger,D,distancetype,n);
save('z_slopes_int_resp_ranger08','slopes')

group = "healthy";
slopes = calc_corrdims(group,timing,smoothwindow,ranger,D,distancetype,n);
save('z_slopes_healthy_ranger08','slopes')

%% analyse sensitivity to other ranger

load('dataBL')

% corrdims=[];
% load('z_slopes_cont_noresp2');
% corrdims=[corrdims; slopes];
% load('z_slopes_cont_resp2');
% corrdims=[corrdims; slopes];
% load('z_slopes_int_noresp2');
% corrdims=[corrdims; slopes];
% load('z_slopes_int_resp2');
% corrdims=[corrdims; slopes];
% load('z_slopes_healthy2');
% corrdims=[corrdims; slopes];
% 
% dataBL.corrdim=corrdims;

corrdimsrange02=[];
load('z_slopes_cont_noresp_ranger02');
corrdimsrange02=[corrdimsrange02; slopes];
load('z_slopes_cont_resp_ranger02');
corrdimsrange02=[corrdimsrange02; slopes];
load('z_slopes_int_noresp_ranger02');
corrdimsrange02=[corrdimsrange02; slopes];
load('z_slopes_int_resp_ranger02');
corrdimsrange02=[corrdimsrange02; slopes];
load('z_slopes_healthy_ranger02');
corrdimsrange02=[corrdimsrange02; slopes];

dataBL.corrdimrange02=corrdimsrange02;

corrdimsrange08=[];
load('z_slopes_cont_noresp_ranger08');
corrdimsrange08=[corrdimsrange08; slopes];
load('z_slopes_cont_resp_ranger08');
corrdimsrange08=[corrdimsrange08; slopes];
load('z_slopes_int_noresp_ranger08');
corrdimsrange08=[corrdimsrange08; slopes];
load('z_slopes_int_resp_ranger08');
corrdimsrange08=[corrdimsrange08; slopes];
load('z_slopes_healthy_ranger08');
corrdimsrange08=[corrdimsrange08; slopes];

dataBL.corrdimrange08=corrdimsrange08;

% boxplots

figure('position',[100,100,700,260])

subplot(1,3,1)
boxplot(dataBL.corrdimrange02,dataBL.healthy)
ylabel('correlation dimension')
xticklabels({'DM1','healthy'})
title('range r = 0.02-0.2')

subplot(1,3,2)

boxplot(dataBL.corrdims,dataBL.healthy)
ylabel('correlation dimension')
xticklabels({'DM1','healthy'})
title('range r = 0.03-0.3')

subplot(1,3,3)
boxplot(dataBL.corrdimrange08,dataBL.healthy)
ylabel('correlation dimension')
xticklabels({'DM1','healthy'})
title('range r = 0.08-0.8')

exportgraphics(gcf,'effect_r.png','Resolution',500)


%% test for full range of r (for ind int_resp A_004)
clear all
close all

% load('int_resp_group');
% groupind = int_resp_group(2);
% timing='b';
% filename="int_resp/" + groupind + "_" + timing  +".csv";

group = "healthy"; 
load('healthy_group.mat');
groupind=healthy_group(9);
filename=group + "/" + groupind + ".csv";



data = readtable(filename);
data1=data.x;
idxact=find(data1>0.05);
data2=data1(idxact(1):idxact(end));
ranger=[0.01 0.8];
D=6;
nrsteps=10;
distancetype='euclidean';
avgwindow=12;              
L=floor(length(data2)/avgwindow);
mvavg2=mean(reshape(data2(1:L*avgwindow),[avgwindow,L]),1)';


rs=logspace(log(ranger(1)*std(mvavg2))/log(10),log(ranger(2)*std(mvavg2))/log(10),nrsteps); %mvavg2 ipv data1
Cds=calc_Cd(mvavg2,D,distancetype,rs);   
figure
hold on
plot(log(rs),log(Cds),'-k','linewidth',1.5)
bound1=std(mvavg2)*0.02;
bound2=std(mvavg2)*0.06;
boundy=floor(log(min(Cds)));
v2 = [log(bound1) boundy; log(bound2) boundy; log(bound2) 0; log(bound1) 0];
f2 = [1 2 3 4];
patch('Faces',f2,'Vertices',v2,'FaceColor','red','FaceAlpha',.2);
xlabel('log(r)')
ylabel('log(Cd)')


%% D=3

load('cont_noresp_group.mat')
load('cont_resp_group.mat')
load('int_noresp_group.mat')
load('int_resp_group.mat')
load('healthy_group.mat')

timing = "b";
smoothwindow = 1; 
avgwindow = 12;
ranger = [0.03 0.3];
D = 3;
distancetype = "euclidean";
n=51840;

group = "cont_noresp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_cont_noresp_d3','slopes')

group = "cont_resp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_cont_resp_d3','slopes')

group = "int_noresp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_int_noresp_d3','slopes')

group = "int_resp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_int_resp_d3','slopes')

group = "healthy";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_healthy_d3','slopes')

%% D=9

load('cont_noresp_group.mat')
load('cont_resp_group.mat')
load('int_noresp_group.mat')
load('int_resp_group.mat')
load('healthy_group.mat')

timing = "b";
smoothwindow = 1; 
ranger = [0.03 0.3];
D = 9;
distancetype = "euclidean";
n=51840;

group = "cont_noresp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_cont_noresp_d9','slopes')

group = "cont_resp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_cont_resp_d9','slopes')

group = "int_noresp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_int_noresp_d9','slopes')

group = "int_resp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_int_resp_d9','slopes')

group = "healthy";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_healthy_d9','slopes')

%% sensitivty test d
load('dataBL')

% corrdims=[];
% load('z_slopes_cont_noresp2');
% corrdims=[corrdims; slopes];
% load('z_slopes_cont_resp2');
% corrdims=[corrdims; slopes];
% load('z_slopes_int_noresp2');
% corrdims=[corrdims; slopes];
% load('z_slopes_int_resp2');
% corrdims=[corrdims; slopes];
% load('z_slopes_healthy2');
% corrdims=[corrdims; slopes];
% 
% dataBL.corrdim=corrdims;

corrdimsd3=[];
load('z_slopes_cont_noresp_d3');
corrdimsd3=[corrdimsd3; slopes];
load('z_slopes_cont_resp_d3');
corrdimsd3=[corrdimsd3; slopes];
load('z_slopes_int_noresp_d3');
corrdimsd3=[corrdimsd3; slopes];
load('z_slopes_int_resp_d3');
corrdimsd3=[corrdimsd3; slopes];
load('z_slopes_healthy_d3');
corrdimsd3=[corrdimsd3; slopes];

dataBL.corrdimd3=corrdimsd3;

corrdimsd9=[];
load('z_slopes_cont_noresp_d9');
corrdimsd9=[corrdimsd9; slopes];
load('z_slopes_cont_resp_d9');
corrdimsd9=[corrdimsd9; slopes];
load('z_slopes_int_noresp_d9');
corrdimsd9=[corrdimsd9; slopes];
load('z_slopes_int_resp_d9');
corrdimsd9=[corrdimsd9; slopes];
load('z_slopes_healthy_d9');
corrdimsd9=[corrdimsd9; slopes];

dataBL.corrdimd9=corrdimsd9;

% boxplots

figure('position',[100,100,700,260])
subplot(1,3,1)
boxplot(dataBL.corrdimd3,dataBL.healthy)
ylabel('correlation dimension')
xticklabels({'DM1','healthy'})
title('d = 3')

subplot(1,3,2)

boxplot(dataBL.corrdims,dataBL.healthy)
ylabel('correlation dimension')
xticklabels({'DM1','healthy'})
title('d = 6')

subplot(1,3,3)
boxplot(dataBL.corrdimd9,dataBL.healthy)
ylabel('correlation dimension')
xticklabels({'DM1','healthy'})
title('d = 9')

exportgraphics(gcf,'effect_D.png','Resolution',500)



%% chebychev distance

load('cont_noresp_group.mat')
load('cont_resp_group.mat')
load('int_noresp_group.mat')
load('int_resp_group.mat')
load('healthy_group.mat')

timing = "b";
smoothwindow = 1;
avgwindow = 12;
ranger = [0.03 0.3];
D = 6;
distancetype = "chebychev";
n=51840;

group = "cont_noresp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_cont_noresp_cheb','slopes')

group = "cont_resp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_cont_resp_cheb','slopes')

group = "int_noresp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_int_noresp_cheb','slopes')

group = "int_resp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_int_resp_cheb','slopes')

group = "healthy";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_healthy_cheb','slopes')

%% analysis cheb vs eucl distance

load('dataBL')

% corrdims=[];
% load('z_slopes_cont_noresp2');
% corrdims=[corrdims; slopes];
% load('z_slopes_cont_resp2');
% corrdims=[corrdims; slopes];
% load('z_slopes_int_noresp2');
% corrdims=[corrdims; slopes];
% load('z_slopes_int_resp2');
% corrdims=[corrdims; slopes];
% load('z_slopes_healthy2');
% corrdims=[corrdims; slopes];
% 
% dataBL.corrdim=corrdims;

corrdimscheb=[];
load('z_slopes_cont_noresp_cheb');
corrdimscheb=[corrdimscheb; slopes];
load('z_slopes_cont_resp_cheb');
corrdimscheb=[corrdimscheb; slopes];
load('z_slopes_int_noresp_cheb');
corrdimscheb=[corrdimscheb; slopes];
load('z_slopes_int_resp_cheb');
corrdimscheb=[corrdimscheb; slopes];
load('z_slopes_healthy_cheb');
corrdimscheb=[corrdimscheb; slopes];

dataBL.corrdimcheb=corrdimscheb;

% boxplots

figure('position',[100,100,433,260])
subplot(1,2,1)
boxplot(dataBL.corrdims,dataBL.healthy)
ylabel('correlation dimension')
xticklabels({'DM1','healthy'})
title('euclidean distance')

subplot(1,2,2)
boxplot(dataBL.corrdimcheb,dataBL.healthy)
ylabel('correlation dimension')
xticklabels({'DM1','healthy'})
title('chebychev distance')

exportgraphics(gcf,'effect_dist.png','Resolution',500)

%% Check for 5 sec bouts
timing = "b";
smoothwindow = 1;
avgwindow=1;% 1 min bouts
ranger = [0.03 0.3];
D = 6;
distancetype = "euclidean";
n=51840;

group = "cont_noresp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_cont_noresp_s','slopes')

group = "cont_resp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_cont_resp_s','slopes')

group = "int_noresp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_int_noresp_s','slopes')

group = "int_resp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_int_resp_s','slopes')

group = "healthy";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_healthy_s','slopes')

%% Check for 1 hour bouts
timing = "b";
smoothwindow = 1;
avgwindow=12*60;% 1 hour bouts
ranger = [0.03 0.3];
D = 6;
distancetype = "euclidean";
n=51840;

group = "cont_noresp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_cont_noresp_hr','slopes')

group = "cont_resp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_cont_resp_hr','slopes')

group = "int_noresp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_int_noresp_hr','slopes')

group = "int_resp";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_int_resp_hr','slopes')

group = "healthy";
slopes = calc_corrdims(group,timing,smoothwindow,avgwindow,ranger,D,distancetype,n);
save('z_slopes_healthy_hr','slopes')

load('dataBL')

% corrdims=[];
% load('z_slopes_cont_noresp3');
% corrdims=[corrdims; slopes];
% load('z_slopes_cont_resp3');
% corrdims=[corrdims; slopes];
% load('z_slopes_int_noresp3');
% corrdims=[corrdims; slopes];
% load('z_slopes_int_resp3');
% corrdims=[corrdims; slopes];
% load('z_slopes_healthy3');
% corrdims=[corrdims; slopes];
% 
% dataBL.corrdim=corrdims;

corrdims_s=[];
load('z_slopes_cont_noresp_s');
corrdims_s=[corrdims_s; slopes];
load('z_slopes_cont_resp_s');
corrdims_s=[corrdims_s; slopes];
load('z_slopes_int_noresp_s');
corrdims_s=[corrdims_s; slopes];
load('z_slopes_int_resp_s');
corrdims_s=[corrdims_s; slopes];
load('z_slopes_healthy_s');
corrdims_s=[corrdims_s; slopes];

dataBL.corrdims_s=corrdims_s;

corrdims_hr=[];
load('z_slopes_cont_noresp_hr');
corrdims_hr=[corrdims_hr; slopes];
load('z_slopes_cont_resp_hr');
corrdims_hr=[corrdims_hr; slopes];
load('z_slopes_int_noresp_hr');
corrdims_hr=[corrdims_hr; slopes];
load('z_slopes_int_resp_hr');
corrdims_hr=[corrdims_hr; slopes];
load('z_slopes_healthy_hr');
corrdims_hr=[corrdims_hr; slopes];

dataBL.corrdims_hr=corrdims_hr;

figure('position',[100,100,700,260])
subplot(1,3,1)
boxplot(dataBL.corrdims_s,dataBL.healthy)
ylabel('correlation dimension')
xticklabels({'DM1','healthy'})
title('scale: 5 seconds')

subplot(1,3,2)
boxplot(dataBL.corrdims,dataBL.healthy)
ylabel('correlation dimension')
xticklabels({'DM1','healthy'})
title('scale: 1 min')

subplot(1,3,3)
boxplot(dataBL.corrdims_hr,dataBL.healthy)
ylabel('correlation dimension')
xticklabels({'DM1','healthy'})
title('scale: 1 hr')
exportgraphics(gcf,'effect_timescale.png','Resolution',500)

%% comparison between metrics
load('dataBL')

figure
hold on
plot(dataBL.avgs(dataBL.healthy==1),dataBL.corrdim(dataBL.healthy==1),'.k','MarkerSize',16)
plot(dataBL.avgs(dataBL.healthy==0),dataBL.corrdim(dataBL.healthy==0),'.r','MarkerSize',16)
xlabel('average activity')
ylabel('correlation dimension')
legend('healthy','DM1')

figure
hold on
plot(dataBL.stds(dataBL.healthy==1),dataBL.corrdim(dataBL.healthy==1),'.k','MarkerSize',16)
plot(dataBL.stds(dataBL.healthy==0),dataBL.corrdim(dataBL.healthy==0),'.r','MarkerSize',16)
xlabel('standard deviation activity')
ylabel('correlation dimension')
legend('healthy','DM1')

figure
hold on
plot(dataBL.CV(dataBL.healthy==1),dataBL.corrdim(dataBL.healthy==1),'.k','MarkerSize',16)
plot(dataBL.CV(dataBL.healthy==0),dataBL.corrdim(dataBL.healthy==0),'.r','MarkerSize',16)
xlabel('CV activity')
ylabel('correlation dimension')
legend('healthy','DM1')

figure
hold on
plot(dataBL.stds(dataBL.healthy==1),dataBL.avgs(dataBL.healthy==1),'.k','MarkerSize',16)
plot(dataBL.stds(dataBL.healthy==0),dataBL.avgs(dataBL.healthy==0),'.r','MarkerSize',16)
xlabel('standard deviation activity')
ylabel('average activity')
legend('healthy','DM1')

dataBL2=dataBL{:,{'avgs','stds','CV','ACs','corrdims','kRA','kAR'}};
xnames = {'avgs','stds','CV','ACs','corrdims','kRA','kAR'};  
f=figure
gplotmatrix(dataBL2,[],dataBL.healthy,[],[],[],[],'variable',xnames)
h = findobj('Tag','legend');
set(h, 'String', {'DM1','Healthy'})
indDELETE=[14 20:21 26:28 32:35 38:42 44:49]+1;
delete(f.Children(indDELETE)) 
set(f.Children(8),'Ytick',[])
set(f.Children(29),'Xtick',[])

% exportgraphics(gcf,'correlations_indicators.png','Resolution',500)

%% attractors
group = "healthy"; 
load('healthy_group.mat');
groupind=healthy_group(14);
filename=group + "/" + groupind + ".csv";
data = readtable(filename);
data1=data.x;
idxact=find(data1>0.05);
data2=data1(idxact(1):idxact(end));
data3=(data2-mean(data2))/std(data2)+mean(data2);

figure
plot(data3(1:end-1),data3(2:end),'-k','linewidth',1.5)
xlabel('data(t)')
ylabel('data(t+1)')
xlim([0 18])
ylim([0 18])

load('int_resp_group');
groupind = int_resp_group(14);
timing='b';
filename="int_resp/" + groupind + "_" + timing  +".csv";
data = readtable(filename);
data1=data.x;
idxact=find(data1>0.05);
data2=data1(idxact(1):idxact(end));
data3=(data2-mean(data2))/std(data2)+mean(data2);

figure
plot(data3(1:end-1),data3(2:end),'-k','linewidth',1.5)
xlabel('data(t)')
ylabel('data(t+1)')
xlim([0 18])
ylim([0 18])

