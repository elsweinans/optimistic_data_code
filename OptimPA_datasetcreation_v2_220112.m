%% CREATE ANALYSIS DATASET FOR OPTIMISTIC PA (INCL. AGE & SEX)
% Written by Jerrald Rector 12-Jan-2022

% Compile all necessary information into one struct/table for analyses
% 1. Table w/ baseline data for binary logistic regresssions only
% 2. Struct w/all data (i.e., pre- and post-intervention) 

clear; close all; clc;

cd 'H:\Onderzoek\Jerrald\Optimistic\Optimistic Full data' % local file path

% Begin with OptimPA dataset created by "transitionTest_v4_220111.m"
load('OptimData5.mat')

% Load minimal dataset for baseline (OptimPA2)
load('transition.mat')

cd 'H:\Onderzoek\Jerrald\Optimistic\Optimistic Full data\actometer'

% Load dataBL.mat
load('dataBL.mat')

% Import patient data
ImportPatientData_220112

% Import healthy data (Excel pw: HvS21)
ImportHealthyData_220112

% Back to main folder
cd 'H:\Onderzoek\Jerrald\Optimistic\Optimistic Full data' % local file path

%% Create examples to translate IDs into a IDmatch key

% % dataBL ID ('grpind')
dataBL.grpind = string(dataBL.grpind); % first change to string
% dataBL.grpind{1}([1,3:5])

% % healthy data
% healthy.ID{1}([1,4:6])

% % patient data
% patients.PatientID1{1}([1:4])

% % OptimPA struct
% OptimPA(50).grpind([1,3:5])

%%  Add kAR and kRA to dataBL (from OptimPA2 table)

% Remove FU data
use = OptimPA2.PrePost == 0;
OptimPA2 = OptimPA2(use,:);

for ii = 1:height(dataBL)
    matches = logical(dataBL.subject(ii) == OptimPA2.subject);
    dataBL.kRA(ii) = OptimPA2.kRA(matches);
    dataBL.kAR(ii) = OptimPA2.kAR(matches);  
end

% Standardize kAR & kRA (non-normal distribution)
kAR = dataBL.kAR;
kRA = dataBL.kRA;
logitkAR = log(kAR./(1 - kAR)); % Logit transformation (make normal)
logitkRA = log(kRA./(1 - kRA));
dataBL.ZkAR = zscore(logitkAR); % Normalize (mean = 0, SD = 1)
dataBL.ZkRA = zscore(logitkRA);


%% Add age and sex from patients and healthy tables

% Pre-allocate
dataBL.age = nan(height(dataBL),1);
dataBL.sex = nan(height(dataBL),1);

% Patients
for j = 1:height(dataBL)
    for jj = 1:height(patients)
        if logical(string(patients.PatientID1{jj}(1:4)) == string(dataBL.grpind{j}([1,3:5])))
            dataBL.age(j) = patients.Age(jj);
            dataBL.sex(j) = patients.Sex(jj);
        end
    end
end

% Healthy
for j = 60:height(dataBL) % healthy's start at 60 (save time)
    for jj = 1:height(healthy)
        if logical(string(healthy.ID{jj}([1,4:6])) == string(dataBL.grpind{j}([1,3:5])))
            dataBL.age(j) = healthy.Age(jj);
            dataBL.sex(j) = healthy.Sex(jj);
        end
    end
end

% Change sex to categorical
dataBL.sex = categorical(dataBL.sex,[1 2],{'Female' 'Male'});

%% Save dataBL

save('dataBL_analysis.mat','dataBL','patients','healthy','OptimPA2')
