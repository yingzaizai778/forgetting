close all, clear all, clc

expt = 'forget_beh';

scriptpath = mfilename('fullpath');
scriptname = mfilename();

mainPath = scriptpath(1:length(scriptpath)-length(scriptname)-length('code')-1);
dataPath = strcat(mainPath,'data/expS1/');
addpath(genpath(fullfile(mainPath, 'toolbox' , 'CircStat2012a')));
addpath(genpath(fullfile(mainPath, 'toolbox' , 'MemToolbox-master')));
%% plot error distributionss
% Load files
disp('Loading Data');
files = dir([dataPath expt '_*.csv']);
alpha = 1:360; clear guessError remError
for isub = 1:numel(files)
    data = readtable([dataPath files(isub).name], 'Delimiter',',');
    data = data{:,:};
    error = data(:,10);
    data(error==999,:) = [];
    target = data(:,7);
    error = data(:,10);

    

    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(target),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(target)))));
    end

    sampError = sort(sampError);
    guessError(isub,1) = sampError(50);

    error(error>guessError(isub,1)) = [];
    pct(isub) = numel(error)/numel(target);

    conError(isub,1) = mean(abs(error));

end


%% Experiment 1
dataPath = strcat(mainPath,'data/exp1/');
% Load files
disp('Loading Data');
files = dir([dataPath expt '_*.csv']);
alpha = 1:360;
for isub = 1:numel(files)
    data = readtable([dataPath files(isub).name], 'Delimiter',',');
    data = data{:,:};
    error = data(:,20);
    data(error==999,:) = [];

    target = data(:,17);
    error = data(:,20);

    % guess for all trials
    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(target),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(target)))));
    end

    sampError = sort(sampError);
    guessErrorAll(isub,1) = sampError(50);

    % data(error>guessErrorAll(isub),:) = [];
    target = data(:,17);
    error = data(:,20);
    whetherCued = data(:,6);
    condition = data(:,4);
    probeRem = data(:,14);
    probeNeutral = data(:,16);
    probeForget = data(:,15);
    blockIdx = data(:,3);
    blockNum = max(data(:,3));

    ineuTargRem = target((whetherCued == 1) .* (condition == 1) .* ((probeNeutral+probeForget) == 1)==1); % neutral item in remember cue trials
    ineuTargForg = target((whetherCued == 1) .* (condition == 2) .* (probeNeutral == 1)==1); % neutral item in forget cue trials
    iremTargRem = target((whetherCued == 1).*(condition == 1).*(probeRem == 1) == 1); % remember item in remember cue trials
    iremTargForg = target((whetherCued == 1).*(condition == 2).*(probeRem == 1) == 1); % remember item in forget cue trials
    ineuTargNeu = target(whetherCued == 0); % item in neutral cue trials

    ineuRespRem = error((whetherCued == 1) .* (condition == 1) .* ((probeNeutral+probeForget) == 1)==1); % neutral item in remember cue trials
    ineuRespForg = error((whetherCued == 1) .* (condition == 2) .* (probeNeutral == 1)==1); % neutral item in forget cue trials
    iremRespRem = error((whetherCued == 1).*(condition == 1).*(probeRem == 1) == 1); % remember item in remember cue trials
    iremRespForg = error((whetherCued == 1).*(condition == 2).*(probeRem == 1) == 1); % remember item in forget cue trials
    ineuRespNeu = error(whetherCued == 0); % item in neutral cue trials

    remError(isub,1) = mean(abs(ineuRespRem)); % neutral
    remError(isub,2) = mean(abs(iremRespRem)); % cued
    fogError(isub,1) = mean(abs(ineuRespForg)); % neutral
    fogError(isub,2) = mean(abs(iremRespForg)); % cued
    neuError(isub,1) = mean(abs(ineuRespNeu));

end

fogError1 = fogError(:,2);
remError1 = remError(:,2);

%% Experiment 2
dataPath = strcat(mainPath,'data/exp2/');
disp('Loading Data');
files = dir([dataPath expt '_*.csv']);
alpha = 1:360;
for isub = 1:numel(files)
    data = readtable([dataPath files(isub).name], 'Delimiter',',');
    data = data{:,:};
    error = data(:,20);
    data(error==999,:) = [];
    
    % simulate guess accuracy for each subjects
    target = data(:,17);
    error = data(:,20);

    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(target),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(target)))));
    end

    % remove guess trials
    sampError = sort(sampError);
    guessError(isub,1) = sampError(50);

    % data(error>guessError(isub),:) = [];
    target = data(:,17);
    pct(isub) = numel(target)/numel(error);

    
    error = data(:,20);
    whetherCued = data(:,6);
    condition = data(:,4);
    probeRem = data(:,14);
    probeNeutral = data(:,16);
    probeForget = data(:,15);

    iremErrorForg = error((whetherCued == 1).*(condition == 2).*(probeRem == 1) == 1);
    ineuErrorRem = error((whetherCued == 1) .* (condition == 1) .* ((probeNeutral+probeForget) == 1)==1);
    iremErrorRem = error((whetherCued == 1).*(condition == 1).*(probeRem == 1) == 1);
    ineuErrorForg = error((whetherCued == 1) .* (condition == 2) .* (probeNeutral == 1)==1);

    remError(isub,1) = mean(abs(ineuErrorRem)); % neutral
    remError(isub,2) = mean(abs(iremErrorRem)); % cued
    fogError(isub,1) = mean(abs(ineuErrorForg)); % neutral
    fogError(isub,2) = mean(abs(iremErrorForg)); % cued

end

fogError2 = fogError(:,2);
remError2 = remError(:,2);

%% Experiment 3
clear remError fogError guessError
%% two session subjects
% Load files
disp('Loading Data');
dataPath = strcat(mainPath,'data/exp3/');
files = dir([dataPath expt '_*.csv']);
alpha = 1:360;
for isub = 1:numel(files)/2
    data1 = readtable([dataPath files((isub-1)*2+1).name], 'Delimiter',',');
    data1 = data1{:,:};
    data2 = readtable([dataPath files((isub-1)*2+2).name], 'Delimiter',',');
    data2 = data2{:,:};
    data = [data1; data2];
    error = data(:,20);
    data(error==999,:) = [];

    % simulate guess accuracy for each subjects
    target = data(:,17);
    error = data(:,20);

    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(target),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(target)))));
    end

    % remove guess trials
    sampError = sort(sampError);
    guessError(isub,1) = sampError(50);

    % data(error>guessError(isub),:) = [];
    target = data(:,17);
    pct(isub) = numel(target)/numel(error);

    error = data(:,20);
    whetherCued = data(:,6);
    condition = data(:,4);
    probeRem = data(:,14);
    probeNeutral = data(:,16);
    probeForget = data(:,15);

    ineuErrorRem = error(((whetherCued == 1).*(condition == 1).*(probeNeutral+probeForget==1))==1); 
    iremErrorRem = error(((whetherCued == 1).*(condition == 1).*(probeRem==1))==1); 
    ineuErrorFog1 = error(((whetherCued == 1).*(condition == 2).*(probeNeutral==1))==1); 
    iremErrorFog1 = error(((whetherCued == 1).*(condition == 2).*(probeRem==1))==1); 
    iremErrorFog2 = error(((whetherCued == 1).*(condition == 3).*(probeRem==1))==1); 

    % error in remember item in remember cue and forget cue trials
    remError(isub,1) = mean(abs(ineuErrorRem)); % neutral
    remError(isub,2) = mean(abs(iremErrorRem)); % cued
    fogError31(isub,1) = mean(abs(ineuErrorFog1)); % neutral
    fogError31(isub,2) = mean(abs(iremErrorFog1)); % cued
    fogError32(isub,1) = mean(abs(iremErrorFog2)); % cued
end

exisingSub = size(remError,1);
% exisingSub = 0;
%% one session subjects
disp('Loading Data');
dataPath = strcat(mainPath,'data/exp3_1/');
files = dir([dataPath expt '_*.csv']);
alpha = 1:360;
for isub = 1:numel(files)
    data = readtable([dataPath files(isub).name], 'Delimiter',',');
    data = data{:,:};
    error = data(:,20);
    data(error==999,:) = [];

    % simulate guess accuracy for each subjects
    target = data(:,17);
    error = data(:,20);

    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(target),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(target)))));
    end

    % remove guess trials
    sampError = sort(sampError);
    guessError(exisingSub+isub,1) = sampError(50);

    data(error>guessError(exisingSub+isub),:) = [];
    target = data(:,17);
    pct(exisingSub+isub) = numel(target)/numel(error);

    whetherCued = data(:,6);
    condition = data(:,4);
    probeRem = data(:,14);
    probeNeutral = data(:,16);
    probeForget = data(:,15);

    ineuErrorRem = error(((whetherCued == 1).*(condition == 1).*(probeNeutral+probeForget==1))==1); 
    iremErrorRem = error(((whetherCued == 1).*(condition == 1).*(probeRem==1))==1); 
    ineuErrorFog1 = error(((whetherCued == 1).*(condition == 2).*(probeNeutral==1))==1); 
    iremErrorFog1 = error(((whetherCued == 1).*(condition == 2).*(probeRem==1))==1); 
    iremErrorFog2 = error(((whetherCued == 1).*(condition == 3).*(probeRem==1))==1); 

    % error in remember item in remember cue and forget cue trials
    remError(exisingSub+isub,1) = mean(abs(ineuErrorRem)); % neutral
    remError(exisingSub+isub,2) = mean(abs(iremErrorRem)); % cued
    fogError31(exisingSub+isub,1) = mean(abs(ineuErrorFog1)); % neutral
    fogError31(exisingSub+isub,2) = mean(abs(iremErrorFog1)); % cued
    fogError32(exisingSub+isub,1) = mean(abs(iremErrorFog2)); % cued
end

fogError3 = fogError31(:,2);
fogError4 = fogError32;
remError3 = remError(:,2);

%% experiment 2 control version

dataPath = strcat(mainPath,'data/exp4/');
disp('Loading Data');
files = dir([dataPath expt '_*.csv']);
alpha = 1:360;
for isub = 1:numel(files)
    data = readtable([dataPath files(isub).name], 'Delimiter',',');
    data = data{:,:};
    error = data(:,20);
    data(error==999,:) = [];
    
    % simulate guess accuracy for each subjects
    target = data(:,17);
    error = data(:,20);

    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(target),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(target)))));
    end

    % remove guess trials
    sampError = sort(sampError);
    guessError(isub,1) = sampError(50);

    % data(error>guessError(isub),:) = [];
    target = data(:,17);
    pct(isub) = numel(target)/numel(error);

    
    error = data(:,20);
    whetherCued = data(:,6);
    condition = data(:,4);
    probeRem = data(:,14);
    probeNeutral = data(:,16);
    probeForget = data(:,15);

    iremErrorForg = error((whetherCued == 1).*(condition == 2).*(probeRem == 1) == 1);
    ineuErrorRem = error((whetherCued == 1) .* (condition == 1) .* ((probeNeutral+probeForget) == 1)==1);
    iremErrorRem = error((whetherCued == 1).*(condition == 1).*(probeRem == 1) == 1);
    ineuErrorForg = error((whetherCued == 1) .* (condition == 2) .* (probeNeutral == 1)==1);

    remError(isub,1) = mean(abs(ineuErrorRem)); % neutral
    remError(isub,2) = mean(abs(iremErrorRem)); % cued
    fogError(isub,1) = mean(abs(ineuErrorForg)); % neutral
    fogError(isub,2) = mean(abs(iremErrorForg)); % cued

end

fogError5 = fogError(:,2);
remError4 = remError(:,2);
%% compare control with Experiments 1-3
[~,p_rem(1),~,stat_rem{1}] = ttest2(conError,remError1);
[~,p_rem(2),~,stat_rem{2}] = ttest2(conError,remError2);
[~,p_rem(3),~,stat_rem{3}] = ttest2(conError,remError3);
[~,p_rem(4),~,stat_rem{4}] = ttest2(conError,remError4);

[~,p_fog(1),~,stat_fog{1}] = ttest2(conError,fogError1);
[~,p_fog(2),~,stat_fog{2}] = ttest2(conError,fogError2);
[~,p_fog(3),~,stat_fog{3}] = ttest2(conError,fogError3);
[~,p_fog(4),~,stat_fog{4}] = ttest2(conError,fogError4);
[~,p_fog(5),~,stat_fog{5}] = ttest2(conError,fogError5);

d(1) = computeCohen_d(conError,fogError1, 'independent');
d(2) = computeCohen_d(conError,fogError2, 'independent');
d(3) = computeCohen_d(conError,fogError3, 'independent');
d(4) = computeCohen_d(conError,fogError4, 'independent');
d(5) = computeCohen_d(conError,fogError5, 'independent');

[BF(1),~] = bf.ttest2(conError,fogError1);
[BF(2),~] = bf.ttest2(conError,fogError2);
[BF(3),~] = bf.ttest2(conError,fogError3);
[BF(4),~] = bf.ttest2(conError,fogError4);
[BF(5),~] = bf.ttest2(conError,fogError5);