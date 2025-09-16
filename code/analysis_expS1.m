close all, clear all, clc

expt = 'forget_beh';

scriptpath = mfilename('fullpath');
scriptname = mfilename();

mainPath = scriptpath(1:length(scriptpath)-length(scriptname)-length('script')-1);
dataPath = strcat(mainPath,'data/data_control/');
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
    conError(isub,1) = mean(abs(error)); % neutral

    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(target),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(target)))));
    end

    guessError(isub,1) = min(sampError);


end

% exclude subjects
idx1 = find(conError(:,1) > guessError(:,1));
idx = unique(idx1);
conError(idx,:) = [];


%% Experiment 1
dataPath = strcat(mainPath,'data/v2/');
%% plot error distributionss
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

    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(ineuTargRem),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(ineuTargRem)))));
    end

    % guessError(isub,1) = min(sampError);
    sampError = sort(sampError);
    guessError(isub,1) = sampError(50);

    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(iremTargRem),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(iremTargRem)))));
    end

    % guessError(isub,2) = min(sampError);
    sampError = sort(sampError);
    guessError(isub,2) = sampError(50);

    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(ineuTargForg),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(ineuTargForg)))));
    end

    % guessError(isub,3) = min(sampError);
    sampError = sort(sampError);
    guessError(isub,3) = sampError(50);

    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(iremTargForg),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(iremTargForg)))));
    end

    % guessError(isub,4) = min(sampError);
    sampError = sort(sampError);
    guessError(isub,4) = sampError(50);


    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(ineuTargNeu),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(ineuTargNeu)))));
    end

    % guessError(isub,5) = min(sampError);
    sampError = sort(sampError);
    guessError(isub,5) = sampError(50);

    ineuRespRem = error((whetherCued == 1) .* (condition == 1) .* ((probeNeutral+probeForget) == 1)==1); % neutral item in remember cue trials
    ineuRespForg = error((whetherCued == 1) .* (condition == 2) .* (probeNeutral == 1)==1); % neutral item in forget cue trials
    iremRespRem = error((whetherCued == 1).*(condition == 1).*(probeRem == 1) == 1); % remember item in remember cue trials
    iremRespForg = error((whetherCued == 1).*(condition == 2).*(probeRem == 1) == 1); % remember item in forget cue trials
    remBlock = unique(blockIdx((whetherCued == 1).*(condition == 1)==1));
    fogBlock = unique(blockIdx((whetherCued == 1).*(condition == 2)==1));
    remIdx = [find(blockIdx==remBlock(1)); find(blockIdx==remBlock(2)); find(blockIdx==remBlock(3)); find(blockIdx==remBlock(4))];
    uncueIdx = find(whetherCued == 0);
    uncueRemIdx = intersect(remIdx, uncueIdx);
    fogIdx = [find(blockIdx==fogBlock(1)); find(blockIdx==fogBlock(2)); find(blockIdx==fogBlock(3)); find(blockIdx==fogBlock(4))];
    uncueIdx = find(whetherCued == 0);
    uncueFogIdx = intersect(fogIdx, uncueIdx);

    ineuRespNeuRem = error(uncueRemIdx); % item in neutral cue trials in remember block
    ineuRespNeuFog = error(uncueFogIdx); % item in neutral cue trials in forget block
    ineuRespNeu = error(whetherCued == 0); % item in neutral cue trials

    remError(isub,1) = mean(abs(ineuRespRem)); % neutral
    remError(isub,2) = mean(abs(iremRespRem)); % cued
    fogError(isub,1) = mean(abs(ineuRespForg)); % neutral
    fogError(isub,2) = mean(abs(iremRespForg)); % cued
    neuError(isub,1) = mean(abs(ineuRespNeu));
    % neuError(isub,1) = mean(abs(ineuRespNeuRem)); % remember cue
    % neuError(isub,2) = mean(abs(ineuRespNeuFog)); % forget cue

end

%%
idx1 = find(fogError(:,1) > guessError(:,3)); idx2 = find(fogError(:,2) > guessError(:,4)); idx3 = find(remError(:,1) > guessError(:,1)); idx4 = find(remError(:,2) > guessError(:,2)); idx5 = find(neuError > guessError(:,5)); 
idx = unique([idx1; idx2; idx3; idx4; idx5]);
remError(idx,:) = []; fogError(idx,:) = []; neuError(idx,:) = [];

fogError1 = fogError(:,2);
remError1 = remError(:,2);

%% Experiment 2
dataPath = strcat(mainPath,'data/v3/');
%% plot error distributions
% Load files
clear remError fogError guessError
disp('Loading Data');
files = dir([dataPath expt '_*.csv']);
alpha = 1:360;
for isub = 1:numel(files)
    data = readtable([dataPath files(isub).name], 'Delimiter',',');
    data = data{:,:};
    error = data(:,20);
    data(error==999,:) = [];
    % exclude trials with error > 50°
    % data(abs(error)> 60,:) = [];
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

    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(ineuTargRem),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(ineuTargRem)))));
    end

    sampError = sort(sampError);
    guessError(isub,1) = sampError(50);

    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(iremTargRem),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(iremTargRem)))));
    end

    sampError = sort(sampError);
    guessError(isub,2) = sampError(50);

    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(ineuTargForg),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(ineuTargForg)))));
    end

    sampError = sort(sampError);
    guessError(isub,3) = sampError(50);


    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(iremTargForg),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(iremTargForg)))));
    end

    sampError = sort(sampError);
    guessError(isub,4) = sampError(50);

    ineuErrorRem = error((whetherCued == 1) .* (condition == 1) .* ((probeNeutral+probeForget) == 1)==1); % uncued item in remember cue trial
    ineuErrorForg = error((whetherCued == 1) .* (condition == 2) .* (probeNeutral == 1)==1); % uncued item in forget cue trial
    iremErrorRem = error((whetherCued == 1).*(condition == 1).*(probeRem == 1) == 1); % remember item in remember cue trial
    iremErrorForg = error((whetherCued == 1).*(condition == 2).*(probeRem == 1) == 1); % remember item in forget cue trial


    remError(isub,1) = mean(abs(ineuErrorRem)); % neutral
    remError(isub,2) = mean(abs(iremErrorRem)); % cued
    fogError(isub,1) = mean(abs(ineuErrorForg)); % neutral
    fogError(isub,2) = mean(abs(iremErrorForg)); % cued

end

%% exclude subjects
idx1 = find(fogError(:,1) > guessError(:,3)); idx2 = find(fogError(:,2) > guessError(:,4)); idx3 = find(remError(:,1) > guessError(:,1)); idx4 = find(remError(:,2) > guessError(:,2)); 
idx = unique([idx1; idx2; idx3; idx4]);
remError(idx,:) = []; fogError(idx,:) = []; 

fogError2 = fogError(:,2);
remError2 = remError(:,2);

%% Experiment 3
clear remError fogError guessError
dataPath = strcat(mainPath,'data/v5/');
disp('Loading Data');
dataPath = strcat(mainPath,'data/v5/');
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
    target = data(:,17);
    error = data(:,20);
    whetherCued = data(:,6);
    condition = data(:,4);
    probeRem = data(:,14);
    probeNeutral = data(:,16);
    probeForget = data(:,15);
    blockIdx = data(:,3);
    blockNum = max(data(:,3));

    ineuTargRem = target(((whetherCued == 1).*(condition == 1).*(probeNeutral+probeForget==1))==1); 
    iremTargRem = target(((whetherCued == 1).*(condition == 1).*(probeRem==1))==1); 
    ineuTargFog1 = target(((whetherCued == 1).*(condition == 2).*(probeNeutral==1))==1); 
    iremTargFog1 = target(((whetherCued == 1).*(condition == 2).*(probeRem==1))==1); 
    iremTargFog2 = target(((whetherCued == 1).*(condition == 3).*(probeRem==1))==1); 

    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(ineuTargRem),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(ineuTargRem)))));
    end

    guessError(isub,1) = min(sampError);
    sampError = sort(sampError);
    guessError(isub,1) = sampError(50);

    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(iremTargRem),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(iremTargRem)))));
    end

    guessError(isub,2) = min(sampError);
    sampError = sort(sampError);
    guessError(isub,2) = sampError(50);

    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(ineuTargFog1),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(ineuTargFog1)))));
    end

    guessError(isub,3) = min(sampError);
    sampError = sort(sampError);
    guessError(isub,3) = sampError(50);


    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(iremTargFog1),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(iremTargFog1)))));
    end
    guessError(isub,4) = min(sampError);
    sampError = sort(sampError);
    guessError(isub,4) = sampError(50);

    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(iremTargFog2),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(iremTargFog2)))));
    end
    guessError(isub,5) = min(sampError);
    sampError = sort(sampError);
    guessError(isub,5) = sampError(50);


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
dataPath = strcat(mainPath,'data/v5_1/');
files = dir([dataPath expt '_*.csv']);
alpha = 1:360;
for isub = 1:numel(files)
    data = readtable([dataPath files(isub).name], 'Delimiter',',');
    data = data{:,:};
    error = data(:,20);
    data(error==999,:) = [];
    target = data(:,17);
    error = data(:,20);
    whetherCued = data(:,6);
    condition = data(:,4);
    probeRem = data(:,14);
    probeNeutral = data(:,16);
    probeForget = data(:,15);
    blockIdx = data(:,3);
    blockNum = max(data(:,3));

    ineuTargRem = target(((whetherCued == 1).*(condition == 1).*(probeNeutral+probeForget==1))==1); 
    iremTargRem = target(((whetherCued == 1).*(condition == 1).*(probeRem==1))==1); 
    ineuTargFog1 = target(((whetherCued == 1).*(condition == 2).*(probeNeutral==1))==1); 
    iremTargFog1 = target(((whetherCued == 1).*(condition == 2).*(probeRem==1))==1); 
    iremTargFog2 = target(((whetherCued == 1).*(condition == 3).*(probeRem==1))==1); 

    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(ineuTargRem),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(ineuTargRem)))));
    end

    guessError(exisingSub+isub,1) = min(sampError);
    sampError = sort(sampError);
    guessError(exisingSub+isub,1) = sampError(50);

    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(iremTargRem),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(iremTargRem)))));
    end

    guessError(exisingSub+isub,2) = min(sampError);
    sampError = sort(sampError);
    guessError(exisingSub+isub,2) = sampError(50);

    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(ineuTargFog1),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(ineuTargFog1)))));
    end

    guessError(exisingSub+isub,3) = min(sampError);
    sampError = sort(sampError);
    guessError(exisingSub+isub,3) = sampError(50);


    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(iremTargFog1),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(iremTargFog1)))));
    end
    guessError(exisingSub+isub,4) = min(sampError);
    sampError = sort(sampError);
    guessError(exisingSub+isub,4) = sampError(50);

    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(iremTargFog2),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(iremTargFog2)))));
    end
    guessError(exisingSub+isub,5) = min(sampError);
    sampError = sort(sampError);
    guessError(exisingSub+isub,5) = sampError(50);


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


%% exclude subjects

idx1 = find(remError(:,1) > guessError(:,1)); 
idx2 = find(remError(:,2) > guessError(:,2)); 
idx3 = find(fogError31(:,1) > guessError(:,3)); 
idx4 = find(fogError31(:,2) > guessError(:,4)); 
idx5 = find(fogError32 > guessError(:,5));
idx = unique([idx1; idx2; idx3; idx4; idx5]);
remError(idx,:) = []; fogError31(idx,:) = []; fogError32(idx,:) = [];


fogError3 = fogError31(:,2);
fogError4 = fogError32;
remError3 = remError(:,2);

%% compare control with Experiments 1-3
[~,p_rem(1),~,stat_rem{1}] = ttest2(conError,remError1);
[~,p_rem(2),~,stat_rem{2}] = ttest2(conError,remError2);
[~,p_rem(3),~,stat_rem{3}] = ttest2(conError,remError3);

[~,p_fog(1),~,stat_fog{1}] = ttest2(conError,fogError1);
[~,p_fog(2),~,stat_fog{2}] = ttest2(conError,fogError2);
[~,p_fog(3),~,stat_fog{3}] = ttest2(conError,fogError3);
[~,p_fog(4),~,stat_fog{4}] = ttest2(conError,fogError4);

d(1) = computeCohen_d(conError,fogError1, 'independent');
d(2) = computeCohen_d(conError,fogError2, 'independent');
d(3) = computeCohen_d(conError,fogError3, 'independent');
d(4) = computeCohen_d(conError,fogError4, 'independent');

[BF(1),~] = bf.ttest2(conError,fogError1);
[BF(2),~] = bf.ttest2(conError,fogError2);
[BF(3),~] = bf.ttest2(conError,fogError3);
[BF(4),~] = bf.ttest2(conError,fogError4);
%% helper functions
function data = exclue_outliers(data)
    while sum(data > mean(data)+2*std(data)) > 0
        data(data > mean(data)+2*std(data)) = [];
    end
end

function bootstrap_CI = get_bootstrap_CI(data, n)
    boots_data = NaN(size(data));
    for ii = 1:n
        idx = randi(numel(data),numel(data),1);
        boots_data(ii) = mean(data(idx));
    end
    boots_data = sort(boots_data);
    bootstrap_CI = [boots_data(round(n*0.025)), boots_data(round(n*0.975))];
end

function bootstrap_CI_all = get_bootstrap_CI_all(data, n)
    bootstrap_CI_all = NaN(2, size(data,2));
    for ii = 1:size(data,2)
        bootstrap_CI = get_bootstrap_CI(data(:,ii), n);
        bootstrap_CI_all(:,ii) = bootstrap_CI';
    end
end


function sig = get_significance(p)
    for i = 1:numel(p)
        if p(i) > 0.05
            sig{i} = ' ';
        elseif p(i) <= 0.05 && p(i) > 0.01
            sig{i} = '*';
        elseif p(i) <= 0.01 && p(i) > 0.001
            sig{i} = '**';
        elseif p(i) <= 0.001
            sig{i} = '***';
        end
    end
end