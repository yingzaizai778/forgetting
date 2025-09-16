close all, clear all, clc

expt = 'forget_beh';

scriptpath = mfilename('fullpath');
scriptname = mfilename();

mainPath = scriptpath(1:length(scriptpath)-length(scriptname)-length('script')-1);

addpath(genpath(fullfile(mainPath, 'toolbox' , 'CircStat2012a')));
addpath(genpath(fullfile(mainPath, 'toolbox' , 'MemToolbox-master')));

%% two session subjects
% Load files
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
    fogError1(isub,1) = mean(abs(ineuErrorFog1)); % neutral
    fogError1(isub,2) = mean(abs(iremErrorFog1)); % cued
    fogError2(isub,1) = mean(abs(iremErrorFog2)); % cued
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
    fogError1(exisingSub+isub,1) = mean(abs(ineuErrorFog1)); % neutral
    fogError1(exisingSub+isub,2) = mean(abs(iremErrorFog1)); % cued
    fogError2(exisingSub+isub,1) = mean(abs(iremErrorFog2)); % cued
end

%% exclude subjects

idx1 = find(remError(:,1) > guessError(:,1)); idx2 = find(remError(:,2) > guessError(:,2)); idx3 = find(fogError1(:,1) > guessError(:,3)); idx4 = find(fogError1(:,2) > guessError(:,4)); idx5 = find(fogError2 > guessError(:,5));
idx = unique([idx1; idx2; idx3; idx4; idx5]);
remError(idx,:) = []; fogError1(idx,:) = []; fogError2(idx,:) = [];


%%
[~,p_neutral,~,stat_neutral] = ttest(remError(:,1),fogError1(:,1));
d_neutral = computeCohen_d(fogError1(:,1),remError(:,1), 'paired');
sig_neutral = get_significance(p_neutral);
[BF_neutral,~] = bf.ttest(fogError1(:,1),remError(:,1));
[~,p_cued(1),~,stat_cued{1}] = ttest(remError(:,2),fogError1(:,2));
[~,p_cued(2),~,stat_cued{2}] = ttest(remError(:,2),fogError2);
[~,p_cued(3),~,stat_cued{3}] = ttest(fogError1(:,2),fogError2);
sig_cued = get_significance(p_cued);
d_cued(1) = computeCohen_d(remError(:,2),fogError1(:,2), 'paired');
d_cued(2) = computeCohen_d(remError(:,2),fogError2, 'paired');
d_cued(3) = computeCohen_d(fogError1(:,2),fogError2, 'paired');
[BF(1),~] = bf.ttest(remError(:,2),fogError1(:,2));
[BF(2),~] = bf.ttest(remError(:,2),fogError2);
[BF(3),~] = bf.ttest(fogError1(:,2),fogError2);
%%
figure(); set(gcf,'color','w', 'Units', 'Normalized', 'OuterPosition', [0 0 0.2 0.5], 'Visible', 'on'); hold on
data = [remError(:,1), fogError1(:,1)]; CI = get_bootstrap_CI_all(data,1000);
b1 = bar(1,mean(data(:,1)),0.6);
b1.FaceColor = '#21409A'; 
b1.EdgeColor = 'k';
b1.LineWidth = 2;

b2 = bar(2,mean(data(:,2)),0.6);
b2.FaceColor = '#21409A'; 
b2.EdgeColor = 'k';
b2.LineWidth = 2;

for isub = 1:size(data,1)
    plot(1:2,data(isub,:),'LineWidth',0.5,'Color',[0.7,0.7,0.7]);
end

errorbar(1:2, mean(data), mean(data)-CI(1,:), CI(2,:)-mean(data), 'CapSize',12,'LineWidth',2,'Color','k', 'LineStyle','none');

line([1 2], [55 55], 'Color', 'k', 'LineWidth', 2); 
line([1 1], [54 55], 'Color', 'k', 'LineWidth', 2); 
line([2 2], [54 55], 'Color', 'k', 'LineWidth', 2); 
text(1.5, 60, ['p=' num2str(round(p_neutral(1),3))], 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 30);
text(1.5, 64, sig_neutral{1}, 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 40, 'FontWeight','bold');

ylabel('Response error(°)')
set(gca,'XLim',[0.5,2.5], 'XTick', 1:2, 'XTickLabel', {'Neutral', 'Neutral'}, 'TickDir', 'out', 'FontName','Arial','FontSize', 30, 'YLim', [0,80], 'LineWidth',2);

saveas(gca,fullfile(mainPath,'results/v5', 'neutral_item1.svg'));
saveas(gca,fullfile(mainPath,'results/v5', 'neutral_item1.png'));

%%
figure(); set(gcf,'color','w', 'Units', 'Normalized', 'OuterPosition', [0 0 0.4 0.5], 'Visible', 'on'); hold on
data = [remError(:,2), fogError1(:,2), fogError2]; CI = get_bootstrap_CI_all(data,1000);
b1 = bar(1,mean(data(:,1)),0.6);
b1.FaceColor = '#009444'; 
b1.EdgeColor = 'k';
b1.LineWidth = 2;

b2 = bar(2,mean(data(:,2)),0.6);
b2.FaceColor = '#009444'; 
b2.EdgeColor = 'k';
b2.LineWidth = 2;

b2 = bar(3,mean(data(:,3)),0.6);
b2.FaceColor = '#009444'; 
b2.EdgeColor = 'k';
b2.LineWidth = 2;

for isub = 1:size(data,1)
    plot(1:3,data(isub,:),'LineWidth',0.5,'Color',[0.7,0.7,0.7]);
end

errorbar(1:3, mean(data), mean(data)-CI(1,:), CI(2,:)-mean(data), 'CapSize',12,'LineWidth',2,'Color','k', 'LineStyle','none');

line([1 1.9], [60 60], 'Color', 'k', 'LineWidth', 2); 
line([1 1], [59 60], 'Color', 'k', 'LineWidth', 2); 
line([1.9 1.9], [59 60], 'Color', 'k', 'LineWidth', 2); 
text(1.5, 64, ['p=' num2str(round(p_cued(1),3))], 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 30);
text(1.5, 67, sig_cued{1}, 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 40, 'FontWeight','bold');

line([2.1 3], [60 60], 'Color', 'k', 'LineWidth', 2); 
line([2.1 2.1], [59 60], 'Color', 'k', 'LineWidth', 2); 
line([3 3], [59 60], 'Color', 'k', 'LineWidth', 2); 
text(2.5, 64, ['p=' num2str(round(p_cued(3),3))], 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 30);
text(2.5, 67, sig_cued{3}, 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 40, 'FontWeight','bold');

line([1 3], [70 70], 'Color', 'k', 'LineWidth', 2); 
line([1 1], [69 70], 'Color', 'k', 'LineWidth', 2); 
line([3 3], [69 70], 'Color', 'k', 'LineWidth', 2); 
text(2, 74, ['p=' num2str(round(p_cued(2),3))], 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 30);
text(2, 77, sig_cued{2}, 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 40, 'FontWeight','bold');

ylabel('Response error(°)')
set(gca,'XLim',[0.5,3.5], 'XTick', 1:3, 'XTickLabel', {'Prioritized', 'Prioritized', 'Prioritized'}, 'TickDir', 'out', 'FontName','Arial','FontSize', 30, 'YLim', [0,80], 'LineWidth',2);

saveas(gca,fullfile(mainPath,'results/v5', 'cued_item1.svg'));
saveas(gca,fullfile(mainPath,'results/v5', 'cued_item1.png'));
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