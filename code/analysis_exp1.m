close all, clear all, clc

expt = 'forget_beh';

scriptpath = mfilename('fullpath');
scriptname = mfilename();

mainPath = scriptpath(1:length(scriptpath)-length(scriptname)-length('script')-1);
dataPath = strcat(mainPath,'data/v2/');
addpath(genpath(fullfile(mainPath, 'toolbox' , 'CircStat2012a')));
addpath(genpath(fullfile(mainPath, 'toolbox' , 'MemToolbox-master')));
addpath(genpath(fullfile(mainPath, 'toolbox' , 'github_repo')));
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

%%

[~,p_rem(1),~,stat_rem{1}] = ttest(remError(:,1),neuError); % neutral vs neutral
[~,p_rem(2),~,stat_rem{2}] = ttest(remError(:,2),neuError); % cued vs neutral
[~,p_rem(3),~,stat_rem{3}] = ttest(remError(:,1),remError(:,2)); % cued vs neutral

d_rem(1) = computeCohen_d(remError(:,1),neuError, 'paired');
d_rem(2) = computeCohen_d(remError(:,2),neuError, 'paired');
d_rem(3) = computeCohen_d(remError(:,1),remError(:,2), 'paired');
[bf_rem(1),~] = bf.ttest(remError(:,1),neuError);
[bf_rem(2),~] = bf.ttest(remError(:,2),neuError);
[bf_rem(3),~] = bf.ttest(remError(:,1),remError(:,2));

sig_rem = get_significance(p_rem);

[~,p_fog(1),~,stat_fog{1}] = ttest(fogError(:,1),neuError);
[~,p_fog(2),~,stat_fog{2}] = ttest(fogError(:,2),neuError);
[~,p_fog(3),~,stat_fog{3}] = ttest(fogError(:,1),fogError(:,2));

d_fog(1) = computeCohen_d(fogError(:,1),neuError, 'paired');
d_fog(2) = computeCohen_d(fogError(:,2),neuError, 'paired');
d_fog(3) = computeCohen_d(fogError(:,1),fogError(:,2), 'paired');
[bf_fog(1),~] = bf.ttest(fogError(:,1),neuError);
[bf_fog(2),~] = bf.ttest(fogError(:,2),neuError);
[bf_fog(3),~] = bf.ttest(fogError(:,1),fogError(:,2));

sig_fog = get_significance(p_fog);

[~,p(1),~,stat{1}] = ttest(fogError(:,1),remError(:,1));
[~,p(2),~,stat{2}] = ttest(fogError(:,2),remError(:,2));

d(1) = computeCohen_d(fogError(:,1),remError(:,1), 'paired');
d(2) = computeCohen_d(fogError(:,2),remError(:,2), 'paired');
[BF(1),~] = bf.ttest(fogError(:,1),remError(:,1));
[BF(2),~] = bf.ttest(fogError(:,2),remError(:,2));

sig = get_significance(p);

[~,p_inter] = ttest(fogError(:,1)-remError(:,1), fogError(:,2)-remError(:,2));

% anova
data = [fogError(:,1), remError(:,1), fogError(:,2), remError(:,2)];
varNames = {'A1_B1','A1_B2','A2_B1','A2_B2'};
t = array2table(data, 'VariableNames', varNames);

within = table();
within.A = categorical([1 1 2 2]');
within.B = categorical([1 2 1 2]');
rm = fitrm(t, 'A1_B1-A2_B2 ~ 1', 'WithinDesign', within);

ranovaResults = ranova(rm, 'WithinModel', 'A*B');


%%
close all
figure(); set(gcf,'color','w', 'Units', 'Normalized', 'OuterPosition', [0 0 0.4 0.5], 'Visible', 'on'); hold on
data = [remError, neuError]; CI = get_bootstrap_CI_all(data,1000);
b1 = bar(1,mean(data(:,1)),0.6);
b1.FaceColor = '#21409A'; 
b1.EdgeColor = 'k';
b1.LineWidth = 2;

b2 = bar(2,mean(data(:,2)),0.6);
b2.FaceColor = '#009444'; 
b2.EdgeColor = 'k';
b2.LineWidth = 2;

b3 = bar(3,mean(data(:,3)),0.6);
b3.FaceColor = [0.5,0.5,0.5]; 
b3.EdgeColor = 'k';
b3.LineWidth = 2;

for isub = 1:size(data,1)
    plot(1:3,data(isub,:),'LineWidth',0.5,'Color',[0.7,0.7,0.7]);
end

errorbar(1:3, mean(data), mean(data)-CI(1,:), CI(2,:)-mean(data), 'CapSize',12,'LineWidth',2,'Color','k', 'LineStyle','none');


% neutral vs baseline
line([1 3], [72 72], 'Color', 'k', 'LineWidth', 2); 
line([1 1], [71 72], 'Color', 'k', 'LineWidth', 2); 
line([3 3], [71 72], 'Color', 'k', 'LineWidth', 2); 
text(2, 75, ['p=' num2str(round(p_rem(1),3))], 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 30);
text(2, 78, sig_rem{1}, 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 40, 'FontWeight','bold');

% cued vs baseline
line([2.1 3], [60 60], 'Color', 'k', 'LineWidth', 2); 
line([2.1 2.1], [59 60], 'Color', 'k', 'LineWidth', 2); 
line([3 3], [59 60], 'Color', 'k', 'LineWidth', 2); 
text(2.5, 63, ['p=' num2str(round(p_rem(2),3))], 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 30);
text(2.5, 66, sig_rem{2}, 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 40, 'FontWeight','bold');

% cued vs neutral
line([1 1.9], [60 60], 'Color', 'k', 'LineWidth', 2); 
line([1.9 1.9], [59 60], 'Color', 'k', 'LineWidth', 2); 
line([1 1], [59 60], 'Color', 'k', 'LineWidth', 2); 
text(1.5, 63, ['p=' num2str(round(p_rem(3),3))], 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 30);
text(1.5, 66, sig_rem{3}, 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 40, 'FontWeight','bold');



ylabel('Response error(°)')
set(gca,'XLim',[0.5,3.5], 'XTick', 1:3, 'XTickLabel', {'Neutral', 'Prioritized', 'Baseline'}, 'TickDir', 'out', 'FontName','Arial','FontSize', 30, 'YLim', [0,80], 'LineWidth',2);

saveas(gca,fullfile(mainPath,'results/v2', 'remember_cue.svg'));
saveas(gca,fullfile(mainPath,'results/v2', 'remember_cue.png'));
%%
figure(); set(gcf,'color','w', 'Units', 'Normalized', 'OuterPosition', [0 0 0.4 0.5], 'Visible', 'on'); hold on
data = [fogError, neuError]; CI = get_bootstrap_CI_all(data,1000);
b1 = bar(1,mean(data(:,1)),0.6);
b1.FaceColor = '#21409A'; 
b1.EdgeColor = 'k';
b1.LineWidth = 2;

b2 = bar(2,mean(data(:,2)),0.6);
b2.FaceColor = '#009444';  
b2.EdgeColor = 'k';
b2.LineWidth = 2;

b3 = bar(3,mean(data(:,3)),0.6);
b3.FaceColor = [0.5,0.5,0.5]; 
b3.EdgeColor = 'k';
b3.LineWidth = 2;

for isub = 1:size(data,1)
    plot(1:3,data(isub,:),'LineWidth',0.5,'Color',[0.7,0.7,0.7]);
end

errorbar(1:3, mean(data), mean(data)-CI(1,:), CI(2,:)-mean(data), 'CapSize',12,'LineWidth',2,'Color','k', 'LineStyle','none');

line([1 3], [72 72], 'Color', 'k', 'LineWidth', 2); 
line([1 1], [71 72], 'Color', 'k', 'LineWidth', 2); 
line([3 3], [71 72], 'Color', 'k', 'LineWidth', 2); 
text(2, 75, ['p=' num2str(round(p_fog(1),3))], 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 30);
text(2, 78, sig_fog{1}, 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 40, 'FontWeight','bold');

line([2.1 3], [60 60], 'Color', 'k', 'LineWidth', 2); 
line([2.1 2.1], [59 60], 'Color', 'k', 'LineWidth', 2); 
line([3 3], [59 60], 'Color', 'k', 'LineWidth', 2); 
text(2.5, 63, ['p=' num2str(round(p_fog(2),3))], 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 30);
text(2.5, 66, sig_fog{2}, 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 40, 'FontWeight','bold');


line([1 1.9], [60 60], 'Color', 'k', 'LineWidth', 2); 
line([1.9 1.9], [59 60], 'Color', 'k', 'LineWidth', 2); 
line([1 1], [59 60], 'Color', 'k', 'LineWidth', 2); 
text(1.5, 63, ['p=' num2str(round(p_fog(3),3))], 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 30);
text(1.5, 66, sig_fog{3}, 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 40, 'FontWeight','bold');


ylabel('Response error(°)')
set(gca,'XLim',[0.5,3.5], 'XTick', 1:3, 'XTickLabel', {'Neutral', 'Prioritized', 'Baseline'}, 'TickDir', 'out', 'FontName','Arial','FontSize', 30, 'YLim', [0,80], 'LineWidth',2);

saveas(gca,fullfile(mainPath,'results/v2', 'forget_cue.svg'));
saveas(gca,fullfile(mainPath,'results/v2', 'forget_cue.png'));
%%

figure(); set(gcf,'color','w', 'Units', 'Normalized', 'OuterPosition', [0 0 0.25 0.5], 'Visible', 'on'); hold on
data = [remError(:,2), fogError(:,2)]; CI = get_bootstrap_CI_all(data,1000);
b1 = bar(1,mean(data(:,1)),0.6);
b1.FaceColor = '#009444'; 
b1.EdgeColor = 'k';
b1.LineWidth = 2;

b2 = bar(2,mean(data(:,2)),0.6);
b2.FaceColor = '#009444'; 
b2.EdgeColor = 'k';
b2.LineWidth = 2;

for isub = 1:size(data,1)
    plot(1:2,data(isub,:),'LineWidth',0.5,'Color',[0.7,0.7,0.7]);
end

errorbar(1:2, mean(data), mean(data)-CI(1,:), CI(2,:)-mean(data), 'CapSize',12,'LineWidth',2,'Color','k', 'LineStyle','none');

line([1 2], [60 60], 'Color', 'k', 'LineWidth', 2); 
line([1 1], [59 60], 'Color', 'k', 'LineWidth', 2); 
line([2 2], [59 60], 'Color', 'k', 'LineWidth', 2); 
text(1.5, 63, ['p=' num2str(round(p(2),3))], 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 30);
text(1.5, 66, sig{2}, 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 40, 'FontWeight','bold');

ylabel('Response error(°)')
set(gca,'XLim',[0.5,2.5], 'XTick', 1:2, 'XTickLabel', {'Prioritized', 'Prioritized'}, 'TickDir', 'out', 'FontName','Arial','FontSize', 30, 'YLim', [0,80], 'LineWidth',2);

saveas(gca,fullfile(mainPath,'results/v2', 'remember_item.svg'));
saveas(gca,fullfile(mainPath,'results/v2', 'remember_item.png'));

figure(); set(gcf,'color','w', 'Units', 'Normalized', 'OuterPosition', [0 0 0.25 0.5], 'Visible', 'on'); hold on
data = [remError(:,1), fogError(:,1)]; CI = get_bootstrap_CI_all(data,1000);
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

line([1 2], [60 60], 'Color', 'k', 'LineWidth', 2); 
line([1 1], [59 60], 'Color', 'k', 'LineWidth', 2); 
line([2 2], [59 60], 'Color', 'k', 'LineWidth', 2); 
text(1.5, 63, ['p=' num2str(round(p(1),3))], 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 30);
text(1.5, 66, sig{1}, 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 40, 'FontWeight','bold');

ylabel('Response error(°)')
set(gca,'XLim',[0.5,2.5], 'XTick', 1:2, 'XTickLabel', {'Neutral', 'Neutral'}, 'TickDir', 'out', 'FontName','Arial','FontSize', 30, 'YLim', [0,80], 'LineWidth',2);

saveas(gca,fullfile(mainPath,'results/v2', 'neutral_item.svg'));
saveas(gca,fullfile(mainPath,'results/v2', 'neutral_item.png'));
%%
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