close all, clear all, clc

expt = 'forget_beh';

scriptpath = mfilename('fullpath');
scriptname = mfilename();

mainPath = scriptpath(1:length(scriptpath)-length(scriptname)-length('code')-1);
dataPath = strcat(mainPath,'data/exp2/');
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
    
    % simulate guess accuracy for each subjects
    target = data(:,17);
    error = data(:,20);

    errorAll(isub,1) = mean(abs(error));

    guess = 1/numel(alpha) * ones(1,numel(alpha));
    for i = 1:1000
        resp = randsample(alpha,numel(target),true,guess);
        sampError(i) = mean(abs(rad2deg(circ_dist(deg2rad(resp'), deg2rad(target)))));
    end

    % remove guess trials
    sampError = sort(sampError);
    guessError(isub,1) = sampError(50);

    % data(abs(error)>guessError(isub),:) = [];
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

%%
[~,p(1),~,stat{1}] = ttest(fogError(:,1),remError(:,1));
[~,p(2),~,stat{2}] = ttest(fogError(:,2),remError(:,2));

d(1) = computeCohen_d(fogError(:,1),remError(:,1), 'paired');
d(2) = computeCohen_d(fogError(:,2),remError(:,2), 'paired');
[BF(1),~] = bf.ttest(fogError(:,1),remError(:,1));
[BF(2),~] = bf.ttest(fogError(:,2),remError(:,2));

sig = get_significance(p);

diff_neu = remError(:,1) - fogError(:,1);
diff_rem = remError(:,2) - fogError(:,2);
[~,p_diff,~,stat_diff] = ttest(diff_neu, diff_rem);

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

figure(); set(gcf,'color','white', 'Units', 'Normalized', 'OuterPosition', [0 0 0.25 0.5], 'Visible', 'on'); hold on
% data = [remError(idx_rem,2), fogError(idx_rem,2)];  CI = get_bootstrap_CI_all(data,1000);
data = [remError(:,2), fogError(:,2)]; 
% data = Cousineau_correction(data);
CI = get_bootstrap_CI_all(data,1000);
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

line([1 2], [55 55], 'Color', 'k', 'LineWidth', 2); 
line([1 1], [54 55], 'Color', 'k', 'LineWidth', 2); 
line([2 2], [54 55], 'Color', 'k', 'LineWidth', 2); 
text(1.5, 60, ['p=' num2str(round(p(2),3))], 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 30);
text(1.5, 64, sig{2}, 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 40, 'FontWeight','bold');

ylabel('Response error(°)')
set(gca,'XLim',[0.5,2.5], 'XTick', 1:2, 'XTickLabel', {'Prioritized', 'Prioritized'}, 'TickDir', 'out', 'FontName','Arial','FontSize', 30, 'YLim', [0,80], 'LineWidth',2);

saveas(gca,fullfile(mainPath,'results/exp2', 'remember_item.svg'));
saveas(gca,fullfile(mainPath,'results/exp2', 'remember_item.png'));

figure(); set(gcf,'color','white', 'Units', 'Normalized', 'OuterPosition', [0 0 0.25 0.5], 'Visible', 'on'); hold on
data = [remError(:,1), fogError(:,1)]; 
% data = Cousineau_correction(data);
CI = get_bootstrap_CI_all(data,1000);
% data = [remError(idx_neu,1), fogError(idx_neu,1)]; CI = get_bootstrap_CI_all(data,1000);
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
text(1.5, 60, ['p=' num2str(round(p(1),3))], 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 30);
text(1.5, 64, sig{1}, 'HorizontalAlignment', 'center','FontName','Arial','FontSize', 40, 'FontWeight','bold');

ylabel('Response error(°)')
set(gca,'XLim',[0.5,2.5], 'XTick', 1:2, 'XTickLabel', {'Neutral', 'Neutral'}, 'TickDir', 'out', 'FontName','Arial','FontSize', 30, 'YLim', [0,80], 'LineWidth',2);

saveas(gca,fullfile(mainPath,'results/exp2', 'neutral_item.svg'));
saveas(gca,fullfile(mainPath,'results/exp2', 'neutral_item.png'));
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

function Z = Cousineau_correction(X)
    Y = X - repmat(mean(X'),size(X,2),1)' + mean(mean(X));
     J = size(X,2);
    Z = sqrt(J/(J-1)) .* ( Y' - repmat(mean(Y), size(Y,1),1)' )' + (repmat(mean(Y), size(Y,1),1)' )' ;
end