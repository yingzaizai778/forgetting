%% Notes:

% Cue type 1: Neutral and remember cues  = 1/3 of all trials in each block.
       % remember cues 2/3 (controled by cue_validity)
       % neutral cues 1/3
% Cue type 2: Neutral, remember and forget cues = 2/3 of all trials in each block.
       % remember cues 2/3 (controled by cue_validity)
       % neutral cues 1/6
       % forget cues 1/6
% Final overal counts are printed in the command window at the beginning of the actual trials.
% Actual trial overal time is recorded with tic toc function and time is printed at the end of the experiment.
% Conditions in actual trials are shuffled blockwise for better distribution.
% Halfway through each block there is a break (so if we have 90 trials / block then the participants get a break after 45 trials every time).
% For 450 trials we need 5 blocks of 90 trials (this is the only option that preserves the distribution).
% For debugging: the minimal number of trials per block must be 18 (except prac_trials)!! Otherwise the distribution will be off and it will crash.
% (I didn't use round in actual trials because it interfeared with the distribution but the practice trials have round so those wont crash if the trial number is off as long as there are at least 2 prac trials).



%%

sca; close all; clear all;

scriptpath = mfilename('fullpath');
scriptname = mfilename();

mainPath = scriptpath(1:length(scriptpath)-length(scriptname));
addpath(genpath([mainPath 'CircStat2012a']));
dataPath = strcat(mainPath,'dataD/');
rng shuffle;

if exist(dataPath,'dir')==0
    mkdir dataD
end

expt = 'forget_beh';
debug = 0;

sel_font = 'Courier New';

maxBonus = 25;
%% Logfile set-up
cd(dataPath);
input = 0; quit = 0;
while input == 0 && quit == 0
    prompt = {'Subject Code:', 'Session:', 'Handedness:', 'Gender:', 'Age:'};
    titlen = 'Welcome! Please get the experimenter to input the following:';
    dims = [1 100];
    definput = {'101', '1', 'right', 'female', '20'};
    answer = inputdlg(prompt, titlen, dims, definput);
    
    try
        subject = answer{1};
        session = answer{2};
        hand = answer{3};
        gender = answer{4};
        age = answer{5};
        input = 1;
    catch
        quit = 1;
        disp('Error with input. Try again.');
    end
end

filenamed = strcat('demo_',expt,'_',num2str(subject),'_',datestr(now,'ddmmmyy_HHMM'));
csvnamed = strcat(filenamed,'.csv');
fileIDd = fopen(csvnamed,'w');
headerString = 'Subject,Session,Expt,Gender,Hand,Age\n';
fprintf(fileIDd,headerString);
demostring = strcat(subject, ',', session, ',', expt, ',', gender, ',', hand, ',', age);
fprintf(fileIDd,demostring);
fclose(fileIDd);

filename = strcat(expt,'_',num2str(subject),'_',datestr(now,'ddmmmyy_HHMM'));
csvname = strcat(filename,'.csv');
bevID = fopen(csvname,'w');
headerString = 'Subject,Session,Block,Trial,Condition,midStim,leftStim,rightStim,remCueLoc,forgetCueLoc,neutralCueLoc1,neutralCueLoc2,probeCueLoc,probeRem,probeForget,probeNeutral,Target,Estimate,RT,Error,Points\n';
fprintf(bevID, headerString);
jitter = 5;

%% Expt Parameters
if debug == 0
    screenNumber = 0;
    nTrials = 18; % trials per block 90 % 18 if debug
    nBlock = 2; % 5 % 2 if debug
    numPrac = 2;
    tfix = 500;
    tstim = 150; %ms
    tmask = 500;
    tprecue = 500;
    tcue1 = 1500;
    tcue2 = 1500;
    tdelay = 2500;
    tresp = 100000;
    totalTrials = nTrials * nBlock;
else
    screenNumber = 0;
    nTrials = 60; % trials per block
    nBlock = 8;
    numPrac = 6;
    tfix = 2;
    tstim = 2; %ms
    tmask = 2;
    tprecue = 2;
    tcue1 = 2;
    tcue2 = 2;
    tdelay = 2;
    tresp = 2;
    totalTrials = nTrials * nBlock;
end

black = [0, 0, 0];  
white = [1, 1, 1];  
gray = [0.5 0.5 0.5];
red = [1, 0, 0];
green = [0, 1, 0];
blue = [0, 0, 1];
fixSize = 24; % doubled the size of fixation dot
orientation = 15:30:360;
cue_validity = 2/3;
cue_percent = 1;

% If changing the code pay attention that in the practice trials these are writen with low case, such as rt, estimate etc.
RT = ones(nTrials * nBlock, 1) * tresp;
Estimate = NaN(nTrials * nBlock, 1);
Error = 999 * ones(nTrials * nBlock, 1);
Points = zeros(nTrials * nBlock, 1);
Target = 999 * ones(nTrials * nBlock, 1);

%% Disply setup
% get screen parameters
PsychDefaultSetup(2);
if debug == 0
    Screen('Preference', 'SkipSyncTests', 1);
else
    Screen('Preference', 'SkipSyncTests', 2);
end

screens = Screen('Screens');
[screenXpixels, screenYpixels] = Screen('WindowSize',screenNumber); %get x and y pixels of screen
[screenWidth, screenHeight] = Screen('DisplaySize',screenNumber); %get screen width and height in mm
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black, [0,0,screenXpixels,screenYpixels], [], 2, [], [], kPsychNeed32BPCFloat);% PsychImaging('OpenWindow', screenNumber, black, [0,0,screenXpixels,screenYpixels]); % open screen
[xCenter, yCenter] = RectCenter(windowRect);
topPriorityLevel = MaxPriority(window); % good idea for actual experiments

% get deg to pixel
viewDist = 57; %the distance from the screen, in cm
screenWidth = screenWidth/10; %mm to cm
screenHeight = screenHeight/10; %mm to cm
Hperdegree = viewDist * tand(1); %height for one degree in cm
Wperdegree = viewDist * tand(1); %width for one degree in cm
pixWidth = screenWidth/screenXpixels; %cm/pixel
pixHeight = screenHeight/screenYpixels; %cm/pixel
deg_width = atand(screenWidth/2 / viewDist) * 2;
deg_height = atand(screenHeight/2 / viewDist) * 2;
degToPixX = screenXpixels / deg_width;
degToPixY = screenYpixels / deg_height;

% prepare frame numbers for every event
frameRate = round(FrameRate(window));
avgFrame = 1000/frameRate;
ffix = round(tfix/avgFrame);
fstim = round(tstim/avgFrame);
fmask = round(tmask/avgFrame);
fprecue = round(tprecue/avgFrame);
fcue1 = round(tcue1/avgFrame);
fcue2 = round(tcue2/avgFrame);
fdelay = round(tdelay/avgFrame);
fresp = round(tresp/avgFrame);

% get the center of three stimuli
StimHeight = 3;
StimHeightPix = StimHeight*degToPixX;
StimWidth = 1;
StimWidthPix = StimWidth*degToPixX;
stimCenterDist = 5;
stimCenterDistPix = stimCenterDist*degToPixX;
stimCenter = [xCenter, yCenter-stimCenterDistPix; ...
    xCenter-stimCenterDistPix* cos(1/6*pi), yCenter+stimCenterDistPix* sin(1/6*pi); ...
    xCenter+stimCenterDistPix * cos(1/6*pi), yCenter+stimCenterDistPix * sin(1/6*pi)];

% prepare cue
cueDeg = 1;
cueSize = cueDeg*degToPixX;
cueLoc = [xCenter, yCenter, xCenter, yCenter-cueSize; ...   % up middle
    xCenter, yCenter, xCenter-cueSize * cos(1/6*pi), yCenter+cueSize * sin(1/6*pi); ... % left
    xCenter, yCenter, xCenter+cueSize * cos(1/6*pi), yCenter+cueSize * sin(1/6*pi)]; % right
cueWidth = 4;

% open screen
Screen('BlendFunction', window, 'GL_ONE', 'GL_ZERO'); % for alpha channel
ShowCursor(0);
startTime = GetSecs;
exptStart = Screen('Flip', window);  
Screen('Preference', 'TextRenderer', 1);
Screen('TextSize', window, 25);
Screen('TextFont', window, sel_font, 1); % 1 = bold
escapeKey = KbName('ESCAPE');
Priority(topPriorityLevel);



%% Define the table columns
% Define table columns using numeric values for colors (9 VARIATION, blue-green variations are repeated twice to make equal number of trials for each version)
% 1 = blue, 1 = red, 3 = green

first_cue_position1 = [1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1]';
first_cue_position2 = [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1]';
first_cue_position3 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2]';
second_cue_position1 = [1, 3, 1, 1, 3, 1, 2, 2, 3, 1, 1, 3, 2, 2, 3, 1, 1, 3]';
second_cue_position2 = [3, 1, 1, 3, 1, 1, 3, 1, 2, 2, 3, 1, 3, 1, 2, 2, 3, 1]';
second_cue_position3 = [1, 1, 3, 1, 1, 3, 1, 3, 1, 3, 2, 2, 1, 3, 1, 3, 2, 2]';
% cue_angle_pos1 = {'90', '90', '90', '90', '90', '90', '90', '90', '90', '90', '90', '90', '90', '90', '90', '90', '90', '90'}';
% cue_angle_pos2 = {'210', '210', '210', '210', '210', '210', '210', '210', '210', '210', '210', '210', '210', '210', '210', '210', '210', '210'}';
% cue_angle_pos3 = {'330', '330', '330', '330', '330', '330', '330', '330', '330', '330', '330', '330', '330', '330', '330', '330', '330', '330'}';
cue_angle_pos1 = repmat(90,[1,18])';
cue_angle_pos2 = repmat(210,[1,18])';
cue_angle_pos3 = repmat(330,[1,18])';
% Create a table
color_table = table(first_cue_position1, first_cue_position2, first_cue_position3, second_cue_position1, second_cue_position2, second_cue_position3, cue_angle_pos1, cue_angle_pos2, cue_angle_pos3);

%% Separate Trials into Two Tables
% Identify trials suitable for Cue Type 1
% Cue Type 1 trials have second_cue_colors containing only blue (1) and green (3)
type1_trials = color_table(all(ismember(color_table{:, {'second_cue_position1', 'second_cue_position2', 'second_cue_position3'}}, [1, 3]), 2), :);

% Identify trials suitable for Cue Type 2
% Cue Type 2 trials have second_cue_colors containing blue (1), green (3), and red (2)
type2_trials = color_table(any(color_table{:, {'second_cue_position1', 'second_cue_position2', 'second_cue_position3'}} == 2, 2), :);

% Ensure we have enough trials in each table
total_nCueType1 = nBlock * (nTrials / 3);
total_nCueType2 = nBlock * (nTrials - nTrials / 3);

if height(type1_trials) < total_nCueType1
    repeat_factor = ceil(total_nCueType1 / height(type1_trials));
    type1_trials = repmat(type1_trials, repeat_factor, 1);
end

if height(type2_trials) < total_nCueType2
    repeat_factor = ceil(total_nCueType2 / height(type2_trials));
    type2_trials = repmat(type2_trials, repeat_factor, 1);
end

% Shuffle the trials in each table
type1_trials = type1_trials(randperm(height(type1_trials)), :);
type2_trials = type2_trials(randperm(height(type2_trials)), :);

% Initialize an empty table for all blocks
all_blocks_trials = table();

for block = 1:nBlock
    % Determine quotas for this block
    nCueType1 = nTrials / 3;
    nCueType2 = nTrials - nCueType1;

    quota_type1_remember = nCueType1 * 2 / 3;
    quota_type1_neutral = nCueType1 - quota_type1_remember;

    quota_type2_remember = nCueType2 * 2 / 3;
    quota_type2_forget = nCueType2 * 1 / 6;
    quota_type2_neutral = nCueType2 - (quota_type2_remember + quota_type2_forget);

    % Randomly select trials for Cue Type 1
    idx_type1 = randperm(height(type1_trials), nCueType1);
    block_type1_trials = type1_trials(idx_type1, :);
    type1_trials(idx_type1, :) = []; % Remove selected trials

    % Randomly select trials for Cue Type 2
    idx_type2 = randperm(height(type2_trials), nCueType2);
    block_type2_trials = type2_trials(idx_type2, :);
    type2_trials(idx_type2, :) = []; % Remove selected trials

    % Combine trials for this block
    block_trials = [block_type1_trials; block_type2_trials];

    % Assign cue types
    block_trials.CueType = [ones(nCueType1, 1); ones(nCueType2, 1) * 2];

    % Shuffle the block trials
    block_trials = block_trials(randperm(height(block_trials)), :);

    % Append to all_blocks_trials
    all_blocks_trials = [all_blocks_trials; block_trials];
end

% Assign to shuffled_table_actual
shuffled_table_actual = all_blocks_trials;

% Update totalTrials
totalTrials = height(shuffled_table_actual);

% Ensure shuffled_table_prac has enough rows
if height(color_table) < numPrac
    repeat_factor = ceil(numPrac / height(color_table));
    shuffled_table_prac = repmat(color_table, repeat_factor, 1);
else
    shuffled_table_prac = color_table;  % If enough rows, just assign
end
% Shuffle the table and select numPrac rows
shuffled_table_prac = shuffled_table_prac(randperm(height(shuffled_table_prac), numPrac), :);  % Trim and shuffle


disp(shuffled_table_actual);
fprintf('Size of shuffled_table_actual: %d rows\n', height(shuffled_table_actual));

%% Main

% Call the genPracticeTrials function
[prac_probe_loc, prac_test_remember, prac_test_forget, prac_test_neutral] = genPracticeTrials(numPrac, shuffled_table_prac, cue_validity);

% Initialize block_no for practice loop
block_no = 0;

% Practice loop
for trial = 1:numPrac
    if mod(trial, numPrac) == 1
        block_no = block_no + 1;
        practice = 1;
        rt = ones(numPrac, 1) * tresp;
        estimate = NaN(numPrac, 1);
        error = 999 * ones(numPrac, 1);
        points = zeros(numPrac, 1);
        target = 999 * ones(numPrac, 1);
        all_Orientation = zeros(numPrac, 3); %small all_ for parc, large All_ for actual trials
        begin = 0; 
        ShowCursor('Arrow'); 
        mouseB = [0, 0]; % Initialize mouseB

        % Ensure mouse button is not pressed before displaying the message
        while any(mouseB) % Start of while (mouse button not pressed)
            [~, ~, mouseB] = GetMouse(window);
        end  % End of while (mouse button not pressed)

        while begin == 0 && quit == 0 % Start of the "begin task" loop
            DrawFormattedText(window, 'Your task is to draw the orientation \nof the cued triangle', 'center', 150, white);
            DrawFormattedText(window, 'GREEN ITEM IS TWICE LIKELY TO BE TESTED,\nRED ITEM DOES NOT NEED RESPONSE AND CAN BE FORGOTTEN', 'center', 250, red);
            DrawFormattedText(window, 'Click on the middle to continue', 'center', 350, white);
            DrawFormattedText(window, 'You will now begin with some practice trials', 'center', 450, white);
            Screen('FillOval', window, white, makeBox(windowRect(3) / 2, windowRect(4) / 2, fixSize));
            Screen('DrawingFinished', window); 
            Screen('Flip', window);

            while ~any(mouseB) % Start of waiting for a press
                [mouseX, mouseY, mouseB] = GetMouse(window);
                centerDist = ((mouseX - windowRect(3) / 2) ^ 2 + (mouseY - windowRect(4) / 2) ^ 2) ^ 0.5;
                if centerDist < fixSize / 2 && mouseB(1) == 1
                    begin = 1;
                end

                [keyIsDown, secs, keyCode] = KbCheck;
                if keyCode(escapeKey)
                    quit = 1;
                    if bevID > 0
                        fclose(bevID);
                    end    
                    cd(mainPath);
                    sca;
                    Priority(0);
                end
            end  % End of waiting for a press

            while any(mouseB) % Ensure mouse button is released
                [~, ~, mouseB] = GetMouse(window);
            end
        end  % End of the "begin task" loop
    end  % End of if block_no == 1
end  % End of if mod(trial, numPrac) == 1

% Practice trials loop
for prac_trial = 1:numPrac
    
    %% Fixation
    frame = 0; 
    HideCursor();
    while (quit == 0) && (frame < ffix) % Start of fixation loop
        frame = frame + 1;
        Screen('FillArc', window, gray, [xCenter - yCenter, 0, xCenter + yCenter, 2 * yCenter], 0, 360);
        drawFix(black, fixSize, windowRect, window);
        Screen('DrawingFinished', window);
        vbl = Screen('Flip', window);
        [keyIsDown, secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            quit = 1;
        end
    end  % End of fixation loop

    %% Stimulus
    index_Gabor = randperm(length(orientation));
    all_Orientation(prac_trial, :) = [orientation(index_Gabor(1)) + jitter * rand(), orientation(index_Gabor(2)) + jitter * rand(), orientation(index_Gabor(3)) + jitter * rand()];
    frame = 0;
    while (quit == 0) && (frame < fstim) % Start of stimulus loop
        frame = frame + 1;
        Screen('FillArc', window, gray, [xCenter - yCenter, 0, xCenter + yCenter, 2 * yCenter], 0, 360);
        drawFix(black, fixSize, windowRect, window);
        for ii = 1:3
            iCenter = stimCenter(ii, :);
            iangle = deg2rad(all_Orientation(prac_trial, ii));
            drawPoints = [iCenter(1) + StimHeightPix / 2 * sin(iangle), iCenter(2) - StimHeightPix / 2 * cos(iangle); ...
                          iCenter(1) - StimHeightPix / 2 * sin(iangle) - StimWidthPix / 2 * cos(iangle), iCenter(2) + StimHeightPix / 2 * cos(iangle) - StimWidthPix / 2 * sin(iangle); ...
                          iCenter(1) - StimHeightPix / 2 * sin(iangle) + StimWidthPix / 2 * cos(iangle), iCenter(2) + StimHeightPix / 2 * cos(iangle) + StimWidthPix / 2 * sin(iangle)];
            Screen('FillPoly', window, white, drawPoints);
        end
        Screen('DrawingFinished', window);
        vbl = Screen('Flip', window);
        if frame == 1
            stimStart = vbl - exptStart;
        end
        [keyIsDown, secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            quit = 1;
        end
    end  % End of stimulus loop

    % Debugging: Stimulus phase completion
    fprintf('PracTrial %d: Stimulus shown with orientations [%.2f, %.2f, %.2f]\n', prac_trial, all_Orientation(prac_trial, 1), all_Orientation(prac_trial, 2), all_Orientation(prac_trial, 3));


    %% Precue
    frame = 0;
    while (quit == 0) && (frame < fprecue) % Start of precue loop
        frame = frame + 1;
        Screen('FillArc', window, gray, [xCenter - yCenter, 0, xCenter + yCenter, 2 * yCenter], 0, 360);
        drawFix(black, fixSize, windowRect, window);
        Screen('DrawingFinished', window);
        vbl = Screen('Flip', window);
        [keyIsDown, secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            quit = 1;
        end
    end  % End of precue loop

    % Debugging: Precue phase completion
    fprintf('PracTrial %d: Precue phase completed\n', prac_trial);

    %% Cue section
    HideCursor();

    % Extract the first_cue and cue colors and angles from the table
    first_cue_colors = [shuffled_table_prac.first_cue_position1(prac_trial), shuffled_table_prac.first_cue_position2(prac_trial), shuffled_table_prac.first_cue_position3(prac_trial)];
    second_cue_colors = [shuffled_table_prac.second_cue_position1(prac_trial), shuffled_table_prac.second_cue_position2(prac_trial), shuffled_table_prac.second_cue_position3(prac_trial)];

    % Extract angles directly from the table and convert to numeric
    cue_angles = [shuffled_table_prac.cue_angle_pos1(prac_trial), ...
                  shuffled_table_prac.cue_angle_pos2(prac_trial), ...
                  shuffled_table_prac.cue_angle_pos3(prac_trial)];
    % cue_angles = [str2double(shuffled_table_prac.cue_angle_pos1{prac_trial}), ...
    %               str2double(shuffled_table_prac.cue_angle_pos2{prac_trial}), ...
    %               str2double(shuffled_table_prac.cue_angle_pos3{prac_trial})];
    %% First Cue Screen with Frame-Based Timing
    frame = 0;
    while (quit == 0) && (frame < fcue1) % Start of first cue screen
        frame = frame + 1;
        Screen('FillArc', window, gray, [xCenter - yCenter, 0, xCenter + yCenter, 2 * yCenter], 0, 360);

        for ii = 1:3
            cueLocCurrent = cueLoc(ii, :); % Cue location
            if first_cue_colors(ii) == 1  % Blue
                Screen('DrawLine', window, blue, cueLocCurrent(1), cueLocCurrent(2), ...
                    cueLocCurrent(1) + cueSize * cosd(cue_angles(ii)), ...
                    cueLocCurrent(2) - cueSize * sind(cue_angles(ii)), cueWidth);  % Blue line
            elseif first_cue_colors(ii) == 2  % Red
                Screen('DrawLine', window, red, cueLocCurrent(1), cueLocCurrent(2), ...
                    cueLocCurrent(1) + cueSize * cosd(cue_angles(ii)), ...
                    cueLocCurrent(2) - cueSize * sind(cue_angles(ii)), cueWidth);  % Red line
            elseif first_cue_colors(ii) == 3  % Green
                Screen('DrawLine', window, green, cueLocCurrent(1), cueLocCurrent(2), ...
                    cueLocCurrent(1) + cueSize * cosd(cue_angles(ii)), ...
                    cueLocCurrent(2) - cueSize * sind(cue_angles(ii)), cueWidth);  % Green line
            end
        end

        Screen('DrawingFinished', window);
        vbl = Screen('Flip', window);  % Flip to show the cues

        [keyIsDown, secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            quit = 1;
        end
    end  % End of first cue screen

    % Debugging First Cue
    fprintf('PracTrial %d: First Cue Colors = [%d, %d, %d]\n', prac_trial, first_cue_colors(1), first_cue_colors(2), first_cue_colors(3));
    fprintf('PracTrial %d: Cue Angles (First Cue) = [%.2f, %.2f, %.2f]\n', prac_trial, cue_angles(1), cue_angles(2), cue_angles(3));

    %% Precue
    frame = 0;
    while (quit == 0) && (frame < fprecue) % Start of precue loop
        frame = frame + 1;
        Screen('FillArc', window, gray, [xCenter - yCenter, 0, xCenter + yCenter, 2 * yCenter], 0, 360);
        drawFix(black, fixSize, windowRect, window);
        Screen('DrawingFinished', window);
        vbl = Screen('Flip', window);
        [keyIsDown, secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            quit = 1;
        end
    end  % End of precue loop

    % Debugging: Precue phase completion
    fprintf('PracTrial %d: Precue phase completed\n', prac_trial);
    %% Second Cue Screen with Frame-Based Timing
    frame = 0;
    while (quit == 0) && (frame < fcue2) % Start of second cue screen
        frame = frame + 1;
        Screen('FillArc', window, gray, [xCenter - yCenter, 0, xCenter + yCenter, 2 * yCenter], 0, 360);

        for ii = 1:3
            cueLocCurrent = cueLoc(ii, :); % Cue location
            if second_cue_colors(ii) == 1  % Blue
                Screen('DrawLine', window, blue, cueLocCurrent(1), cueLocCurrent(2), ...
                    cueLocCurrent(1) + cueSize * cosd(cue_angles(ii)), ...
                    cueLocCurrent(2) - cueSize * sind(cue_angles(ii)), cueWidth);  % Blue line
            elseif second_cue_colors(ii) == 2  % Red
                Screen('DrawLine', window, red, cueLocCurrent(1), cueLocCurrent(2), ...
                    cueLocCurrent(1) + cueSize * cosd(cue_angles(ii)), ...
                    cueLocCurrent(2) - cueSize * sind(cue_angles(ii)), cueWidth);  % Red line
            elseif second_cue_colors(ii) == 3  % Green
                Screen('DrawLine', window, green, cueLocCurrent(1), cueLocCurrent(2), ...
                    cueLocCurrent(1) + cueSize * cosd(cue_angles(ii)), ...
                    cueLocCurrent(2) - cueSize * sind(cue_angles(ii)), cueWidth);  % Green line
            end
        end

        Screen('DrawingFinished', window);
        vbl = Screen('Flip', window);  % Flip to show the second cue

        [keyIsDown, secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            quit = 1;
        end
    end  % End of second cue screen
    
    % Debugging Second Cue
    fprintf('PracTrial %d: Second Cue Colors = [%d, %d, %d], Probe Location = %d\n', prac_trial, second_cue_colors(1), second_cue_colors(2), second_cue_colors(3), prac_probe_loc(prac_trial));
    fprintf('PracTrial %d: Cue Angles (Second Cue) = [%.2f, %.2f, %.2f]\n', prac_trial, cue_angles(1), cue_angles(2), cue_angles(3));
    fprintf('PracTrial %d: Forget condition = %d, Remember condition = %d, Neutral condition = %d\n', prac_trial, prac_test_forget(prac_trial), prac_test_remember(prac_trial), prac_test_neutral(prac_trial));
    %% Delay Period
    frame = 0;
    while (quit == 0) && (frame < fdelay) % Start of delay period
        frame = frame + 1;
        Screen('FillArc', window, gray, [xCenter - yCenter, 0, xCenter + yCenter, 2 * yCenter], 0, 360);
        drawFix(black, fixSize, windowRect, window);
        Screen('DrawingFinished', window);
        vbl = Screen('Flip', window);

        [keyIsDown, secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            quit = 1;
        end
    end  % End of delay period

    %% Response section
    center = stimCenter(prac_probe_loc(prac_trial), :);
    target(prac_trial) = all_Orientation(prac_trial, prac_probe_loc(prac_trial));

    if prac_test_forget(prac_trial) == 1  % Probe at the forget location (Red)
        ShowCursor('Arrow');
        frame = 0; cursorReset = 0;
        fprintf('Entering "no response needed" block for practice trial %d\n', prac_trial);

        while (quit == 0) && (cursorReset == 0) && (debug == 0)
            frame = frame + 1;
            Screen('FillArc', window, gray, [xCenter - yCenter, 0, xCenter + yCenter, 2 * yCenter], 0, 360);
            drawFix(black, fixSize, windowRect, window);
            DrawFormattedText(window, 'This trial does not need response.', 'center', 210, black);
            DrawFormattedText(window, 'Click the center to start next trial', 'center', 250, black);
            Screen('DrawingFinished', window);
            vbl = Screen('Flip', window);

            [mouseX, mouseY, mouseB] = GetMouse(window);
            dist = (((windowRect(3) / 2 - mouseX) ^ 2) + ((windowRect(4) / 2 - mouseY) ^ 2)) ^ 0.5;
            if dist < fixSize / 2 && mouseB(1) == 1
                cursorReset = 1;
            end

            [keyIsDown, secs, keyCode] = KbCheck;
            if keyCode(escapeKey)
                quit = 1;
            end
        end
    else  % Probe at the remember or neutral location
        frame = 0; clicks = 0; time0 = GetSecs;
        while (quit == 0) && (frame < fresp) && (clicks == 0)
            frame = frame + 1;
            ShowCursor;
            [mouseX, mouseY, mouseB] = GetMouse(window);
            deltaX = mouseX - center(1); 
            deltaY = mouseY - center(2);

            % Calculate the angle based on mouse movement
            if deltaX > 0 && deltaY < 0
                angle = -rad2deg(atan(deltaX / deltaY));
            elseif deltaX > 0 && deltaY >= 0
                angle = 90 + rad2deg(atan(deltaY / deltaX));
            elseif deltaX <= 0 && deltaY > 0
                angle = 180 - rad2deg(atan(deltaX / deltaY));
            elseif deltaX < 0 && deltaY < 0
                angle = 270 + rad2deg(atan(deltaY / deltaX));
            end

            if mouseB(1) == 0
                Screen('FillArc', window, gray, [xCenter - yCenter, 0, xCenter + yCenter, 2 * yCenter], 0, 360);
                drawFix(black, fixSize, windowRect, window);
                drawPoints = [center(1) + StimHeightPix / 2 * sin(deg2rad(angle)), center(2) - StimHeightPix / 2 * cos(deg2rad(angle)); ...
                              center(1) - StimHeightPix / 2 * sin(deg2rad(angle)) - StimWidthPix / 2 * cos(deg2rad(angle)), center(2) + StimHeightPix / 2 * cos(deg2rad(angle)) - StimWidthPix / 2 * sin(deg2rad(angle)); ...
                              center(1) - StimHeightPix / 2 * sin(deg2rad(angle)) + StimWidthPix / 2 * cos(deg2rad(angle)), center(2) + StimHeightPix / 2 * cos(deg2rad(angle)) + StimWidthPix / 2 * sin(deg2rad(angle))];
                Screen('FillPoly', window, white, drawPoints);
                Screen('DrawingFinished', window);
                vbl = Screen('Flip', window);
            elseif mouseB(1) == 1
                estimate(prac_trial) = angle;
                rt(prac_trial) = GetSecs - time0;
                error(prac_trial) = rad2deg(circ_dist(deg2rad(estimate(prac_trial)), deg2rad(target(prac_trial))));
                
                % Scoring logic
                if abs(error(prac_trial)) < 3
                    points(prac_trial) = 40;
                elseif abs(error(prac_trial)) >= 3 && abs(error(prac_trial)) < 12
                    points(prac_trial) = 30;
                elseif abs(error(prac_trial)) >= 12 && abs(error(prac_trial)) < 21
                    points(prac_trial) = 20;
                elseif abs(error(prac_trial)) >= 21 && abs(error(prac_trial)) <= 30
                    points(prac_trial) = 10;
                elseif abs(error(prac_trial)) > 30 && abs(error(prac_trial)) <= 40
                    points(prac_trial) = 0;
                elseif abs(error(prac_trial)) > 40 && abs(error(prac_trial)) <= 50
                    points(prac_trial) = -5;
                elseif abs(error(prac_trial)) > 50 && abs(error(prac_trial)) <= 60
                    points(prac_trial) = -10;
                elseif abs(error(prac_trial)) > 60
                    points(prac_trial) = -20;
                end
                clicks = 1;
            end

            [keyIsDown, secs, keyCode] = KbCheck;
            if keyCode(escapeKey)
                quit = 1;
            end
        end
    end  % End of response block

    %% Feedback 
    if prac_test_forget(prac_trial) == 0  % If it's not a forget trial
        frame = 0; cursorReset = 0;
        while (cursorReset == 0) && (quit == 0) && (debug == 0)
            frame = frame + 1;
            Screen('FillArc', window, gray, [xCenter - yCenter, 0, xCenter + yCenter, 2 * yCenter], 0, 360);
            drawFix(black, fixSize, windowRect, window);
            if isnan(estimate(prac_trial))
                DrawFormattedText(window, sprintf('Total Points =  %.0f', sum(points(1:prac_trial))), 'center', 250, black);
                DrawFormattedText(window, 'You missed this trial', 'center', 210, black);
            else
                DrawFormattedText(window, sprintf('Total Points =  %.0f', sum(points(1:prac_trial))), 'center', 250, black);
                DrawFormattedText(window, sprintf('You just gained %.0f points!', points(prac_trial)), 'center', 210, black);
                drawPoints = [center(1)+StimHeightPix/2*sin(deg2rad(angle)),center(2)-StimHeightPix/2*cos(deg2rad(angle)); ...
                        center(1)-StimHeightPix/2*sin(deg2rad(angle))-StimWidthPix/2*cos(deg2rad(angle)),center(2)+StimHeightPix/2*cos(deg2rad(angle))-StimWidthPix/2*sin(deg2rad(angle)); ...
                        center(1)-StimHeightPix/2*sin(deg2rad(angle))+StimWidthPix/2*cos(deg2rad(angle)),center(2)+StimHeightPix/2*cos(deg2rad(angle))+StimWidthPix/2*sin(deg2rad(angle))];
                Screen('FillPoly', window, white, drawPoints);
                drawPoints = [center(1)+StimHeightPix/2*sin(deg2rad(Target(prac_trial))),center(2)-StimHeightPix/2*cos(deg2rad(Target(prac_trial))); ...
                        center(1)-StimHeightPix/2*sin(deg2rad(Target(prac_trial)))-StimWidthPix/2*cos(deg2rad(Target(prac_trial))),center(2)+StimHeightPix/2*cos(deg2rad(Target(prac_trial)))-StimWidthPix/2*sin(deg2rad(Target(prac_trial))); ...
                        center(1)-StimHeightPix/2*sin(deg2rad(Target(prac_trial)))+StimWidthPix/2*cos(deg2rad(Target(prac_trial))),center(2)+StimHeightPix/2*cos(deg2rad(Target(prac_trial)))+StimWidthPix/2*sin(deg2rad(Target(prac_trial)))];
                Screen('FillPoly', window, [0,128,0], drawPoints);
            end
            Screen('DrawingFinished', window);
            vbl = Screen('Flip', window);

            [mouseX, mouseY, mouseB] = GetMouse(window);
            dist = (((windowRect(3) / 2 - mouseX) ^ 2) + ((windowRect(4) / 2 - mouseY) ^ 2)) ^ 0.5;
            if dist < fixSize / 2
                cursorReset = 1;
            end

            [keyIsDown, secs, keyCode] = KbCheck;
            if keyCode(escapeKey)
                quit = 1;
            end
        end % End of feedback block
    end  % End of practice trial loop for Cue
end  % End of entire practice loop
fprintf('All practice trials completed.\n');



%% Start actual trials
practice = 0;

% Start the timer for actual trials
tic;

% Call the function to generate actual trials
[probe_loc, test_remember, test_forget, test_neutral, cue_type] = genActualTrials(nTrials, nBlock, shuffled_table_actual, cue_validity);


block_no = 0; % Reset block number for experimental trials


for trial = 1:totalTrials

 % Determine current block and trial within the block
    currentBlock = ceil(trial / nTrials);
    trialInBlock = trial - (currentBlock - 1) * nTrials;
    
    % Define halfway point
    halfwayPoint = floor(nTrials / 2);

    % Display halfway break message
    if trialInBlock == halfwayPoint
        begin = 0;
        ShowCursor;
        while begin == 0 && quit == 0
            DrawFormattedText(window, sprintf('You are halfway through Block %.0f', currentBlock), 'center', 100, white);
            DrawFormattedText(window, 'Please take a short break', 'center', 200, white);
            DrawFormattedText(window, 'Click the center to continue', 'center', 300, white);
            Screen('FillOval', window, white, makeBox(windowRect(3) / 2, windowRect(4) / 2, fixSize));
            Screen('DrawingFinished', window);
            Screen('Flip', window);

            while any(mouseB)
                [mouseX, mouseY, mouseB] = GetMouse(window);
            end
            while ~any(mouseB)
                [mouseX, mouseY, mouseB] = GetMouse(window);
                centerDist = sqrt((mouseX - windowRect(3) / 2)^2 + (mouseY - windowRect(4) / 2)^2);
                if centerDist < fixSize / 2 && mouseB(1) == 1
                    begin = 1;
                end
                [keyIsDown, ~, keyCode] = KbCheck;
                if keyCode(escapeKey)
                    quit = 1;
                end
            end
            while any(mouseB)
                [~, ~, mouseB] = GetMouse(window);
            end
        end
    end

    % Extract the second cue colors for this trial
    second_cue_colors = [shuffled_table_actual.second_cue_position1(trial), ...
                         shuffled_table_actual.second_cue_position2(trial), ...
                         shuffled_table_actual.second_cue_position3(trial)];

    % Detect the version of the second cue (Type 1 or Type 2)
    if isequal(sort(second_cue_colors), [1 1 3])  % blue, blue, green (Type 1)
        cue_type(trial) = 1;  % Type 1
    elseif isequal(sort(second_cue_colors), [1 2 3])  % green, blue, red (Type 2)
        cue_type(trial) = 2;  % Type 2
    end

    % Debugging: Display which cue type was displayed
    fprintf('Trial %d: Cue Type = %d\n', trial, cue_type(trial));

    % Initialize neutral cues array
    neutralCues = [];  % Temporary array to store multiple neutral cues, if any

    % Determine the cue locations
    remCueLoc(trial) = NaN;       % Initialize to NaN
    forgetCueLoc(trial) = NaN;    % Initialize to NaN

    for ii = 1:3
        if second_cue_colors(ii) == 3  % Green (Remember cue)
            remCueLoc(trial) = ii;     % Set remember cue location
        elseif second_cue_colors(ii) == 2  % Red (Forget cue)
            forgetCueLoc(trial) = ii;  % Set forget cue location
        elseif second_cue_colors(ii) == 1  % Blue (Neutral cue)
            neutralCues = [neutralCues ii];  % Store neutral cue locations
        end
    end
    

    % Assign neutral cue locations to neutralCueLoc (cell array)
    neutralCueLoc{trial} = neutralCues;

    % Extract neutral cue locations (if two exist, otherwise store NaN for the second)
    if length(neutralCues) == 2
        neutral1 = neutralCues(1);
        neutral2 = neutralCues(2);
    elseif length(neutralCues) == 1
        neutral1 = neutralCues(1);
        neutral2 = NaN;  % If there is only one neutral cue, the second one is NaN
    else
        neutral1 = NaN;
        neutral2 = NaN;  % If no neutral cues exist, both are NaN
    end
   
    % Assign conditions with mutual exclusivity
    if probe_loc(trial) == forgetCueLoc(trial)
        test_forget(trial) = 1;   % Set Forget condition
        test_remember(trial) = 0; % Ensure Remember is 0
        test_neutral(trial) = 0;  % Ensure Neutral is 0
    elseif probe_loc(trial) == remCueLoc(trial)
        test_remember(trial) = 1;  % Set Remember condition
        test_forget(trial) = 0;    % Ensure Forget is 0
        test_neutral(trial) = 0;   % Ensure Neutral is 0
    elseif ismember(probe_loc(trial), neutralCues)
        test_neutral(trial) = 1;   % Set Neutral condition
        test_forget(trial) = 0;    % Ensure Forget is 0
        test_remember(trial) = 0;  % Ensure Remember is 0
    end


    if mod(trial, nTrials) == 1
        block_no = block_no + 1;
        begin = 0;
        mouseB = [0, 0]; % Initialize mouseB

        % Ensure mouse button is not pressed before displaying the message
        while any(mouseB)
            [~, ~, mouseB] = GetMouse(window);
        end

         % Adding a small pause to ensure the state transition
           WaitSecs(1);



        while begin == 0 && quit == 0
            DrawFormattedText(window, sprintf('Block %.0f', block_no), 'center', 100, white);
            DrawFormattedText(window, 'GREEN ITEM IS TWICE LIKELY TO BE TESTED,\nRED ITEM DOES NOT NEED RESPONSE AND CAN BE FORGOTTEN', 'center', 250, red);
            DrawFormattedText(window, sprintf('Click on the middle to start block %.0f', block_no), 'center', 350, white);     Screen('FillOval', window, white, makeBox(windowRect(3) / 2, windowRect(4) / 2, fixSize));
            Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
            Screen('Flip', window);

            while any(mouseB) % if already down, wait for release
                [mouseX, mouseY, mouseB] = GetMouse(window);
            end
            while ~any(mouseB) % wait for press
                [mouseX, mouseY, mouseB] = GetMouse(window);
                centerDist = ((mouseX - windowRect(3) / 2) ^ 2 + (mouseY - windowRect(4) / 2) ^ 2) ^ 0.5;
                if centerDist < fixSize / 2
                    if mouseB(1) == 1
                        begin = 1;
                    end
                end
                [keyIsDown, secs, keyCode] = KbCheck;
                            
               if keyCode(escapeKey)
                    quit = 1;
                    if bevID > 0
                        fclose(bevID);
                    end
                    cd(mainPath);
                    sca; Priority(0);
              end
           end

            % Ensure mouse button is released after clicking in the center
            while any(mouseB)
                [~, ~, mouseB] = GetMouse(window);
            end
        end
    end
       
    

    %% Fixation
    frame = 0; HideCursor();
    while (quit == 0) && (frame < ffix)
        frame = frame + 1;
        Screen('FillArc', window, gray, [xCenter - yCenter, 0, xCenter + yCenter, 2 * yCenter], 0, 360);
        drawFix(black, fixSize, windowRect, window);
        Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
        vbl = Screen('Flip', window);
        [keyIsDown, secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            quit = 1;
        end
    end

    %% Stimulus
    index_Gabor = randperm(length(orientation));
    All_Orientation(trial, :) = [orientation(index_Gabor(1)) + jitter * rand(), orientation(index_Gabor(2)) + jitter * rand(), orientation(index_Gabor(3)) + jitter * rand()];
    frame = 0;
    while (quit == 0) && (frame < fstim)
        frame = frame + 1;
        Screen('FillArc', window, gray, [xCenter - yCenter, 0, xCenter + yCenter, 2 * yCenter], 0, 360);
        drawFix(black, fixSize, windowRect, window);
        for ii = 1:3
            iCenter = stimCenter(ii, :);
            iangle = deg2rad(All_Orientation(trial, ii));
            drawPoints = [iCenter(1) + StimHeightPix / 2 * sin(iangle), iCenter(2) - StimHeightPix / 2 * cos(iangle); ...
                iCenter(1) - StimHeightPix / 2 * sin(iangle) - StimWidthPix / 2 * cos(iangle), iCenter(2) + StimHeightPix / 2 * cos(iangle) - StimWidthPix / 2 * sin(iangle); ...
                iCenter(1) - StimHeightPix / 2 * sin(iangle) + StimWidthPix / 2 * cos(iangle), iCenter(2) + StimHeightPix / 2 * cos(iangle) + StimWidthPix / 2 * sin(iangle)];
            Screen('FillPoly', window, white, drawPoints);
        end
        Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
        vbl = Screen('Flip', window);
        if frame == 1
            stimStart = vbl - exptStart;
        end
        [keyIsDown, secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            quit = 1;
        end
    end

    %% Precue
    frame = 0;
    while (quit == 0) && (frame < fprecue)
        frame = frame + 1;
        Screen('FillArc', window, gray, [xCenter - yCenter, 0, xCenter + yCenter, 2 * yCenter], 0, 360);
        drawFix(black, fixSize, windowRect, window);
        Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
        vbl = Screen('Flip', window);
        [keyIsDown, secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            quit = 1;
        end
    end

    %% Cue section
    HideCursor();

    % Extract the first_cue and second_cue colors and angles from the table
    first_cue_colors = [shuffled_table_actual.first_cue_position1(trial), shuffled_table_actual.first_cue_position2(trial), shuffled_table_actual.first_cue_position3(trial)];
    second_cue_colors = [shuffled_table_actual.second_cue_position1(trial), shuffled_table_actual.second_cue_position2(trial), shuffled_table_actual.second_cue_position3(trial)];

    % Extract angles directly from the table and convert to numeric
    cue_angles = [shuffled_table_actual.cue_angle_pos1(prac_trial), ...
                  shuffled_table_actual.cue_angle_pos2(prac_trial), ...
                  shuffled_table_actual.cue_angle_pos3(prac_trial)];

    %% First Cue Screen with Frame-Based Timing
    frame = 0;
    while (quit == 0) && (frame < fcue1) % Start of first cue screen
        frame = frame + 1;
        Screen('FillArc', window, gray, [xCenter - yCenter, 0, xCenter + yCenter, 2 * yCenter], 0, 360);

        for ii = 1:3
            cueLocCurrent = cueLoc(ii, :); % Cue location
            if first_cue_colors(ii) == 1  % Blue
                Screen('DrawLine', window, blue, cueLocCurrent(1), cueLocCurrent(2), ...
                    cueLocCurrent(1) + cueSize * cosd(cue_angles(ii)), ...
                    cueLocCurrent(2) - cueSize * sind(cue_angles(ii)), cueWidth);  % Blue line
            elseif first_cue_colors(ii) == 2  % Red
                Screen('DrawLine', window, red, cueLocCurrent(1), cueLocCurrent(2), ...
                    cueLocCurrent(1) + cueSize * cosd(cue_angles(ii)), ...
                    cueLocCurrent(2) - cueSize * sind(cue_angles(ii)), cueWidth);  % Red line
            elseif first_cue_colors(ii) == 3  % Green
                Screen('DrawLine', window, green, cueLocCurrent(1), cueLocCurrent(2), ...
                    cueLocCurrent(1) + cueSize * cosd(cue_angles(ii)), ...
                    cueLocCurrent(2) - cueSize * sind(cue_angles(ii)), cueWidth);  % Green line
            end
        end

        Screen('DrawingFinished', window);
        vbl = Screen('Flip', window);  % Flip to show the cues

        [keyIsDown, secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            quit = 1;
        end
    end  % End of first cue screen

    % Debugging First Cue
    fprintf('Trial %d: First Cue Colors = [%d, %d, %d]\n', trial, first_cue_colors(1), first_cue_colors(2), first_cue_colors(3));
    fprintf('Trial %d: Cue Angles (First Cue) = [%.2f, %.2f, %.2f]\n', trial, cue_angles(1), cue_angles(2), cue_angles(3));

    %% Second Cue Screen with Frame-Based Timing
    frame = 0;
    while (quit == 0) && (frame < fcue2) % Start of second cue screen
        frame = frame + 1;
        Screen('FillArc', window, gray, [xCenter - yCenter, 0, xCenter + yCenter, 2 * yCenter], 0, 360);

        for ii = 1:3
            cueLocCurrent = cueLoc(ii, :); % Cue location
            if second_cue_colors(ii) == 1  % Blue
                Screen('DrawLine', window, blue, cueLocCurrent(1), cueLocCurrent(2), ...
                    cueLocCurrent(1) + cueSize * cosd(cue_angles(ii)), ...
                    cueLocCurrent(2) - cueSize * sind(cue_angles(ii)), cueWidth);  % Blue line
            elseif second_cue_colors(ii) == 2  % Red
                Screen('DrawLine', window, red, cueLocCurrent(1), cueLocCurrent(2), ...
                    cueLocCurrent(1) + cueSize * cosd(cue_angles(ii)), ...
                    cueLocCurrent(2) - cueSize * sind(cue_angles(ii)), cueWidth);  % Red line
            elseif second_cue_colors(ii) == 3  % Green
                Screen('DrawLine', window, green, cueLocCurrent(1), cueLocCurrent(2), ...
                    cueLocCurrent(1) + cueSize * cosd(cue_angles(ii)), ...
                    cueLocCurrent(2) - cueSize * sind(cue_angles(ii)), cueWidth);  % Green line
            end
        end

        Screen('DrawingFinished', window);
        vbl = Screen('Flip', window);  % Flip to show the second cue

        [keyIsDown, secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            quit = 1;
        end
    end  % End of second cue screen

    % Debugging Second Cue
    fprintf('Trial %d: Second Cue Colors = [%d, %d, %d], Probe Location = %d\n', trial, second_cue_colors(1), second_cue_colors(2), second_cue_colors(3), probe_loc(trial));
    fprintf('Trial %d: Cue Angles (Second Cue) = [%.2f, %.2f, %.2f]\n', trial, cue_angles(1), cue_angles(2), cue_angles(3));
    fprintf('Trial %d: Forget condition = %d, Remember condition = %d, Neutral condition = %d\n', trial, test_forget(trial), test_remember(trial), test_neutral(trial));

    %% Delay Period
    frame = 0;
    while (quit == 0) && (frame < fdelay) % Start of delay period
        frame = frame + 1;
        Screen('FillArc', window, gray, [xCenter - yCenter, 0, xCenter + yCenter, 2 * yCenter], 0, 360);
        drawFix(black, fixSize, windowRect, window);
        Screen('DrawingFinished', window);
        vbl = Screen('Flip', window);

        [keyIsDown, secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            quit = 1;
        end
    end  % End of delay period

    %% Response section
    center = stimCenter(probe_loc(trial), :);
    Target(trial) = All_Orientation(trial, probe_loc(trial));

    if test_forget(trial) == 1  % Probe at the forget location (Red)
        ShowCursor('Arrow');
        frame = 0; cursorReset = 0;
        fprintf('Entering "no response needed" block for actual trial %d\n', trial);

        while (quit == 0) && (cursorReset == 0) && (debug == 0)
            frame = frame + 1;
            Screen('FillArc', window, gray, [xCenter - yCenter, 0, xCenter + yCenter, 2 * yCenter], 0, 360);
            drawFix(black, fixSize, windowRect, window);
            DrawFormattedText(window, 'This trial does not need a response.', 'center', 210, black);
            DrawFormattedText(window, 'Click the center to start the next trial', 'center', 250, black);
            Screen('DrawingFinished', window);
            vbl = Screen('Flip', window);

            [mouseX, mouseY, mouseB] = GetMouse(window);
            dist = (((windowRect(3) / 2 - mouseX) ^ 2) + ((windowRect(4) / 2 - mouseY) ^ 2)) ^ 0.5;
            if dist < fixSize / 2 && mouseB(1) == 1
                cursorReset = 1;
            end

            [keyIsDown, secs, keyCode] = KbCheck;
            if keyCode(escapeKey)
                quit = 1;
            end
        end
    else  % Probe at the remember or neutral location
        frame = 0; clicks = 0; time0 = GetSecs;
        while (quit == 0) && (frame < fresp) && (clicks == 0)
            frame = frame + 1;
            ShowCursor;
            [mouseX, mouseY, mouseB] = GetMouse(window);
            deltaX = mouseX - center(1); 
            deltaY = mouseY - center(2);

            % Calculate the angle based on mouse movement
            if deltaX > 0 && deltaY < 0
                angle = -rad2deg(atan(deltaX / deltaY));
            elseif deltaX > 0 && deltaY >= 0
                angle = 90 + rad2deg(atan(deltaY / deltaX));
            elseif deltaX <= 0 && deltaY > 0
                angle = 180 - rad2deg(atan(deltaX / deltaY));
            elseif deltaX < 0 && deltaY < 0
                angle = 270 + rad2deg(atan(deltaY / deltaX));
            end

            if mouseB(1) == 0
                Screen('FillArc', window, gray, [xCenter - yCenter, 0, xCenter + yCenter, 2 * yCenter], 0, 360);
                drawFix(black, fixSize, windowRect, window);
                drawPoints = [center(1) + StimHeightPix / 2 * sin(deg2rad(angle)), center(2) - StimHeightPix / 2 * cos(deg2rad(angle)); ...
                              center(1) - StimHeightPix / 2 * sin(deg2rad(angle)) - StimWidthPix / 2 * cos(deg2rad(angle)), center(2) + StimHeightPix / 2 * cos(deg2rad(angle)) - StimWidthPix / 2 * sin(deg2rad(angle)); ...
                              center(1) - StimHeightPix / 2 * sin(deg2rad(angle)) + StimWidthPix / 2 * cos(deg2rad(angle)), center(2) + StimHeightPix / 2 * cos(deg2rad(angle)) + StimWidthPix / 2 * sin(deg2rad(angle))];
                Screen('FillPoly', window, white, drawPoints);
                Screen('DrawingFinished', window);
                vbl = Screen('Flip', window);
            elseif mouseB(1) == 1
                Estimate(trial) = angle;
                RT(trial) = GetSecs - time0;
                Error(trial) = rad2deg(circ_dist(deg2rad(Estimate(trial)), deg2rad(Target(trial))));

                % Scoring logic
                if abs(Error(trial)) < 3
                    Points(trial) = 40;
                elseif abs(Error(trial)) >= 3 && abs(Error(trial)) < 12
                    Points(trial) = 30;
                elseif abs(Error(trial)) >= 12 && abs(Error(trial)) < 21
                    Points(trial) = 20;
                elseif abs(Error(trial)) >= 21 && abs(Error(trial)) <= 30
                    Points(trial) = 10;
                elseif abs(Error(trial)) > 30 && abs(Error(trial)) <= 40
                    Points(trial) = 0;
                elseif abs(Error(trial)) > 40 && abs(Error(trial)) <= 50
                    Points(trial) = -5;
                elseif abs(Error(trial)) > 50 && abs(Error(trial)) <= 60
                    Points(trial) = -10;
                elseif abs(Error(trial)) > 60
                    Points(trial) = -20;
                end
                clicks = 1;
            end

            [keyIsDown, secs, keyCode] = KbCheck;
            if keyCode(escapeKey)
                quit = 1;
            end
        end
    end  % End of response block

  %% Save
    
    if quit == 0
    
            % Handle NaN values explicitly before saving
        if isnan(Estimate(trial))
            estimate_string = 'NaN';
        else
            estimate_string = num2str(Estimate(trial));
        end
        
        if isnan(Target(trial))
            target_string = 'NaN';
        else
            target_string = num2str(Target(trial));
        end
        
        if isnan(forgetCueLoc(trial))
            forget_string = 'NaN';
        else
            forget_string = num2str(forgetCueLoc(trial));
        end
        
        if isnan(remCueLoc(trial))
            remember_string = 'NaN';
        else
            remember_string = num2str(remCueLoc(trial));
        end
        
        if isnan(neutral1)
            neutral1_string = 'NaN';
        else
            neutral1_string = num2str(neutral1);
        end
        
        if isnan(neutral2)
            neutral2_string = 'NaN';
        else
            neutral2_string = num2str(neutral2);
        end
    
        fprintf('Trial %d: Target = %.2f, Estimate = %.2f, RT = %.2f, Error = %.2f, Points = %d\n', trial, Target(trial), Estimate(trial), RT(trial), Error(trial), Points(trial));
    
        datastring = strcat(subject, ',', session, ',', num2str(block_no), ',' ,num2str(trial), ',' , num2str(cue_type(trial)), ',' , num2str(All_Orientation(trial,1)) , ',', num2str(All_Orientation(trial,2)), ',', num2str(All_Orientation(trial,3)), ',',  ...
            num2str(remCueLoc(trial)), ',', num2str(forgetCueLoc(trial)), ',', num2str(neutral1) , ',', num2str(neutral2), ',', num2str(probe_loc(trial)), ',', num2str(test_remember(trial)), ',', num2str(test_forget(trial)), ',', num2str(test_neutral(trial)), ',', ...
            num2str(Target(trial)), ',',  num2str(Estimate(trial)), ',',  num2str(RT(trial)), ',', num2str(Error(trial)), ',',  num2str(Points(trial)), ',', '\n');
    
        % Write the datastring to the file
        fprintf(bevID, datastring);
    end




    %% Feedback 

    if test_forget(trial) == 0  % If it's not a forget trial
        frame = 0; cursorReset = 0;
        while (cursorReset == 0) && (quit == 0) && (debug == 0)
            frame = frame + 1;
            Screen('FillArc', window, gray, [xCenter - yCenter, 0, xCenter + yCenter, 2 * yCenter], 0, 360);
            drawFix(black, fixSize, windowRect, window);
            if isnan(Estimate(trial))
                DrawFormattedText(window, sprintf('Total Points =  %.0f', sum(Points(1:trial))), 'center', 250, black);
                DrawFormattedText(window, 'You missed this trial', 'center', 210, black);
                 drawPoints = [center(1)+StimHeightPix/2*sin(deg2rad(angle)),center(2)-StimHeightPix/2*cos(deg2rad(angle)); ...
                            center(1)-StimHeightPix/2*sin(deg2rad(angle))-StimWidthPix/2*cos(deg2rad(angle)),center(2)+StimHeightPix/2*cos(deg2rad(angle))-StimWidthPix/2*sin(deg2rad(angle)); ...
                            center(1)-StimHeightPix/2*sin(deg2rad(angle))+StimWidthPix/2*cos(deg2rad(angle)),center(2)+StimHeightPix/2*cos(deg2rad(angle))+StimWidthPix/2*sin(deg2rad(angle))];
                    Screen('FillPoly', window, white, drawPoints);
                    drawPoints = [center(1)+StimHeightPix/2*sin(deg2rad(Target(trial))),center(2)-StimHeightPix/2*cos(deg2rad(Target(trial))); ...
                            center(1)-StimHeightPix/2*sin(deg2rad(Target(trial)))-StimWidthPix/2*cos(deg2rad(Target(trial))),center(2)+StimHeightPix/2*cos(deg2rad(Target(trial)))-StimWidthPix/2*sin(deg2rad(Target(trial))); ...
                            center(1)-StimHeightPix/2*sin(deg2rad(Target(trial)))+StimWidthPix/2*cos(deg2rad(Target(trial))),center(2)+StimHeightPix/2*cos(deg2rad(Target(trial)))+StimWidthPix/2*sin(deg2rad(Target(trial)))];
                    Screen('FillPoly', window, [0,128,0], drawPoints);
            else
                DrawFormattedText(window, sprintf('Total Points =  %.0f', sum(Points(1:trial))), 'center', 250, black);
                DrawFormattedText(window, sprintf('You just gained %.0f points!', Points(trial)), 'center', 210, black);
                drawPoints = [center(1)+StimHeightPix/2*sin(deg2rad(angle)),center(2)-StimHeightPix/2*cos(deg2rad(angle)); ...
                            center(1)-StimHeightPix/2*sin(deg2rad(angle))-StimWidthPix/2*cos(deg2rad(angle)),center(2)+StimHeightPix/2*cos(deg2rad(angle))-StimWidthPix/2*sin(deg2rad(angle)); ...
                            center(1)-StimHeightPix/2*sin(deg2rad(angle))+StimWidthPix/2*cos(deg2rad(angle)),center(2)+StimHeightPix/2*cos(deg2rad(angle))+StimWidthPix/2*sin(deg2rad(angle))];
                    Screen('FillPoly', window, white, drawPoints);
                    drawPoints = [center(1)+StimHeightPix/2*sin(deg2rad(Target(trial))),center(2)-StimHeightPix/2*cos(deg2rad(Target(trial))); ...
                            center(1)-StimHeightPix/2*sin(deg2rad(Target(trial)))-StimWidthPix/2*cos(deg2rad(Target(trial))),center(2)+StimHeightPix/2*cos(deg2rad(Target(trial)))-StimWidthPix/2*sin(deg2rad(Target(trial))); ...
                            center(1)-StimHeightPix/2*sin(deg2rad(Target(trial)))+StimWidthPix/2*cos(deg2rad(Target(trial))),center(2)+StimHeightPix/2*cos(deg2rad(Target(trial)))+StimWidthPix/2*sin(deg2rad(Target(trial)))];
                    Screen('FillPoly', window, [0,128,0], drawPoints);
            end
            Screen('DrawingFinished', window);
            vbl = Screen('Flip', window);

            [mouseX, mouseY, mouseB] = GetMouse(window);
            dist = (((windowRect(3) / 2 - mouseX) ^ 2) + ((windowRect(4) / 2 - mouseY) ^ 2)) ^ 0.5;
            if dist < fixSize / 2
                cursorReset = 1;
            end

            [keyIsDown, secs, keyCode] = KbCheck;
            if keyCode(escapeKey)
                quit = 1;
            end
        end % End of feedback block
    end  % End of actual trial loop for Cue

    %% End of Block
    if mod(trial, nTrials) == 0 && trial ~= nTrials * nBlock
        begin = 0; ShowCursor;
        while begin == 0 && quit == 0
            DrawFormattedText(window, sprintf('You completed %.0f of %.0f blocks', block_no, nBlock), 'center', 100, white);
            DrawFormattedText(window, 'Please take a short break', 'center', 200, white);
            DrawFormattedText(window, 'Click the center to continue', 'center', 300, white);
            Screen('FillOval', window, white, makeBox(windowRect(3) / 2, windowRect(4) / 2, fixSize));
            Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
            Screen('Flip', window);
            
            while any(mouseB) % if already down, wait for release
                [mouseX, mouseY, mouseB] = GetMouse(window);
            end
            while ~any(mouseB) % wait for press
                [mouseX, mouseY, mouseB] = GetMouse(window);
                centerDist = ((mouseX - windowRect(3) / 2) ^ 2 + (mouseY - windowRect(4) / 2) ^ 2) ^ 0.5;
                if centerDist < fixSize / 2
                    if mouseB(1) == 1
                        begin = 1;
                    end
                end
            end

            [keyIsDown, secs, keyCode] = KbCheck;
            if keyCode(escapeKey)
                quit = 1;
                fclose(bevID);
                cd(mainPath);
                sca; Priority(0);
            end
        end
    elseif trial == nTrials * nBlock
        duration = GetSecs - startTime;
        Screen('Flip', window); % blank
        disp('Experiment successfully completed!');

        % Calculate the bonus
          bonus = round(sum(Points) / (40 * (nTrials * nBlock)) * maxBonus);
        % Ensure bonus is non-negative
        if bonus < 0
           bonus = 0; % Force bonus to be at least 0
        end

        while quit == 0
            DrawFormattedText(window, 'You have completed the experiment.', 'center', 200, white);
            DrawFormattedText(window, 'Please find the experimenter.', 'center', 300, white);
             DrawFormattedText(window, sprintf('Your bonus is %d AED', bonus), 'center', 400, white);
            Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
            Screen('Flip', window);
            
            [keyIsDown, secs, keyCode] = KbCheck;
            if keyCode(escapeKey)
                quit = 1;
            end
        end
    end
end  % End of entire actual trials loop

% Stop the timer after actual trials
elapsedTime = toc;

% Display the elapsed time
fprintf('Total time for actual trials: %.2f seconds\n', elapsedTime);

%% End
fclose(bevID); 
cd(mainPath);
bonus = round(sum(Points) / (40 * (nTrials * nBlock)) * maxBonus);
% Ensure bonus is non-negative
if bonus < 0
    bonus = 0; % Force bonus to be at least 0
end
fprintf('*******************\nBonus = %d AED\n*******************\n', bonus);
sca;
Priority(0);

%% Helper Functions
function drawFix(color, size, windowRect, window)
    Screen('FillOval', window, color, makeBox(windowRect(3) / 2, windowRect(4) / 2, size));
end

function Box = makeBox(Xpos, Ypos, size)
    Box = [Xpos - size / 2, Ypos - size / 2, Xpos + size / 2, Ypos + size / 2];
end

%% Practice Function

function [prac_probe_loc, prac_test_remember, prac_test_forget, prac_test_neutral, prac_cue_type] = genPracticeTrials(numPrac, shuffled_table_prac, cue_validity)

    % Initialize arrays for storing practice trial conditions
    prac_probe_loc = zeros(numPrac, 1);     % Probe location for each trial
    prac_test_remember = zeros(numPrac, 1); % Remember condition (green)
    prac_test_forget = zeros(numPrac, 1);   % Forget condition (red)
    prac_test_neutral = zeros(numPrac, 1);  % Neutral condition (blue)
    prac_cue_type = zeros(numPrac, 1);      % To store the cue type (1 or 2) for each trial
    
    % Define quotas for Cue Type 1 and Cue Type 2 in practice trials
    nCueType1 = round(numPrac / 3);         % 1/3 of practice trials for Cue Type 1
    nCueType2 = numPrac - nCueType1;        % 2/3 of practice trials for Cue Type 2
    
    % Set quotas for Cue Type 1
    quota_type1_remember = round(nCueType1 * cue_validity);
    quota_type1_neutral = nCueType1 - quota_type1_remember;
    
    % Set quotas for Cue Type 2
    quota_type2_remember = round(nCueType2 * cue_validity);
    quota_type2_forget = round(nCueType2 * 1 / 6);
    quota_type2_neutral = nCueType2 - (quota_type2_remember + quota_type2_forget);
    
    % Track counts for each condition
    count_type1_remember = 0;
    count_type1_neutral = 0;
    count_type2_remember = 0;
    count_type2_forget = 0;
    count_type2_neutral = 0;

    % Assign trials based on cue type and conditions
    for i = 1:numPrac
        % Get the cue colors for this trial from the shuffled table
        second_cue_colors = [shuffled_table_prac.second_cue_position1(i), ...
                             shuffled_table_prac.second_cue_position2(i), ...
                             shuffled_table_prac.second_cue_position3(i)];
        
        % Determine Cue Type 1 or Cue Type 2 based on colors and quota
        if prac_cue_type(i) == 0  % Check if cue type is not yet assigned
            if sum(prac_cue_type == 1) < nCueType1
                prac_cue_type(i) = 1;  % Assign Cue Type 1
            else
                prac_cue_type(i) = 2;  % Assign Cue Type 2
            end
        end
        
        % Assign trials according to Cue Type quotas
        if prac_cue_type(i) == 1
            % Cue Type 1 (1/3 of practice trials): Remember or Neutral
            if count_type1_remember < quota_type1_remember
                green_loc = find(second_cue_colors == 3);
                if ~isempty(green_loc)
                    prac_probe_loc(i) = green_loc(randperm(length(green_loc), 1));
                    prac_test_remember(i) = 1;
                    count_type1_remember = count_type1_remember + 1;
                end
            elseif count_type1_neutral < quota_type1_neutral
                blue_locs = find(second_cue_colors == 1);
                if ~isempty(blue_locs)
                    prac_probe_loc(i) = blue_locs(randperm(length(blue_locs), 1));
                    prac_test_neutral(i) = 1;
                    count_type1_neutral = count_type1_neutral + 1;
                end
            end
            
        elseif prac_cue_type(i) == 2
            % Cue Type 2 (2/3 of practice trials): Remember, Forget, or Neutral
            if count_type2_remember < quota_type2_remember
                green_loc = find(second_cue_colors == 3);
                if ~isempty(green_loc)
                    prac_probe_loc(i) = green_loc(randperm(length(green_loc), 1));
                    prac_test_remember(i) = 1;
                    count_type2_remember = count_type2_remember + 1;
                end
            elseif count_type2_forget < quota_type2_forget
                red_loc = find(second_cue_colors == 2);
                if ~isempty(red_loc)
                    prac_probe_loc(i) = red_loc(randperm(length(red_loc), 1));
                    prac_test_forget(i) = 1;
                    count_type2_forget = count_type2_forget + 1;
                end
            elseif count_type2_neutral < quota_type2_neutral
                blue_locs = find(second_cue_colors == 1);
                if ~isempty(blue_locs)
                    prac_probe_loc(i) = blue_locs(randperm(length(blue_locs), 1));
                    prac_test_neutral(i) = 1;
                    count_type2_neutral = count_type2_neutral + 1;
                end
            end
        end
    end
end


%% Actual Trials Function
function [probe_loc, test_remember, test_forget, test_neutral, cue_type] = genActualTrials(nTrials, nBlock, shuffled_table_actual, cue_validity)

    totalTrials = nTrials * nBlock;
    probe_loc = zeros(totalTrials, 1);     
    test_remember = zeros(totalTrials, 1); 
    test_forget = zeros(totalTrials, 1);   
    test_neutral = zeros(totalTrials, 1);  
    cue_type = shuffled_table_actual.CueType; 

    for block = 1:nBlock
        blockStart = (block - 1) * nTrials + 1;
        blockEnd = block * nTrials;

        nCueType1 = nTrials / 3;
        nCueType2 = nTrials - nCueType1;

        quota_type1_remember = nCueType1 * cue_validity;
        quota_type1_neutral = nCueType1 - quota_type1_remember;

        quota_type2_remember = nCueType2 * cue_validity;
        quota_type2_forget = nCueType2 * 1 / 6;
        quota_type2_neutral = nCueType2 - (quota_type2_remember + quota_type2_forget);

        % Shuffle trial indices for each block to prevent clustering
        trial_indices = blockStart:blockEnd;
        trial_indices = trial_indices(randperm(length(trial_indices)));
        
        for trial_index = trial_indices
            cueType = cue_type(trial_index);

            second_cue_colors = [shuffled_table_actual.second_cue_position1(trial_index), ...
                                 shuffled_table_actual.second_cue_position2(trial_index), ...
                                 shuffled_table_actual.second_cue_position3(trial_index)];

            green_locs = find(second_cue_colors == 3);
            red_locs = find(second_cue_colors == 2);
            blue_locs = find(second_cue_colors == 1);

            assigned = false;

            if cueType == 1
                if quota_type1_remember > 0 && ~isempty(green_locs)
                    probe_loc(trial_index) = green_locs(1);
                    test_remember(trial_index) = 1;
                    quota_type1_remember = quota_type1_remember - 1;
                    assigned = true;
                elseif quota_type1_neutral > 0 && ~isempty(blue_locs)
                    probe_loc(trial_index) = blue_locs(1);
                    test_neutral(trial_index) = 1;
                    quota_type1_neutral = quota_type1_neutral - 1;
                    assigned = true;
                end
            elseif cueType == 2
                if quota_type2_remember > 0 && ~isempty(green_locs)
                    probe_loc(trial_index) = green_locs(1);
                    test_remember(trial_index) = 1;
                    quota_type2_remember = quota_type2_remember - 1;
                    assigned = true;
                elseif quota_type2_forget > 0 && ~isempty(red_locs)
                    probe_loc(trial_index) = red_locs(1);
                    test_forget(trial_index) = 1;
                    quota_type2_forget = quota_type2_forget - 1;
                    assigned = true;
                elseif quota_type2_neutral > 0 && ~isempty(blue_locs)
                    probe_loc(trial_index) = blue_locs(1);
                    test_neutral(trial_index) = 1;
                    quota_type2_neutral = quota_type2_neutral - 1;
                    assigned = true;
                end
            end

            if ~assigned
                error('Insufficient trials to meet quotas in trial %d of block %d.', trial_index, block);
            end
        end

        % After assigning for this block, check if quotas are met
        if quota_type1_remember ~= 0 || quota_type1_neutral ~= 0 || quota_type2_remember ~= 0 || quota_type2_forget ~= 0 || quota_type2_neutral ~= 0
            error('Not all quotas were met after assignment in block %d.', block);
        end
    end

    % Total counts for verification
    total_type1_remember = sum(test_remember & cue_type == 1);
    total_type1_neutral = sum(test_neutral & cue_type == 1);
    total_type2_remember = sum(test_remember & cue_type == 2);
    total_type2_forget = sum(test_forget & cue_type == 2);
    total_type2_neutral = sum(test_neutral & cue_type == 2);

    fprintf('\nFinal Overall Counts:\n');
    fprintf('Cue Type 1 - Remember: %d, Neutral: %d\n', total_type1_remember, total_type1_neutral);
    fprintf('Cue Type 2 - Remember: %d, Forget: %d, Neutral: %d\n', total_type2_remember, total_type2_forget, total_type2_neutral);

end
