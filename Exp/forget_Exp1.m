sca; close all; clear all;

scriptpath = mfilename('fullpath');
scriptname = mfilename();

mainPath = scriptpath(1:length(scriptpath)-length(scriptname));
addpath(genpath([mainPath 'CircStat2012a']));
dataPath = strcat(mainPath,'Data/');

if exist(dataPath,'dir')==0
    mkdir Data
end

expt = 'forget_beh';
debug = 0;

sel_font = 'Courier New';

maxBonus = 25;
%% Logfile set-up
cd(dataPath);
input = 0;
while input  == 0
    prompt = {'Subject Code:','Session:','Handedness:','Gender:','Age:'};
    titlen = 'Welcome! Please get the experimenter to input the following:';
    dims = [1 100];
    definput = {'101','1','right', 'female', '20'};
    answer = inputdlg(prompt,titlen,dims,definput);
    
    try
        subject = answer{1};
        session = answer{2};
        hand = answer{3};
        gender = answer{4};
        age = answer{5};
        input = 1;
    catch
        disp('Error with input. Try again.')
    end
end

if mod(str2double(subject),2) == 1 % odd subID - ABBA; even subID - BAAB
    order = [1 1 2 2 2 2 1 1]; %1 : remember; 2: forget
else
    order = [2 2 1 1 1 1 2 2];
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
headerString = 'Subject,Session,Block,Condition,Trial,whetherCued,leftStim,rightStim,midStim,remCueLoc,fogCueLoc,neutralCueLoc,probeCueLoc,probeRem,probeForget,probeNeutral,Target,Estimate,RT,Error,Points\n';
fprintf(bevID,headerString);
jitter = 5;
%% Expt Parameters
if debug == 0
    screenNumber = 0;
    nTrials = 60; % trials per block
    nBlock = 8;
    numPrac = 8;
    tfix = 500;
    tstim = 150; %ms
    tmask = 500;
    tprecue = 500;
    tcue = 500;
    tdelay = 1500;
    tresp = 100000;
else
    screenNumber = 1;
    nTrials = 60; % trials per block
    nBlock = 8;
    numPrac = 8;
    tfix = 2;
    tstim = 2; %ms
    tmask = 2;
    tprecue = 2;
    tcue = 2;
    tdelay = 2;
    tresp = 2;
end


white = 1.0; black = 0.0;
gray = [0.5 0.5 0.5]; yellow = [255,165,0]./255;
red = [1,0,0];
fixSize = 12;
orientation = 15:30:360;

cue_percent = 3/5;
cue_validity = 2/3;

% record
All_Orientation = zeros(nTrials*nBlock,3);
remember_loc = zeros(nTrials*nBlock,1);
remember_Orientation = zeros(nTrials*nBlock,1);
forget_loc = zeros(nTrials*nBlock,1);
forget_Orientation = zeros(nTrials*nBlock,1);
neutral_loc = zeros(nTrials*nBlock,1);
neutral_Orientation = zeros(nTrials*nBlock,1);
probe_Orientation = zeros(nTrials*nBlock,1);

RT = ones(nTrials*nBlock,1)*tresp;
Estimate = NaN(nTrials*nBlock,1);
Error = 999*ones(nTrials*nBlock,1);
Points = zeros(nTrials*nBlock,1);
Target = 999*ones(nTrials*nBlock,1);

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
viewDist = 88; %the distance from the screen, in cm
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

% prepare frame numbers for every even
frameRate = round(FrameRate(window));
avgFrame = 1000/frameRate;
ffix = round(tfix/avgFrame);
fstim = round(tstim/avgFrame);
fmask = round(tmask/avgFrame);
fprecue = round(tprecue/avgFrame);
fcue = round(tcue/avgFrame);
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

% prepar cue
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
exptStart = Screen('Flip', window);  % let's mark this as the start of the expt
Screen('Preference', 'TextRenderer', 1);
Screen('TextSize', window, 25);
Screen('TextFont', window, sel_font,1); % 1 = bold
escapeKey = KbName('ESCAPE');
Priority(topPriorityLevel);

% get task map
[cueCol_rem, probe_loc_rem, test_remember_rem, test_forget_rem, test_neutral_rem, cue_trials_rem, cue_loc_index_all_rem] = genTrials(cue_percent, cue_validity, nTrials*nBlock/2, 1, 0); % remember trials
[cueCol_forg, probe_loc_forg, test_remember_forg, test_forget_forg, test_neutral_forg, cue_trials_forg, cue_loc_index_all_forg] = genTrials(cue_percent, cue_validity, nTrials*nBlock/2, 2, 0); % remember trials
if mod(str2double(subject),2) == 1 % odd subID - ABBA; even subID - BAAB
    cueCol(1:nTrials*2,:) = cueCol_rem(1:nTrials*2, :);
    probe_loc(1:nTrials*2) = probe_loc_rem(1:nTrials*2);
    test_remember(1:nTrials*2) = test_remember_rem(1:nTrials*2);
    test_forget(1:nTrials*2) = test_forget_rem(1:nTrials*2);
    test_neutral(1:nTrials*2) = test_neutral_rem(1:nTrials*2);
    cue_trials(1:nTrials*2) = cue_trials_rem(1:nTrials*2);
    cue_loc_idx(1:nTrials*2, :) = cue_loc_index_all_rem(1:nTrials*2, :);
    cueCol(nTrials*2+1:nTrials*6,:) = cueCol_forg;
    probe_loc(nTrials*2+1:nTrials*6) = probe_loc_forg;
    test_remember(nTrials*2+1:nTrials*6) = test_remember_forg;
    test_forget(nTrials*2+1:nTrials*6) = test_forget_forg;
    test_neutral(nTrials*2+1:nTrials*6) = test_neutral_forg;
    cue_trials(nTrials*2+1:nTrials*6) = cue_trials_forg;
    cue_loc_idx(nTrials*2+1:nTrials*6,:) = cue_loc_index_all_forg;
    cueCol(nTrials*6+1:nTrials*8,:) = cueCol_rem(nTrials*2+1:end, :);
    probe_loc(nTrials*6+1:nTrials*8) = probe_loc_rem(nTrials*2+1:end);
    test_remember(nTrials*6+1:nTrials*8) = test_remember_rem(nTrials*2+1:end);
    test_forget(nTrials*6+1:nTrials*8) = test_forget_rem(nTrials*2+1:end);
    test_neutral(nTrials*6+1:nTrials*8) = test_neutral_rem(nTrials*2+1:end);
    cue_trials(nTrials*6+1:nTrials*8) = cue_trials_rem(nTrials*2+1:end);
    cue_loc_idx(nTrials*6+1:nTrials*8, :) = cue_loc_index_all_rem(nTrials*2+1:end, :);
else
    cueCol(1:nTrials*2,:) = cueCol_forg(1:nTrials*2, :);
    probe_loc(1:nTrials*2) = probe_loc_forg(1:nTrials*2);
    test_remember(1:nTrials*2) = test_remember_forg(1:nTrials*2);
    test_forget(1:nTrials*2) = test_forget_forg(1:nTrials*2);
    test_neutral(1:nTrials*2) = test_neutral_forg(1:nTrials*2);
    cue_trials(1:nTrials*2) = cue_trials_forg(1:nTrials*2);
    cue_loc_idx(1:nTrials*2, :) = cue_loc_index_all_forg(1:nTrials*2, :);
    cueCol(nTrials*2+1:nTrials*6,:) = cueCol_rem;
    probe_loc(nTrials*2+1:nTrials*6) = probe_loc_rem;
    test_remember(nTrials*2+1:nTrials*6) = test_remember_rem;
    test_forget(nTrials*2+1:nTrials*6) = test_forget_rem;
    test_neutral(nTrials*2+1:nTrials*6) = test_neutral_rem;
    cue_trials(nTrials*2+1:nTrials*6) = cue_trials_rem;
    cue_loc_idx(nTrials*2+1:nTrials*6,:) = cue_loc_index_all_rem;
    cueCol(nTrials*6+1:nTrials*8,:) = cueCol_forg(nTrials*2+1:end, :);
    probe_loc(nTrials*6+1:nTrials*8) = probe_loc_forg(nTrials*2+1:end);
    test_remember(nTrials*6+1:nTrials*8) = test_remember_forg(nTrials*2+1:end);
    test_forget(nTrials*6+1:nTrials*8) = test_forget_forg(nTrials*2+1:end);
    test_neutral(nTrials*6+1:nTrials*8) = test_neutral_forg(nTrials*2+1:end);
    cue_trials(nTrials*6+1:nTrials*8) = cue_trials_forg(nTrials*2+1:end);
    cue_loc_idx(nTrials*6+1:nTrials*8, :) = cue_loc_index_all_forg(nTrials*2+1:end, :);
end
%% Main
quit = 0; block_no = 0; mouseB = [0, 0];
for trial = 1:nTrials*nBlock
    if mod(trial,nTrials) == 1
        block_no = block_no + 1;
        if block_no == 1 || block_no == 3 %|| block_no == 7
            practice = 1;
            rt = ones(numPrac,1)*tresp;
            estimate = NaN(numPrac,1);
            error = 999*ones(numPrac,1);
            points = zeros(numPrac,1);
            target = 999*ones(numPrac,1);
            all_Orientation = zeros(numPrac,3);
            begin = 0; ShowCursor('Arrow'); 
            while begin == 0 && quit == 0
                DrawFormattedText(window,'Your task is to draw the orientation of the cued triangle', 'center', 150, white);
                if order(block_no) == 1
                    DrawFormattedText(window,'GREEN ITEM IS TWICE LIKELY TO BE TESTED', 'center', 250, red);
                elseif order(block_no) == 2
                    DrawFormattedText(window,'GREEN ITEM IS TWICE LIKELY TO BE TESTED, RED ITEM DOES NOT NEED RESPONSE AND CAN BE FORGOT', 'center', 250, red);
                end
                DrawFormattedText(window,'Click on the middle to continue', 'center', 350, white);
                DrawFormattedText(window,'You will now begin with some practice trials', 'center', 450, white);
                Screen('FillOval', window, white, makeBox(windowRect(3)/2,windowRect(4)/2,fixSize))
                Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
                Screen('Flip', window);
    
                while any(mouseB) % if already down, wait for release
                    [mouseX,mouseY,mouseB] = GetMouse(window);
                end
                while ~any(mouseB) % wait for press
                    [mouseX,mouseY,mouseB] = GetMouse(window);
                    centerDist = ((mouseX-windowRect(3)/2)^2 + (mouseY-windowRect(4)/2)^2)^0.5;
                    if centerDist < fixSize/2
                        if mouseB(1) == 1
                            begin = 1;
                        end
                    end
                end


                [keyIsDown,secs, keyCode] = KbCheck;
                if keyCode(escapeKey)
                    quit = 1;
                    fclose(bevID);
                    cd(mainPath);
                    sca; Priority(0);
                end
            end
    
            % practice
            [prec_cueCol, prec_probe_loc, ~, prac_test_forget, ~, ~, ~] = genTrials(3/4, cue_validity, numPrac, order(block_no),practice);
            for prac_trial = 1:numPrac
                %% Fix
                frame = 0; HideCursor();
                while  (quit == 0) && (frame < ffix)
                    frame = frame + 1;
                    Screen('FillArc',window,gray,[xCenter-yCenter, 0, xCenter+yCenter, 2*yCenter],0,360);
                    drawFix(black,fixSize,windowRect,window);
                    Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
                    vbl = Screen('Flip', window);
                    [keyIsDown,secs, keyCode] = KbCheck;
                    if keyCode(escapeKey)
                        quit = 1;
                    end
                end
                %% Stim
                index_Gabor = randperm(length(orientation));
                all_Orientation(prac_trial,:) = [orientation(index_Gabor(1))+jitter*rand(), orientation(index_Gabor(2))+jitter*rand(), orientation(index_Gabor(3))+jitter*rand()];
                frame = 0;
                while  (quit == 0) && (frame < fstim)
                    frame = frame + 1;
                    Screen('FillArc',window,gray,[xCenter-yCenter, 0, xCenter+yCenter, 2*yCenter],0,360);
                    drawFix(black,fixSize,windowRect,window);
                    for ii = 1:3
                        iCenter = stimCenter(ii,:);
                        iangle = deg2rad(all_Orientation(prac_trial,ii));
                        drawPoints = [iCenter(1)+StimHeightPix/2*sin(iangle),iCenter(2)-StimHeightPix/2*cos(iangle); ...
                                iCenter(1)-StimHeightPix/2*sin(iangle)-StimWidthPix/2*cos(iangle),iCenter(2)+StimHeightPix/2*cos(iangle)-StimWidthPix/2*sin(iangle); ...
                                iCenter(1)-StimHeightPix/2*sin(iangle)+StimWidthPix/2*cos(iangle),iCenter(2)+StimHeightPix/2*cos(iangle)+StimWidthPix/2*sin(iangle)];
                        Screen('FillPoly', window, white, drawPoints);
                        % Screen('DrawTexture', window, gaborTex,[], gabLoc(ii,:), all_Orientation(prac_trial,ii));
                    end
                    Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
                    
                    vbl = Screen('Flip', window);
                    if frame == 1
                        stimStart  = vbl - exptStart;
                    end
                    
                    [keyIsDown,secs, keyCode] = KbCheck;
                    if keyCode(escapeKey)
                        quit = 1;
                    end
                end

                %% precue
                frame = 0;
                while  (quit == 0) && (frame < fprecue)
                    frame = frame + 1;
                    Screen('FillArc',window,gray,[xCenter-yCenter, 0, xCenter+yCenter, 2*yCenter],0,360);
                    drawFix(black,fixSize,windowRect,window);
                    Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
                    vbl = Screen('Flip', window);
                    [keyIsDown,secs, keyCode] = KbCheck;
                    if keyCode(escapeKey)
                        quit = 1;
                    end
                end
    
                %% cue
                frame = 0;
                while  (quit == 0) && (frame < fcue)
                    frame = frame + 1;
                    Screen('FillArc',window,gray,[xCenter-yCenter, 0, xCenter+yCenter, 2*yCenter],0,360);
                    for ii = 1:3
                        Screen('DrawLine',window,prec_cueCol{prac_trial}(ii,:),cueLoc(ii,1), cueLoc(ii,2), cueLoc(ii,3), cueLoc(ii,4), cueWidth);
                    end
                    Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
                    vbl = Screen('Flip', window);
                    [keyIsDown,secs, keyCode] = KbCheck;
                    if keyCode(escapeKey)
                        quit = 1;
                    end
                end
    
                %% delay
                frame = 0;
                while  (quit == 0) && (frame < fdelay)
                    frame = frame + 1;
                    Screen('FillArc',window,gray,[xCenter-yCenter, 0, xCenter+yCenter, 2*yCenter],0,360);
                    drawFix(black,fixSize,windowRect,window);
                    Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
                    vbl = Screen('Flip', window);
                    [keyIsDown,secs, keyCode] = KbCheck;
                    if keyCode(escapeKey)
                        quit = 1;
                    end
                end
    
                %% response
                center = stimCenter(prec_probe_loc(prac_trial),:);
                target(prac_trial) = all_Orientation(prac_trial,prec_probe_loc(prac_trial));
                ShowCursor('Arrow');
                if (order(block_no) == 2) && (prac_test_forget(prac_trial) == 1) % forget trial in forget block
                    frame = 0; cursorReset = 0;
                    while  (quit == 0) && (cursorReset ==0) && (debug == 0)
                        frame = frame + 1;
                        Screen('FillArc',window,gray,[xCenter-yCenter, 0, xCenter+yCenter, 2*yCenter],0,360);
                        drawFix(black,fixSize,windowRect,window);
                        DrawFormattedText(window, 'This trial does not need response.', 'center', 210, black);
                        DrawFormattedText(window, 'Clik the center to start next trial', 'center', 250, black);
                        Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
                        vbl = Screen('Flip', window);
                        [mouseX,mouseY,mouseB] = GetMouse(window);
                        dist = (((windowRect(3)/2 - mouseX)^2) + ((windowRect(4)/2 - mouseY)^2))^0.5;
                        if dist < fixSize/2 && mouseB(1) == 1
                            cursorReset = 1;
                        end
                        
                        [keyIsDown,secs, keyCode] = KbCheck;
                        if keyCode(escapeKey)
                            quit = 1;
                        end

                    end
                else 
                    frame = 0; clicks = 0; time0=GetSecs;
                    while (quit ==0 ) && (frame < fresp) && (clicks ==0)
    
                        frame = frame+1;
                        ShowCursor;
                        [mouseX,mouseY,mouseB] = GetMouse(window);
                        deltaX = mouseX-center(1); deltaY = mouseY-center(2);
                        if deltaX > 0 && deltaY < 0
                            angle = -rad2deg(atan(deltaX/deltaY));
                        elseif deltaX > 0 && deltaY >= 0
                            angle = 90+rad2deg(atan(deltaY/deltaX));
                        elseif deltaX <= 0 && deltaY > 0
                            angle = 180-rad2deg(atan(deltaX/deltaY));
                        elseif deltaX < 0 && deltaY < 0
                            angle = 270+rad2deg(atan(deltaY/deltaX));
                        end
            
                        if mouseB(1) == 0
                            Screen('FillArc',window,gray,[xCenter-yCenter, 0, xCenter+yCenter, 2*yCenter],0,360);
                            drawFix(black,fixSize,windowRect,window);
                            drawPoints = [center(1)+StimHeightPix/2*sin(deg2rad(angle)),center(2)-StimHeightPix/2*cos(deg2rad(angle)); ...
                                    center(1)-StimHeightPix/2*sin(deg2rad(angle))-StimWidthPix/2*cos(deg2rad(angle)),center(2)+StimHeightPix/2*cos(deg2rad(angle))-StimWidthPix/2*sin(deg2rad(angle)); ...
                                    center(1)-StimHeightPix/2*sin(deg2rad(angle))+StimWidthPix/2*cos(deg2rad(angle)),center(2)+StimHeightPix/2*cos(deg2rad(angle))+StimWidthPix/2*sin(deg2rad(angle))];
                            Screen('FillPoly', window, white, drawPoints);
                            % Screen('DrawTexture', window, gaborTex,[], gabLoc(prec_probe_loc(prac_trial),:), angle);
                            Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
                            vbl = Screen('Flip', window);
            
                        elseif mouseB(1) == 1
                            estimate(prac_trial) = angle;
                            rt(prac_trial) = GetSecs-time0;
                            error(prac_trial) = rad2deg(circ_dist(deg2rad(estimate(prac_trial)), deg2rad(target(prac_trial))));
                            if abs(error(prac_trial)) < 3
                                points(prac_trial) = 40;
                            elseif abs(error(prac_trial)) >= 3 && abs(error(prac_trial)) < 12
                                points(prac_trial) = 30;
                            elseif abs(error(prac_trial)) >= 12 && abs(error(prac_trial)) < 21
                                points(prac_trial) = 20;
                            elseif abs(error(prac_trial)) >= 21 && abs(error(prac_trial)) <= 30
                                points(prac_trial) = 10;
                            else
                                points(prac_trial) = 0;
                            end
                            clicks = 1;
                            % HideCursor;
                        end
    
                        [keyIsDown,secs, keyCode] = KbCheck;
                        if keyCode(escapeKey)
                            quit = 1;
                        end
    
                    end
                end
    
                %% feedback
    
                if (order(block_no) == 2) && (prac_test_forget(prac_trial) == 1) % forget trial in forget block
                else
                    frame = 0; cursorReset = 0;
                    while  (cursorReset == 0) && (quit == 0) && (debug == 0)
                        frame = frame + 1;
                        Screen('FillArc',window,gray,[xCenter-yCenter, 0, xCenter+yCenter, 2*yCenter],0,360);
                        drawFix(black,fixSize,windowRect,window);
                            if isnan(estimate(prac_trial))
                                DrawFormattedText(window, sprintf('Total Points =  %.0f', sum(points(1:prac_trial))), 'center',250, black);
                                DrawFormattedText(window, 'You missed this trial', 'center', 210, black);
                                drawPoints = [center(1)+StimHeightPix/2*sin(deg2rad(angle)),center(2)-StimHeightPix/2*cos(deg2rad(angle)); ...
                                        center(1)-StimHeightPix/2*sin(deg2rad(angle))-StimWidthPix/2*cos(deg2rad(angle)),center(2)+StimHeightPix/2*cos(deg2rad(angle))-StimWidthPix/2*sin(deg2rad(angle)); ...
                                        center(1)-StimHeightPix/2*sin(deg2rad(angle))+StimWidthPix/2*cos(deg2rad(angle)),center(2)+StimHeightPix/2*cos(deg2rad(angle))+StimWidthPix/2*sin(deg2rad(angle))];
                                Screen('FillPoly', window, white, drawPoints);
                                drawPoints = [center(1)+StimHeightPix/2*sin(deg2rad(target(prac_trial))),center(2)-StimHeightPix/2*cos(deg2rad(target(prac_trial))); ...
                                        center(1)-StimHeightPix/2*sin(deg2rad(target(prac_trial)))-StimWidthPix/2*cos(deg2rad(target(prac_trial))),center(2)+StimHeightPix/2*cos(deg2rad(target(prac_trial)))-StimWidthPix/2*sin(deg2rad(target(prac_trial))); ...
                                        center(1)-StimHeightPix/2*sin(deg2rad(target(prac_trial)))+StimWidthPix/2*cos(deg2rad(target(prac_trial))),center(2)+StimHeightPix/2*cos(deg2rad(target(prac_trial)))+StimWidthPix/2*sin(deg2rad(target(prac_trial)))];
                                Screen('FillPoly', window, [0,128,0], drawPoints);
                            else
                                DrawFormattedText(window, sprintf('Total Points =  %.0f', sum(points(1:prac_trial))), 'center',250, black);
                                % DrawFormattedText(window, sprintf('Your error is %.0f deg.', error(prac_trial)), 'center', 140, black);
                                DrawFormattedText(window, sprintf('You just gained %.0f points!', points(prac_trial)), 'center', 210, black);
                                drawPoints = [center(1)+StimHeightPix/2*sin(deg2rad(angle)),center(2)-StimHeightPix/2*cos(deg2rad(angle)); ...
                                        center(1)-StimHeightPix/2*sin(deg2rad(angle))-StimWidthPix/2*cos(deg2rad(angle)),center(2)+StimHeightPix/2*cos(deg2rad(angle))-StimWidthPix/2*sin(deg2rad(angle)); ...
                                        center(1)-StimHeightPix/2*sin(deg2rad(angle))+StimWidthPix/2*cos(deg2rad(angle)),center(2)+StimHeightPix/2*cos(deg2rad(angle))+StimWidthPix/2*sin(deg2rad(angle))];
                                Screen('FillPoly', window, white, drawPoints);
                                drawPoints = [center(1)+StimHeightPix/2*sin(deg2rad(target(prac_trial))),center(2)-StimHeightPix/2*cos(deg2rad(target(prac_trial))); ...
                                        center(1)-StimHeightPix/2*sin(deg2rad(target(prac_trial)))-StimWidthPix/2*cos(deg2rad(target(prac_trial))),center(2)+StimHeightPix/2*cos(deg2rad(target(prac_trial)))-StimWidthPix/2*sin(deg2rad(target(prac_trial))); ...
                                        center(1)-StimHeightPix/2*sin(deg2rad(target(prac_trial)))+StimWidthPix/2*cos(deg2rad(target(prac_trial))),center(2)+StimHeightPix/2*cos(deg2rad(target(prac_trial)))+StimWidthPix/2*sin(deg2rad(target(prac_trial)))];
                                Screen('FillPoly', window, [0,128,0], drawPoints);
                            end
                        Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
                        vbl = Screen('Flip', window);
                        [mouseX,mouseY,mouseB] = GetMouse(window);
                        dist = (((windowRect(3)/2 - mouseX)^2) + ((windowRect(4)/2 - mouseY)^2))^0.5;
                        if dist < fixSize/2
                            cursorReset = 1;
                        end
                        
                        [keyIsDown,secs, keyCode] = KbCheck;
                        if keyCode(escapeKey)
                            quit = 1;
                        end
                    end
                end
             end
            
        end

    %% start actual trial
        practice = 0;
        ShowCursor('Arrow');
        begin = 0;
        while begin == 0 && quit == 0
                if order(block_no) == 1
                    DrawFormattedText(window,'GREEN ITEM IS TWICE LIKELY TO BE TESTED', 'center', 250, red);
                elseif order(block_no) == 2
                    DrawFormattedText(window,'GREEN ITEM IS TWICE LIKELY TO BE TESTED, RED ITEM DOES NOT NEED RESPONSE AND CAN BE FORGOT', 'center', 250, red);
                end
            DrawFormattedText(window, sprintf('Click on the middle to start block %.0f', block_no), 'center', 350, white);
            % DrawFormattedText(window,sprintf('You completed %.0f of %.0f blocks', block_no-1, nBlock), 'center', 450, white);
            Screen('FillOval', window, white, makeBox(windowRect(3)/2,windowRect(4)/2,fixSize))
            Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
            Screen('Flip', window);

            while any(mouseB) % if already down, wait for release
                [mouseX,mouseY,mouseB] = GetMouse(window);
            end
            while ~any(mouseB) % wait for press
                [mouseX,mouseY,mouseB] = GetMouse(window);
                centerDist = ((mouseX-windowRect(3)/2)^2 + (mouseY-windowRect(4)/2)^2)^0.5;
                if centerDist < fixSize/2
                    if mouseB(1) == 1
                        begin = 1;
                    end
                end
            end


            [keyIsDown,secs, keyCode] = KbCheck;
            if keyCode(escapeKey)
                quit = 1;
                fclose(bevID);
                cd(mainPath);
                sca; Priority(0);
            end
        end

    end


    
    %% Fix
    frame = 0; HideCursor();
    while  (quit == 0) && (frame < ffix)
        frame = frame + 1;
        Screen('FillArc',window,gray,[xCenter-yCenter, 0, xCenter+yCenter, 2*yCenter],0,360);
        drawFix(black,fixSize,windowRect,window);
        Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
        vbl = Screen('Flip', window);
        [keyIsDown,secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            quit = 1;
        end
    end
    %% Stim
    index_Gabor = randperm(length(orientation));
    All_Orientation(trial,:) = [orientation(index_Gabor(1))+jitter*rand(), orientation(index_Gabor(2))+jitter*rand(), orientation(index_Gabor(3))+jitter*rand()];
    frame = 0;
    while  (quit == 0) && (frame < fstim)
        frame = frame + 1;
        Screen('FillArc',window,gray,[xCenter-yCenter, 0, xCenter+yCenter, 2*yCenter],0,360);
        drawFix(black,fixSize,windowRect,window);
        for ii = 1:3
            iCenter = stimCenter(ii,:);
            iangle = deg2rad(All_Orientation(trial,ii));
            drawPoints = [iCenter(1)+StimHeightPix/2*sin(iangle),iCenter(2)-StimHeightPix/2*cos(iangle); ...
                    iCenter(1)-StimHeightPix/2*sin(iangle)-StimWidthPix/2*cos(iangle),iCenter(2)+StimHeightPix/2*cos(iangle)-StimWidthPix/2*sin(iangle); ...
                    iCenter(1)-StimHeightPix/2*sin(iangle)+StimWidthPix/2*cos(iangle),iCenter(2)+StimHeightPix/2*cos(iangle)+StimWidthPix/2*sin(iangle)];
            Screen('FillPoly', window, white, drawPoints);
            % Screen('DrawTexture', window, gaborTex,[], gabLoc(ii,:), All_Orientation(trial,ii));
        end
        Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
        
        vbl = Screen('Flip', window);
        if frame == 1
            stimStart  = vbl - exptStart;
        end
        
        [keyIsDown,secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            quit = 1;
        end
    end

    %% precue
    frame = 0;
    while  (quit == 0) && (frame < fprecue)
        frame = frame + 1;
        Screen('FillArc',window,gray,[xCenter-yCenter, 0, xCenter+yCenter, 2*yCenter],0,360);
        drawFix(black,fixSize,windowRect,window);
        Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
        vbl = Screen('Flip', window);
        [keyIsDown,secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            quit = 1;
        end
    end

    %% cue
    frame = 0;
    while  (quit == 0) && (frame < fcue)
        frame = frame + 1;
        Screen('FillArc',window,gray,[xCenter-yCenter, 0, xCenter+yCenter, 2*yCenter],0,360);
        for ii = 1:3
            Screen('DrawLine',window,cueCol{trial}(ii,:),cueLoc(ii,1), cueLoc(ii,2), cueLoc(ii,3), cueLoc(ii,4), cueWidth);
        end
        Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
        vbl = Screen('Flip', window);
        [keyIsDown,secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            quit = 1;
        end
    end

    %% delay
    frame = 0;
    while  (quit == 0) && (frame < fdelay)
        frame = frame + 1;
        Screen('FillArc',window,gray,[xCenter-yCenter, 0, xCenter+yCenter, 2*yCenter],0,360);
        drawFix(black,fixSize,windowRect,window);
        Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
        vbl = Screen('Flip', window);
        [keyIsDown,secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            quit = 1;
        end
    end

    %% response
    center = stimCenter(probe_loc(trial),:);
    Target(trial) = All_Orientation(trial,probe_loc(trial));
    if (order(block_no) == 2) && (test_forget(trial) == 1) % forget trial in forget block
        frame = 0; cursorReset = 0;
        while  (quit == 0) && (debug == 0) && (cursorReset ==0)
            frame = frame + 1;
            ShowCursor;
            Screen('FillArc',window,gray,[xCenter-yCenter, 0, xCenter+yCenter, 2*yCenter],0,360);
            drawFix(black,fixSize,windowRect,window);
            DrawFormattedText(window, 'This trial does not need response.', 'center', 210, black);
            DrawFormattedText(window, 'Clik the center to start next trial', 'center', 250, black);
            Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
            vbl = Screen('Flip', window);
            [keyIsDown,secs, keyCode] = KbCheck;
            [mouseX,mouseY,mouseB] = GetMouse(window);
            dist = (((windowRect(3)/2 - mouseX)^2) + ((windowRect(4)/2 - mouseY)^2))^0.5;
            if dist < fixSize/2 && mouseB(1) == 1
                cursorReset = 1;
            end

            if keyCode(escapeKey)
                quit = 1;
            end
        end
    else 
        frame = 0; clicks = 0; time0=GetSecs;
        while (quit ==0 ) && (frame < fresp) && (clicks ==0)

            frame = frame+1;
            ShowCursor;
            [mouseX,mouseY,mouseB] = GetMouse(window);
            deltaX = mouseX-center(1); deltaY = mouseY-center(2);
            if deltaX > 0 && deltaY < 0
                angle = -rad2deg(atan(deltaX/deltaY));
            elseif deltaX > 0 && deltaY >= 0
                angle = 90+rad2deg(atan(deltaY/deltaX));
            elseif deltaX <= 0 && deltaY > 0
                angle = 180-rad2deg(atan(deltaX/deltaY));
            elseif deltaX < 0 && deltaY < 0
                angle = 270+rad2deg(atan(deltaY/deltaX));
            end

            if mouseB(1) == 0
                Screen('FillArc',window,gray,[xCenter-yCenter, 0, xCenter+yCenter, 2*yCenter],0,360);
                drawFix(black,fixSize,windowRect,window);
                drawPoints = [center(1)+StimHeightPix/2*sin(deg2rad(angle)),center(2)-StimHeightPix/2*cos(deg2rad(angle)); ...
                        center(1)-StimHeightPix/2*sin(deg2rad(angle))-StimWidthPix/2*cos(deg2rad(angle)),center(2)+StimHeightPix/2*cos(deg2rad(angle))-StimWidthPix/2*sin(deg2rad(angle)); ...
                        center(1)-StimHeightPix/2*sin(deg2rad(angle))+StimWidthPix/2*cos(deg2rad(angle)),center(2)+StimHeightPix/2*cos(deg2rad(angle))+StimWidthPix/2*sin(deg2rad(angle))];
                Screen('FillPoly', window, white, drawPoints);
                % Screen('DrawTexture', window, gaborTex,[], gabLoc(probe_loc(trial),:), angle);
                Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
                vbl = Screen('Flip', window);

            elseif mouseB(1) == 1
                Estimate(trial) = angle;
                RT(trial) = GetSecs-time0;
                Error(trial) = rad2deg(circ_dist(deg2rad(Estimate(trial)), deg2rad(Target(trial)))); % positive: cw estimate, negative: ccw estimate
                if abs(Error(trial)) < 3
                    Points(trial) = 40;
                elseif abs(Error(trial)) >= 3 && abs(Error(trial)) < 12
                    Points(trial) = 30;
                elseif abs(Error(trial)) >= 12 && abs(Error(trial)) < 21
                    Points(trial) = 20;
                elseif abs(Error(trial)) >= 21 && abs(Error(trial)) <= 30
                    Points(trial) = 10;
                else
                    Points(trial) = 0;
                end
                clicks = 1;
                % HideCursor;
            end

            [keyIsDown,secs, keyCode] = KbCheck;
            if keyCode(escapeKey)
                quit = 1;
            end

        end
    end
    %% Save
    if quit == 0
        datastring = strcat(subject, ',', session, ',', num2str(block_no), ',' ,num2str(order(block_no)), ',' ,num2str(trial), ',' , num2str(cue_trials(trial)), ',' , num2str(All_Orientation(trial,1)) , ',', num2str(All_Orientation(trial,2)), ',', num2str(All_Orientation(trial,3)), ',',  ...
            num2str(cue_loc_idx(trial,1)), ',', num2str(cue_loc_idx(trial,2)), ',', num2str(cue_loc_idx(trial,3)) , ',', num2str(probe_loc(trial)), ',', num2str(test_remember(trial)), ',', num2str(test_forget(trial)), ',', num2str(test_neutral(trial)), ',', ...
            num2str(Target(trial)), ',',  num2str(Estimate(trial)), ',',  num2str(RT(trial)), ',', num2str(Error(trial)), ',',  num2str(Points(trial)),',',  '\n');
        fprintf(bevID,datastring);
    end
    %% feedback
    if (order(block_no) == 2) && (test_forget(trial) == 1) % forget trial in forget block
    else
        frame = 0; cursorReset = 0;
        while  (cursorReset == 0) && (quit == 0) && (debug == 0)
            frame = frame + 1;
            Screen('FillArc',window,gray,[xCenter-yCenter, 0, xCenter+yCenter, 2*yCenter],0,360);
            drawFix(black,fixSize,windowRect,window);
                if isnan(Estimate(trial))
                    DrawFormattedText(window, sprintf('Total Points =  %.0f', sum(Points(1:trial))), 'center',250, black);
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
                    DrawFormattedText(window, sprintf('Total Points =  %.0f', sum(Points(1:trial))), 'center',250, black);
                    % DrawFormattedText(window, sprintf('Your error is %.0f', Error(trial)), 'center', 140, black);
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
            Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
            vbl = Screen('Flip', window);
    
            [mouseX,mouseY,mouseB] = GetMouse(window);
            dist = (((windowRect(3)/2 - mouseX)^2) + ((windowRect(4)/2 - mouseY)^2))^0.5;
            if dist < fixSize/2
                cursorReset = 1;
            end
            
            [keyIsDown,secs, keyCode] = KbCheck;
            if keyCode(escapeKey)
                quit = 1;
            end
        end
    end
    %% end of block
    if mod(trial,nTrials) ==0 && trial ~= nTrials*nBlock
        begin = 0;  ShowCursor;
        while begin == 0 && quit == 0

            DrawFormattedText(window,sprintf('You completed %.0f of %.0f blocks', block_no, nBlock), 'center', 100, white);
            DrawFormattedText(window,'Please take a short break', 'center', 200, white);
            DrawFormattedText(window,'Click the center to continue', 'center', 300, white);
            Screen('FillOval', window, white, makeBox(windowRect(3)/2,windowRect(4)/2,fixSize))
            Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
            Screen('Flip', window);
            
            while any(mouseB) % if already down, wait for release
                [mouseX,mouseY,mouseB] = GetMouse(window);
            end
            while ~any(mouseB) % wait for press
                [mouseX,mouseY,mouseB] = GetMouse(window);
                centerDist = ((mouseX-windowRect(3)/2)^2 + (mouseY-windowRect(4)/2)^2)^0.5;
                if centerDist < fixSize/2
                    if mouseB(1) == 1
                        begin = 1;
                    end
                end
            end

        
            [keyIsDown,secs, keyCode] = KbCheck;
            if keyCode(escapeKey)
                quit = 1;
                fclose(bevID);
                cd(mainPath);
                sca; Priority(0);
            end
        end
    elseif trial == nTrials*nBlock
        duration = GetSecs - startTime
        Screen('Flip', window); % blank
        disp('Experiment successfully completed!')
        while quit ==0
            DrawFormattedText(window,'You have completed the experiment.', 'center', 200, white);
            DrawFormattedText(window,'Please find the experimenter.', 'center', 300, white);
            Screen('DrawingFinished', window);  % Tell PTB no more drawing commands will be issued until the next flip
            Screen('Flip', window);
            
            [keyIsDown,secs, keyCode] = KbCheck;
            if keyCode(escapeKey)
                quit = 1;
            end
            
        end
    end

end

%% End
fclose(bevID); 
cd(mainPath);
%bonus = min(maxBonus, round(50 * sum(Points) / (40*(nTrials*nBlock))/5)*5);
bonus = round(sum(points)/(40*nTrials*nBlock)*maxBonus);
fprintf('*******************\nBonus = %d AED\n*******************\n', bonus);
sca;
Priority(0);


%% helper functions
function drawFix(color,fixSize,windowRect,window)
    Screen('FillArc', window, color, [windowRect(3)/2 - fixSize, windowRect(4)/2 - fixSize , windowRect(3)/2 + fixSize, windowRect(4)/2 + fixSize], 0, 360);
end

function box = makeBox(x,y,r)
% makes a square rect with xy in the center, 'radius' (half-height) of r
box = [x-r, y-r, x+r, y+r];
end


function [GaborP, Gabor, mask] = genGabor(degToPix)
    % make gabor
    GaborP.contrast = 0.6;
    GaborP.SF = 0.8; % cyc/deg
    GaborP.circleDiameterDeg = 6; % degree
    GaborP.circleDiameterPix = GaborP.circleDiameterDeg * degToPix;
    GaborP.roughness = 1; % pixels
    GaborP.phaseRad = 2; % radians
    GaborP.theta = 0; % starts from upright
    GaborP.sigmaDig = GaborP.circleDiameterDeg/5; 
    GaborP.sigmaPix = GaborP.sigmaDig * degToPix;
    Gabor = MakeGaborTex(GaborP.contrast,degToPix,GaborP.SF,GaborP.circleDiameterPix,GaborP.roughness,GaborP.phaseRad,GaborP.theta,GaborP.sigmaPix);
    mask = concentricCirclesTexture(0.5,GaborP.circleDiameterPix,GaborP.SF,GaborP.sigmaPix,degToPix,GaborP.roughness);
end


function [cueCol, probe_loc, test_remember, test_forget, test_neutral, cue_trials, cue_loc_index_all] = genTrials(cue_percent, cue_validity, nTrials, condition, practice)
cue_trials = [ones(1,nTrials*cue_percent), zeros(1,round(nTrials*(1-cue_percent)))];
cue_loc_index = [1,2,3; 1,3,2; 2,1,3; 2,3,1; 3,1,2; 3,2,1]; % [remember, forget, neutral], 1: left, 2: right, 3: middle; 
if practice == 0
    cue_loc_index = repmat(cue_loc_index, [nTrials*cue_percent/size(cue_loc_index,1),1]);
else
    cue_loc_index = cue_loc_index(1:nTrials*cue_percent, :);
end
cue_loc_index_nocue = zeros(round(nTrials*(1-cue_percent)),3);
cue_loc_index_all = [cue_loc_index; cue_loc_index_nocue];

rand_idx = randperm(nTrials*cue_percent);
test_remember_idx = rand_idx(1:nTrials*cue_percent*cue_validity); % if cue is valid, test in remember item
if practice == 0
    test_forget_idx = rand_idx(nTrials*cue_percent*cue_validity+1:...
        nTrials*cue_percent*cue_validity+nTrials*cue_percent*(1-cue_validity)/2); % if cue is not valid, sometime test forget item
    test_neutral_idx = rand_idx(nTrials*cue_percent*cue_validity+nTrials*cue_percent*(1-cue_validity)/2+1:end); % sometime test neutral item
else
    test_forget_idx = rand_idx(nTrials*cue_percent*cue_validity+1:...
        nTrials*cue_percent*cue_validity+round(nTrials*cue_percent*(1-cue_validity)/2)); % if cue is not valid, sometime test forget item
    test_neutral_idx = rand_idx(nTrials*cue_percent*cue_validity+round(nTrials*cue_percent*(1-cue_validity)/2)+1:end); % sometime test neutral item

end
test_remember = zeros(nTrials, 1); % put index = 1 for testing remember trials
test_forget = zeros(nTrials, 1); % put index = 1 for testing forget trials
test_neutral = zeros(nTrials, 1); % put index = 1 for testing neutral trials
for ii = 1:nTrials*cue_percent
    if sum(test_remember_idx == ii) ==1
        probe_loc(ii) = cue_loc_index_all(ii,1); % probe location is at the remember item
        test_remember(ii) = 1;
    elseif sum(test_forget_idx == ii) ==1
        probe_loc(ii) = cue_loc_index_all(ii,2); % probe location is at the forget item
        test_forget(ii) = 1;
    elseif sum(test_neutral_idx == ii) ==1
        probe_loc(ii) = cue_loc_index_all(ii,3);  % probe location is at the forget item
        test_neutral(ii) = 1;
    end
end

% no cue trials, probe at random position
if practice == 0
    nocue_probe_loc = repmat([1,2,3], [1, round(nTrials*(1-cue_percent)/3)]);
else
    nocue_probe_loc = randi(3,1,nTrials*(1-cue_percent));
end
probe_loc(nTrials*cue_percent+1:nTrials) = nocue_probe_loc(randperm(numel(nocue_probe_loc)));

% randomize trials
rand_idx = randperm(nTrials);
cue_trials = cue_trials(rand_idx);
probe_loc = probe_loc(rand_idx);
test_remember = test_remember(rand_idx);
test_neutral = test_neutral(rand_idx);
test_forget = test_forget(rand_idx);
cue_loc_index_all = cue_loc_index_all(rand_idx,:);


% cue colors
colors = [91, 255, 46; 0,0,255; 255,0,0]./255; % remember, neutral, forget
cueCol = cell(nTrials,1);
for ii = 1:nTrials
    if cue_trials(ii) == 1
        remember_loc = cue_loc_index_all(ii,1);
        forget_loc = cue_loc_index_all(ii,2);
        neutral_loc = cue_loc_index_all(ii,3);
        cueCol{ii}(remember_loc,:) = colors(1,:);
        if condition == 1 % only remember cue
            cueCol{ii}(forget_loc,:) = colors(2,:);
            cueCol{ii}(neutral_loc,:) = colors(2,:);
        elseif condition == 2 % remember and one forget cue
            cueCol{ii}(forget_loc,:) = colors(3,:);
            cueCol{ii}(neutral_loc,:) = colors(2,:);
        end
    elseif cue_trials(ii) == 0
        cueCol{ii} = [colors(2,:); colors(2,:); colors(2,:)];
    end
end

end

function mainTexture = concentricCirclesTexture(contrast,circleDiameter, spatialFrequency,sigma, degToPix,roughness)

    lambda = degToPix/(spatialFrequency); 
    trim = 0.05;
    X = 1:roughness:circleDiameter; 
    [xx yy] = meshgrid(X-circleDiameter/2, X-circleDiameter/2); 

    B= double(((xx).^2+(yy).^2 >(circleDiameter/2)^2));
    
    R = sqrt((xx).^2+(yy).^2); %radius 
    k=2*pi/lambda; % wave vector 2*pi/lambda where lambda = distance/numOfCircles
    
    phi=pi/2; % phase 
    
    Z = contrast*sin(k*R + phi);
    B(B == 0) = Z(B == 0);

    
    X = 1:roughness:circleDiameter; 
    X0 = (X / circleDiameter) - .5; 
    [Xm Ym] = meshgrid(X0, X0); 
    
    s = sigma / circleDiameter;
    gauss = exp( -((((Xm).^2)+((Ym).^2)) ./ (2* s^2)) );
    gauss(gauss < trim) = 0; % trim around edges (for 8-bit colour displays)
    %Gauss = ones(size(gauss));
    %(B==0) = gauss(B==0);
    
    BG = B.*gauss;   
    BG = BG + 1;
    
    mainTexture = BG*127.5;%cat(3,Z,Z,Z);
    mainTexture(:,:,2) = gauss.*255;
end 

function mainTexture = MakeGaborTex(contrast,degToPixX,spatialFrequency,circleDiameterPix,roughness,phaseRad,theta,sigma)

    trim = 0.05;
    
    %make gabor
    lamda = degToPixX/(spatialFrequency); % wavelength (number of pixels per degree)
    X = 1:roughness:circleDiameterPix;                           % X is a vector from 1 to imageSize
    X0 = (X / circleDiameterPix) - .5;                 % rescale X -> -.5 to .5
    freq = circleDiameterPix/lamda;                    % compute frequency from wavelength
    Xf = X0 * freq * 2*pi;                  % convert X to radians: 0 -> ( 2*pi * frequency)
  
    sinX = sin( Xf + phaseRad) ;            % make phase-shifted sinewave

    [Xm Ym] = meshgrid(X0, X0);             % 2D matrices
    Xf = Xm * freq * 2*pi;

    
    Xt = Xm * sind(theta);     % compute proportion of Xm for given orientation
    Yt = Ym * cosd(theta);     % compute proportion of Ym for given orientation
    XYt = [ Xt + Yt ];                      % sum X and Y components
    XYf = XYt * freq * 2*pi;                % convert to radians and scale by frequency
    grating = contrast * sin( XYf + phaseRad);                   % make 2D sinewave

    s = sigma / circleDiameterPix;                     % gaussian width as fraction of imageSize
    gauss = exp( -(((Xm.^2)+(Ym.^2)) ./ (2* s^2)) ); % formula for 2D gaussian

    gauss(gauss < trim) = 0;                 % trim around edges (for 8-bit colour displays)

    gabor = grating .* gauss;                % use .* dot-product

    mainTexture = (gabor+1).*127.5;
    mainTexture(:,:,2) = gauss.*255;

    
end
