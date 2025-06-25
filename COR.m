%% Orientation and Confidence Judgments (OCJ)
% Last updated: May 23, 2025
%
% To-do list
% 1. Gamma correction (compare with pldaps code)
% 2. Stimulus parameter unit. Currently, pixels, frames amd sf affect Bv
% 3. Stimulus uncertainty manipulation
% 4. Sound
% 5. Show the number of remaining blocks during inter-block interval
% 6. Spatially mapped orientation demo?
% 7. Check center of the monitor
% Estimates orientation in dat file

%% Start clean
close all;
clearvars;

% Random seed
rng shuffle

%% Use EyeLink?
dummymode = 0; % variable for eyelink usage. (1 for not using eye-tracking)

%% Initialize the script
% Experiment information
% TODO: take subject name and uniquely generate the session no
answer = str2double(inputdlg({'Subject number:','Session number:'},'Experiment Information')); % subject number = 0 for pilot
if sum(isnan(answer))
    error('Invalid input.')
end
subjectNum = answer(1);
sessionNum = answer(2);

% Name the .mat file
matFile = ['COR' num2str(subjectNum,'%.2d') '.mat'];

% Check for incorrect input
if subjectNum ~= 0 % not a pilot
    if sessionNum == 1 
        % This is the first session but the file exists
        if ~isempty(dir(matFile)) 
            error('Incorrect subject or session number.')
        end
    else
        % This is not the first session but the file does not exist
        if isempty(dir(matFile)) 
            error('You''ve skipped the first session.')
        end
        % Load table from the previous session
        load(matFile)
        % There exists data from this session
        if sum( dat.session == sessionNum )
            error('You''ve already done this session.')
        end
    end
end

% TODO: maybe create a different file for different session
% Add instruction screen

%% Table
% Extreme values
% Contrast 0.015, duration 0.06, loc 4.5, spread 5
% Contrast 0.018, duration 0.06, loc 4.5, spread 45 
interTrlInterval = 1;                                                      % Inter trial interval in seconds
fixationDur = 0.5;                                                         % Fixation duration in seconds
stimOrientations = linspace(0, 179, 10);                                   % 
stimLoc_x = 0;                                                             % Stimulus location in visual field degrees
stimLoc_y = 4.5;                                                             % Stimulus location in visual field degrees
stimDur = [0.2, 0.2];  %0.15 0.5. 0.14 is  the minimum                                                    % Stimulus duration in seconds
stimSpread = [45, 45]; % 45                                                       % Stimulus spread in degrees
% stimContrast = [0.015, 0.05];                                                 % Stimulus contrast levels
stimContrast = [0.018, 0.018]; %0.018
respMaxDur = 5;                                                            % Maximum allowed time for user to respond (2 seconds)
respSuccessWaitDur = 0.5;
numBlocks   = 2;                                                           % Number of blocks 
nTrialsPerBlock = numel(stimOrientations)*numel(stimSpread)*numel(stimContrast)*numel(stimDur);    % Assuming each trial takes max of 5 second, a block should take ~8 minutes
nTrials = numBlocks*nTrialsPerBlock;                                       % Total number of trials to run in this session
durFeedback = 1;
beeperDur = 0.15;
durFeedbackFixBreak = .5;
checkFixationMaxDur = 2;
respScreenGazeHoldDur = 0.2;
trlStatus = 0;                                                             % 0 for not completed yet, 1 for successful completion

% TODO: handle scenario when the same subject runs session the next day

% Simple vectors
subjectVec    = subjectNum*ones(nTrials,1);
sessionVec    = sessionNum*ones(nTrials,1); % fill in sesion number when each session starts
blockVec      = reshape(repmat(1:numBlocks, nTrialsPerBlock, 1), nTrials, 1);
trialVec      = repmat((1:nTrialsPerBlock)', numBlocks, 1);
itiVec        = interTrlInterval*ones(nTrials,1);
fixDurVec     = fixationDur*ones(nTrials, 1);
respMaxDurVec = respMaxDur*ones(nTrials, 1);
stimLocVec_x  = repmat(stimLoc_x, nTrials, 1);
stimLocVec_y  = repmat(stimLoc_y, nTrials, 1);
trlStatusVec  = trlStatus*ones(nTrials, 1);

% Factorial design for stimulus
[d, s, c] = ndgrid(stimDur, stimSpread, stimContrast);
combinations = [d(:), s(:), c(:)];  % 8x3

comboRepeated = repmat(combinations, length(stimOrientations), 1);  % Repeat the combinations for each orientation (10 times)
orientationsRepeated = repelem(stimOrientations(:), size(combinations, 1));  % Repeat orientations to match combination matrix
stimMatrix = [orientationsRepeated, comboRepeated]; % Final matrix: each row = [orientation, duration, spread, contrast]
finalStimMatrix = repmat(stimMatrix, numBlocks, 1);

% Shuffle the stim matrix across all the trials (note: shuffling is not within a single block)
shuffledIdx = randperm(size(finalStimMatrix, 1));
finalStimMatrix = finalStimMatrix(shuffledIdx, :);

stimOriVec        = rand(nTrials, 1)*179; % Replace fixed orientations with random orientations sampled from uniform distribution (0 - 179)
% stimOriVec        = finalStimMatrix(:, 1); % Stim orientation vector for all the trials  in this block
stimDurVec        = finalStimMatrix(:, 2); % Stim duration vector for all the trials  in this block
stimDispersionVec = finalStimMatrix(:, 3); % Stim dispersion vector for all the trials  in this block
stimContrastVec   = finalStimMatrix(:, 4); % Stim contrast dispersion vector for all the trials  in this block

stimSampleMeanOriVec   = nan(nTrials,1);
stimSampleSpreadVec    = nan(nTrials,1);
reportedOriVec         = nan(nTrials,1);
reportedConfVec        = nan(nTrials,1);
rewardVec              = nan(nTrials,1);
respTimeVec            = nan(nTrials,1);



% Finalize the table
variableNames = {...
    'subject','session','block','trial', 'ITI',... % TODO: experiment information
    'fixationDur', 'stimLocX', 'stimLocY', 'stimDur', 'stimOri', 'stimSpread', 'stimContrast',... % TODO: independent and confounding variables
    'stimSampleMeanOri', 'stimSampleSpread', ...
    'respMaxDur', 'reportedOri','reportedConf','reward','reactionTime', ... % TODO: dependent variables
    'trialStatus'
    };

dat = table( ...
        subjectVec, sessionVec, blockVec, trialVec, itiVec, ...
        fixDurVec, stimLocVec_x, stimLocVec_y, stimDurVec, stimOriVec, stimDispersionVec, stimContrastVec, ...
        stimSampleMeanOriVec, stimSampleSpreadVec, ...
        respMaxDurVec, reportedOriVec, reportedConfVec, rewardVec, respTimeVec, ...
        trlStatusVec, ...
        'VariableNames', variableNames ...
    );

description = table(categorical(...
    {...
    'subject number';...
    'session number (within each subject)';...
    'block number (within each session)';...
    'trial number (within each block)';...
    'Inmter trial Interval (in seconds)';...
    'fixation duration (in seconds)';...
    'Stimulus duration (in seconds)';...
    'stimulus location: (x) location in visual degrees';...
    'stimulus location (y) location in visual degrees'
    'Stimulus orientation (in degree)';...
    'Stimulus spread (in degrees) - dispersion';...
    'stimulus RMS contrast [sigma]';...
    'sample average orientation (in degree) - sample avg on each trial might be slightly different from actual stimulus orientation used by generative process';...
    'sample orientation spread (in degree) - sample spread on each trial might be slightly different from actual spread used by generative process';...
    'Maximum duration for which response window is shown (in s)';...
    'Reported orientation (in deg)';...
    'Reported confidence (1: Confident (HC), 0: Not confident(LC))';...
    'Reward earned on this trial';...
    'Reaction time (time interval between onset of response screen and choice commit) (in s)';...
    'Trial status - whether trial is completed or not'; ...
    }),...
    'VariableNames',{'description'},'RowNames',variableNames');

%% Reward
baselineReward = 10; % If done over multiple session, then this should probably be reloaded from the previous session
currentTotalReward = baselineReward; % This will not be baseline if a session file for this subject already exist

%%
tSessionStart = GetSecs;
try
    %% Initialize EyeLink
    edfFile = initializeEyeLink(subjectNum, dummymode);
    
    % Start timer
    tStartSession = GetSecs;
    
    % Initialize Psychtoolbox
    psychToolBoxConfig = initPsychToolBox();   
    
    %% Main experiment block goes here

    % Show trial and instruction before the experiment starts
    % Spectral env should be generated before
    % Before the experiment block starts - do following
    % Obtain fixation window (No need to do this in every trial). Just
    % obtain once and be done with it
    fixationWinCfg = initializeFixationWin(psychToolBoxConfig.ppd, ...
        psychToolBoxConfig.xCenter, psychToolBoxConfig.yCenter, ...
        psychToolBoxConfig.black, psychToolBoxConfig.white);

    % Maybe generating the stimuli movies for the factorial design
    % conditions is not such a good idea! Briefly, it would mean that
    % across conditions differ by not just the value of the corresponding
    % factors but also by random phase. Where within conditions phase
    % remains the same across trials. If we find difference across 
    % conditions, it will be hard to discern whether the difference was
    % dur to random phase of the value of the factors itself. One way to
    % rule out this would be to have random phases within condition
    % itself. If certain effect is observed across conditions, then it
    % probably arised becaused of difference in factor values and NOT
    % random phase.

    % Generate stimuli corresponding to stim duration and stimSpread
    % TODO: This needs a better way of handling. Factorial design might 
    % get more complicated.
    % Also, note that using pregenerated stimCfg will only help to avoid
    % any potential issues with non-uniform spacing of the spatial
    % frequency grid.
    oriRefStim = 45; % 45 degrees
    % stimCfgs = cell(numel(stimDur), numel(stimSpread));  % preallocate cell array
    % 
    % for durIdx=1:numel(stimDur)
    %     for spreadIdx=1:numel(stimSpread)
    %         stimDuration_ = stimDur(durIdx);
    %         stimSpread_   = stimSpread(spreadIdx);
    % 
    %         fprintf("Dur: %.2f, spread: %.2f \n", stimDuration_, stimSpread_);
    % 
    %         stimCfg = generateStimuli(psychToolBoxConfig, stimDuration_, stimSpread_, oriRefStim);
    %         stimCfgs{durIdx, spreadIdx} = stimCfg;  % store in cell
    % 
    %         % % Can't use a random seed phase because it seem not possible to
    %         % % have same random phase for all the combinations in factorial
    %         % % design. One of the factor is stim duration, which affects the
    %         % % number of frames in the stimuli. Random phase initialization
    %         % % for the stimuli depends upon the number of frames as well.
    %         % movie = randomCloud( stimCfg.envSpectral, false, [], false, stimCfg.counterphase, [] );
    % 
    %     end
    % end

    for blockIDx=1:numBlocks
        
        % Calibrate eyelink at start of each block
        eyeUsed = setupEyeLink(psychToolBoxConfig.w);
        
        % Wait for user to start the experiment
        beginExpScreen(psychToolBoxConfig)
        Eyelink('Message','START');

        % Dummy value - for ITI calculations
        tEndPrevTrl = GetSecs;
        
        % for trialIDx=1:nTrialsPerBlock
        for trialIDx=1:10
        
            % === Fetch the trial configuration from the dat file ===
            trlCfgIdx = find( (dat.session == sessionNum) & (dat.block == blockIDx) & (dat.trial == trialIDx) );
            trlCfg = dat( (dat.session == sessionNum) & (dat.block == blockIDx) & (dat.trial == trialIDx), :);
            
            % === Generate stimuli ===
            % This done for each trial But do it even before
            % during the inter trial interval
            % Generate stimuli on trial by trial basis
            % How to avoid problem of non uniform orientation grids!!! -
            % Preallocate at 45 degrees i.e. at one angle
            % stimCfg = generateStimuli(psychToolBoxConfig, trlCfg);
            % Get the preallocated stimCfg for this condition
            % idxStimDur = find(stimDur == trlCfg.stimDur, 1);
            % idxStimSpread = find(stimSpread == trlCfg.stimSpread, 1);
            % stimCfg_ = stimCfgs{idxStimDur, idxStimSpread};
            
            % In order to simulate what internal noise does to a stimuli in
            % external noise, the stimuli sample needs to be generated in
            % each trial. Every sample will have its own sample dispersion
            % and mean. Note: the stimuli is still generated at reference
            % stimuli to avoid problems associated with non uniform grid
            % spacing.
            stimCfg = generateStimuli(psychToolBoxConfig, trlCfg.stimDur, trlCfg.stimSpread, oriRefStim);
            
            % Draw contrast dependent stimuli
            % Bote: the stim param has reference orientation of 45 degree
            stimParam = drawStimuli(stimCfg, trlCfg, psychToolBoxConfig);
            % Calculate dispersion and mean orientation for the sample
            % stimuli
            stimMetrics = calcStimMetrics(stimParam.movie);
            fprintf("Actual Stim: angle = %.2f, stdev = %.2f \n", trlCfg.stimOri, trlCfg.stimSpread);
            
            % Inter trial interval - Do whateven preperation needs to be
            % done before this ITI
            WaitSecs( trlCfg.ITI - ( GetSecs - tEndPrevTrl ) );
            
            % === Fixation point ====
            fixnData = showFixation(psychToolBoxConfig, fixationWinCfg);
            
            % Wait untill subject maintains fixation
            fixationHeld = true;
            
            [fixationHeld, eyeUsed] = checkFixation(psychToolBoxConfig, eyeUsed, fixationWinCfg.fixThreshold, ...
                fixnData.tStart, trlCfg.fixationDur, ...
                psychToolBoxConfig.xCenter, psychToolBoxConfig.yCenter, psychToolBoxConfig.quitKey, ...
                checkFixationMaxDur);
            
            if ~fixationHeld
                % TODO: set info in dat file
                % Maye show trial aborted screen as well
                abortTrial(psychToolBoxConfig.w, durFeedbackFixBreak, beeperDur);
                % Everytime trial is aborted, update the timer for ITI calculation
                [~, tEndPrevTrl] = Screen('Flip', psychToolBoxConfig.w); 
                continue;
            end
            
            % === Stimulus Phase ===
            fixationHeld = true;
            
            % Loop over frames
            for iF = 1:length(stimParam.movietex)
                Screen('FillOval', psychToolBoxConfig.w, fixationWinCfg.fixColor, fixationWinCfg.fixRect);
                Screen('DrawTextures', psychToolBoxConfig.w, stimParam.movietex(iF), ...
                    [], stimParam.movieRect, oriRefStim - trlCfg.stimOri);
                Screen('Flip', psychToolBoxConfig.w);
                
                % fprintf("%.2f= ", trlCfg.stimOri - stimMetrics.meanAngle)
                % For each frame make sure fixation is maintained
                % Abort trial if subject breaks fixation
                % Send event to eyelink
                if iF == 1
                    Eyelink('Message','STIMULUS');
                end
                
                % Keep this inside the loop to avoid any timing issues
                % Check fixation
                if Eyelink('NewFloatSampleAvailable') > 0
                    evt = Eyelink('NewestFloatSample');
                    gx = evt.gx(eyeUsed+1); % [left eye gaze x, right eye gaze x] +1 because this is matlab
                    gy = evt.gy(eyeUsed+1);
                    % Abort trial if subject break fixation
                    if norm([psychToolBoxConfig.xCenter, psychToolBoxConfig.yCenter]-[gx,gy]) > fixationWinCfg.fixThreshold
                        fixationHeld = false; 
                        break
                    end
                end
            end
            
            if ~fixationHeld
                % TODO: set info in dat file
                abortTrial(psychToolBoxConfig.w, durFeedbackFixBreak, beeperDur);
                % Everytime trial is aborted, update the timer for ITI calculation
                [~, tEndPrevTrl] = Screen('Flip', psychToolBoxConfig.w);
                continue;
            end
            
            % === Response screen ====
            % Screen('Flip', psychToolBoxConfig.w);
            respData = showResponseScreen(psychToolBoxConfig, fixationWinCfg, ...
                trlCfg, eyeUsed, respScreenGazeHoldDur, beeperDur, respSuccessWaitDur);
            
            if ~respData.responseGiven
                % TODO: set info in dat file
                DrawFormattedText(psychToolBoxConfig.w, 'Timeout occured while waiting for response', ...
                    'center', 'center',[255 0 0]);
                Screen('Flip', psychToolBoxConfig.w);
                
                % All aboted trial will be tried again towards the end
                abortTrial(psychToolBoxConfig.w, durFeedbackFixBreak, beeperDur);
                % Everytime trial is aborted, update the timer for ITI calculation
                [~, tEndPrevTrl] = Screen('Flip', psychToolBoxConfig.w);
                continue;
            end
            
            % This trial was a success. Give reward, show feedback
            % sampleTrueOri = trlCfg.stimOri;
            reportedOri   = respData.reportedAngle;
            confReport    = respData.reportedConf;
            
            % Calculate reward
            reward = calcReward(trlCfg.stimOri, reportedOri, confReport);
            currentTotalReward = currentTotalReward + reward;
            
            if reward == 0
                DrawFormattedText(psychToolBoxConfig.w,['+' num2str(reward,'%.2f')],'center', ...
                    psychToolBoxConfig.yCenter - 30,[255 0 0]);
            else
                DrawFormattedText(psychToolBoxConfig.w,['+' num2str(reward,'%.3f')],'center', ...
                        psychToolBoxConfig.yCenter - 30,[0 255 0]);
            end
            
            DrawFormattedText(psychToolBoxConfig.w, ['$' num2str(currentTotalReward, '%.2f')], ...
                'center', psychToolBoxConfig.yCenter + 30, 255);
            Screen('Flip', psychToolBoxConfig.w);
            Eyelink('Message','FEEDBACK');
            WaitSecs(durFeedback);
            
            % Start timer for inter-trial interval
            [~, tEndPrevTrl] = Screen('Flip', psychToolBoxConfig.w);
            
            % Record everything - Like reaction time, choice and other data
            % in dat file.
            
            % Save info in file
            dat.stimSampleMeanOri(trlCfgIdx)  = trlCfg.stimOri + (stimMetrics.meanAngle - oriRefStim);
            dat.stimSampleSpread(trlCfgIdx)   = stimMetrics.stdAngle;
            dat.reportedOri(trlCfgIdx)        = reportedOri;
            dat.reportedConf(trlCfgIdx)       = confReport;
            dat.reward(trlCfgIdx)             = reward;
            dat.reactionTime(trlCfgIdx)       = respData.reactionTime;
            dat.trialStatus(trlCfgIdx)        = 1; % This trial was completed
            
        end

        % End of block related operations
        
        % Send event to eyelink
        Eyelink('Message', 'END');
    
        % Stop recording
        Eyelink('StopRecording');
    
        % Close the .edf file
        WaitSecs(1);
        Eyelink('CloseFile');

        % TODO: save info and print block stats
        endOfBlockScreen(psychToolBoxConfig)

    end
    
    % Announce the end of the session
    DrawFormattedText(psychToolBoxConfig.w,'The end', 'center', 'center', psychToolBoxConfig.white);
    Screen('Flip', psychToolBoxConfig.w);
    WaitSecs(3);

    % Save the .edf file
    fprintf('\nReceiving data file ''%s''\n', edfFile);
    Eyelink('ReceiveFile', edfFile, pwd, 1);
    fprintf('Data file ''%s'' received\n', edfFile);
        
    % Save the .mat file
    save(matFile,'dat', 'description');
    
    % Clean up
    quitExperiment();
    
    fprintf("Session duration: %.2f", GetSecs - tSessionStart)
    
catch ME

    % Notify error
    disp([newline...
        '====================' newline...
        '====== Caught ======' newline...
        '====================' newline])
    
    % Save everything
    save('Caught.mat');
    
    % Clean up
    quitExperiment();
    
end

% %% Test script
% stimCfg = stimCfgsDebug{2, 2}; % generateStimuli(psychToolBoxConfig, 0.5, 5, 45);
% env = stimCfg.envSpectral;
% envV = env(:);
% [NX, NY, NF] = size(env);
% 
% df = 1;
% x = (-NX/df/2:df:NX/df/2-1) + mod(NX/df,2) * df/2;
% y = (-NY/df/2:df:NX/df/2-1) + mod(NY/df,2) * df/2;
% t = (-NF/1/2:1:NF/1/2-1) + mod(NF/1,2) * 1/2;
% 
% [fx, fy, ft] = meshgrid(x,y,t);
% 
% fxV = fx(:);
% fyV = fy(:);
% angl = atan2(fyV, fxV);
% Ae = sum( envV .* exp( 1j * 2 * angl) );
% meanAngl = 0.5 * rad2deg(angle(Ae));
% R = abs(Ae) / sum(envV);
% stdAngl = rad2deg(sqrt(-2*log(R))) / 2;
% 
% fprintf("Stim metrics: sampleAngle = %.2f, sampleStdDev = %.2f \n", meanAngl, stdAngl)
% 
% movie = randomCloud( env, false, [], false, stimCfg.counterphase );
% movie = rectify( movie, 'RMS', 0.12, false ) * (psychToolBoxConfig.white-psychToolBoxConfig.black) + psychToolBoxConfig.black;
% 
% Z = fftshift(fftn(movie));
% A = abs(Z); 
% A(round(size(A, 1)/2), round(size(A, 2)/2), round(size(A, 3)/2)) = 0; % Remove DC component
% 
% [N_X, N_Y, N_frame] = size(Z);
% [fx,fy,ft] = getGrids(N_X, N_Y, N_frame);
% 
% ori = atan2(fy, fx);
% Ae_itheta = sum( A(:) .* exp(1j* 2 * ori(:)) );
% R = abs(Ae_itheta) / sum(A(:));
% 
% meanAngle = 0.5 * rad2deg(angle(Ae_itheta));
% stdAngle = rad2deg(sqrt(-2 * log(R))) / 2;  % This is the approxiamtion for von-misses distribution
% 
% fprintf("Stim metrics: sampleAngle = %.2f, sampleStdDev = %.2f \n", meanAngle, stdAngle)
%     
% disp("done")

%% ----------------------------- FUNCTION ----------------------------- %%

%% Experiment begin screen
function beginExpScreen(psychToolBoxConfig)

% Before starting the experiment wait for user to press space key
% to begin the experiment.
DrawFormattedText(psychToolBoxConfig.w, 'Press SPACE whenever you are ready to begin experiment...', ...
    'center', 'center', 255);
Screen('Flip', psychToolBoxConfig.w);

% Wait for keypress
while true
    [~,~,keyCode] = KbCheck(-1);
    if keyCode(psychToolBoxConfig.quitKey)
        quitExperiment();
        return
    end
    
    [keyIsDown, ~, keyCode] = KbCheck;
    if keyIsDown
        if keyCode(psychToolBoxConfig.spaceKey)
            break;  % exit loop if space is pressed
        end
    end
end

KbReleaseWait; % (Optional) Clear key buffer before continuing

end

%% End of bloc                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                k screen
function endOfBlockScreen(psychToolBoxConfig)

timeoutSecs = 60;  % maximum wait time
tStart = GetSecs;

% Display end-of-block message
DrawFormattedText(psychToolBoxConfig.w, ...
    'End of block.\n\nPress SPACE to continue...\n\n(Automatically continuing in 60 seconds)', ...
    'center', 'center', 255);
Screen('Flip', psychToolBoxConfig.w);

while true
    % For emergency exit
    [~,~,keyCode] = KbCheck(-1);
    if keyCode(psychToolBoxConfig.quitKey)
        quitExperiment();
        return
    end
    
    [keyIsDown, ~, keyCode] = KbCheck;

    % Check if SPACE is pressed
    if keyIsDown && keyCode(psychToolBoxConfig.spaceKey)
        break;
    end

    % Check for timeout
    if GetSecs - tStart > timeoutSecs
        break;
    end
end

KbReleaseWait;  % Wait for all keys to be released before continuing

end

%% Initialize EyeLink
function edfFile = initializeEyeLink(subjectNum, dummymode)

% Name the .edf file
edfFile = ['COR' num2str(subjectNum,'%.2d')];

% Exit program if this fails
if ~EyelinkInit(dummymode, 1) 
    error('Eyelink Init failed.');
end

% Open file to record data
res = Eyelink('Openfile', edfFile);
if res ~= 0
    error('Cannot create EDF file.');
end

% Select which events are saved in the EDF file or available online. Include everything just in case
Eyelink('Command','file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
Eyelink('Command','link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');

% Select which sample data is saved in EDF file or available online. Include everything just in case
Eyelink('Command','file_sample_data = LEFT,RIGHT,GAZE,HREF,RAW,AREA,GAZERES,BUTTON,STATUS,INPUT');
Eyelink('Command','link_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT');
    
end

function eyeUsed = setupEyeLink(w)

% Eyelink setting
el = EyelinkInitDefaults(w);

% Calibration & validation
EyelinkDoTrackerSetup(el);

% Start recording
WaitSecs(0.05);
Eyelink('StartRecording');

% Check which eye is available
eyeUsed = Eyelink('EyeAvailable'); % 0 = left, 1 = right, 2 = binocular
% Get samples from right eye if binocular
if eyeUsed == 2
    eyeUsed = 1;
end
    
end

%% Initialize Psychtoolbox settings

function psychToolBoxConfig = initPsychToolBox()

% Some preferences
Screen('Preference','SkipSyncTests', 1);
Screen('Preference', 'VisualDebugLevel', 3);
Screen('Preference', 'SuppressAllWarnings', 1);
screenNumber = max(Screen('Screens'));
doublebuffer = 1;

% Luminance
white = 255;
black = 1;
gray = ( white - black ) / 2;

% For color correction
S = load('/Users/gorislab/Desktop/psychophysics/Calib/2018/humanRig20180612.mat');  % TODO: set according to the system - change this to a common config variable
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask','FinalFormatting','DisplayColorCorrection','SimpleGamma');

% Open window (use PsychImaging('OpenWindow',...)?)
[w, rect] = PsychImaging('OpenWindow', screenNumber, ...
    gray, [], 32, doublebuffer+1, [], 6);

% Color correction
PsychColorCorrection('SetEncodingGamma', w, S.gam.power, 'FinalFormatting');

% Get the width and height of the window in pixels
[screenWidth, screenHeight] = Screen('WindowSize', w);

% Get the center coordinate of the window
[xCenter, yCenter] = RectCenter(rect);

% Get the refresh rate of the screen
frameRate = 75; % [Hz] % This affects the dispersion computation

% Darkroom spec [cm]
monitorWidth = 36;  % TODO: set according to the system - change this to a common config variable
distance = 81;
ppd = screenWidth/(2*atand(monitorWidth/2/distance));

% Intial flip
Screen('Flip',w);

% Set up blending for drawing smooth points and lines
Screen('BlendFunction',w,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

% Set up audio
InitializePsychSound();
PsychPortAudio('Open');

% Set text size
Screen('TextSize',w,25);

% Keyboard setting 
KbName('UnifyKeyNames');
calibrateKey = KbName('c');
spaceKey = KbName('space');
quitKey = KbName('q');
zKey = KbName('z');
xKey = KbName('x');
periodKey = KbName('.>');
slashKey = KbName('/?');
leftArrowKey = KbName('LeftArrow');
rightArrowKey = KbName('RightArrow');
RestrictKeysForKbCheck([spaceKey quitKey zKey xKey periodKey slashKey leftArrowKey rightArrowKey calibrateKey]);

% Provide against interference
HideCursor;
ListenChar(2);
PL = MaxPriority(w);
Priority(PL);

psychToolBoxConfig.w = w;
psychToolBoxConfig.rect = rect;
psychToolBoxConfig.white = white;
psychToolBoxConfig.black = black;
psychToolBoxConfig.gray  = gray;
psychToolBoxConfig.screenWidth  = screenWidth;
psychToolBoxConfig.screenHeight  = screenHeight;
psychToolBoxConfig.xCenter  = xCenter;
psychToolBoxConfig.yCenter  = yCenter;
psychToolBoxConfig.frameRate  = frameRate;
psychToolBoxConfig.monitorWidth  = monitorWidth;
psychToolBoxConfig.distance  = distance;
psychToolBoxConfig.ppd  = ppd;
psychToolBoxConfig.quitKey = quitKey;
psychToolBoxConfig.spaceKey = spaceKey;
psychToolBoxConfig.calibrateKey = calibrateKey;

end

%% Quit experiment
function quitExperiment()

Eyelink('ShutDown'); 
Screen('CloseAll');
Priority(0);
ShowCursor;
ListenChar(0);
            
end


%% Initialize fixation
function fixationWinCfg = initializeFixationWin(ppd, xCenter, yCenter, black, white)

fixThreshold = .75 * ppd; % radius of the fixation window
fixDiameter = .2 * ppd; % diameter of the fixation point
fixOuterRect = CenterRectOnPointd(2*fixDiameter*[0 0 1 1], xCenter, yCenter);
fixInnerRect = CenterRectOnPointd(fixDiameter*[0 0 1 1], xCenter, yCenter);
fixRect = [fixOuterRect; fixInnerRect]';
fixColor = repmat([black white], 3, 1);

fixationWinCfg.fixThreshold = fixThreshold;
fixationWinCfg.fixDiameter = fixDiameter;
fixationWinCfg.fixRect = fixRect;
fixationWinCfg.fixColor = fixColor;

end

%% Show fixation window
function fixnData = showFixation(psychToolBoxConfig, fixationWinCfg )

Screen('FillOval', psychToolBoxConfig.w, fixationWinCfg.fixColor, fixationWinCfg.fixRect);
[~,tStart] = Screen('Flip', psychToolBoxConfig.w);
Eyelink('Message','FIXATION');

fixnData.tStart = tStart;

end

%% Stim configuration
function stimCfg = generateStimuli(psychToolBoxConfig, stimDuration, stimSpread, stimOri)

% !!!!!!!!!!
% Beware: this function is completely deterministic. No randomness!
% Changing this will have implications!
% !!!!!!!!!!

% Either generate stimuli of each combination (in factorial design) or
% generate it on the fly during each iteration - check if this has any
% implication for timing

% Stimulus location
ppd = psychToolBoxConfig.ppd;
frameRate = psychToolBoxConfig.frameRate;
white = psychToolBoxConfig.white;

% Stimulus diameter
stimDiameter = round( 3 * ppd );
stimDiameter = stimDiameter + ~mod(stimDiameter,2); % Why does stim diameter have to be even?
envSpatial = white * envelopeSpatial( stimDiameter, 'raised-cosine', 0.75 ); % Is this like aperture?

% Stimulus parameters
nFrames = ceil( stimDuration * frameRate );
sf = 1.5 * stimDiameter/ppd;                          % cpd * p/ppd = cycle
Bsf = sf * ( 2^(.5*1) - 1 );                          % 0.1; 1 % cycle * ( 2^(.5*octave) - 1 ) = cycle. Bsf / (diameter/ppd) is in cpd
Bv = 1.5 * nFrames/frameRate / (stimDiameter/ppd);    % dps * f/fps / (p/ppd) = unitless - What is this velocity? Is it bandwidth? More standard would be to define temporal frequency (Hz) or velocity in pixels/frame, not deg/s.
Vx = 0;                                               % 1 * ppd / frameRate; % dps * ppd / (1/s) = pixels/frame (velocity)
Vy = 0;                                               % 0 * ppd / frameRate;
Bth = deg2rad(stimSpread);                            % Orientation bandwidth - dispersion
theta = deg2rad(stimOri);                             % Orietation
sftScale = 1;                                         % spatiotemporal scaling factor
logGabor = true;                                      % (logical) if true it uses a log-Gabor kernel
alpha = 1;                                            % exponent for the color envelope
ft0 = 2;                                              % 2 Hz contrast modulation
Bft = 0.1;
counterphase = false;


% Spectral envelope
% Spatial and temporal frequency grids (Fourier space)
[fx,fy,ft] = getGrids( stimDiameter, stimDiameter, nFrames );
envSpectral = envelopeSpectral( fx, fy, ft, Vx, Vy, Bv, sf, Bsf, ft0, Bft, sftScale, logGabor, theta, Bth, alpha, counterphase );

stimCfg.counterphase = counterphase;
stimCfg.stimDiameter = stimDiameter;
stimCfg.envSpatial   = envSpatial;
stimCfg.envSpectral  = envSpectral;
stimCfg.nFrames      = nFrames;

end

%% Generate stimuli
function stimParam = drawStimuli(stimCfg, trlCfg, psychToolBoxConfig)
    
    movie = randomCloud( stimCfg.envSpectral, false, [], false, stimCfg.counterphase );
    movie = rectify( movie, 'RMS', trlCfg.stimContrast, false ) * (psychToolBoxConfig.white-psychToolBoxConfig.black) + psychToolBoxConfig.black;
    movietex = nan(1, stimCfg.nFrames);
    for iF = 1:length(movietex)
        movietex(iF) = Screen('MakeTexture', psychToolBoxConfig.w, cat(3,movie(:,:,iF), stimCfg.envSpatial));
    end
    
    stimX = trlCfg.stimLocX * psychToolBoxConfig.ppd;
    stimY = trlCfg.stimLocY * psychToolBoxConfig.ppd;
    movieRect = CenterRectOnPointd(stimCfg.stimDiameter*[0 0 1 1], psychToolBoxConfig.xCenter + stimX, psychToolBoxConfig.yCenter + stimY);
    
    stimParam.movie = movie;
    stimParam.movietex = movietex;
    stimParam.movieRect = movieRect;
    
end

% Calculate mean orientation and dispersion from the generated stimuli
function stimMetrics = calcStimMetrics(stimMovie)

Z = fftshift(fftn(stimMovie));
A = abs(Z); % Amplitude spectrum
% A(round(size(A, 1)/2), round(size(A, 2)/2), round(size(A, 3)/2) + 1) = 0; % Remove DC component - depends upon frameRate
A(floor(size(A, 1)/2) + 1, floor(size(A, 2)/2) + 1, floor(size(A, 3)/2) + 1) = 0; 

[N_X, N_Y, N_frame] = size(Z);
[fx,fy,ft] = getGrids(N_X, N_Y, N_frame);

ori = atan2(fy, fx);

% Compute circular mean and standard deviation
% Note: we are working with orientations which are pi-periodic. Hence,
% orientations needs to be multiplied by 2 in order to get 2pi periodicity.
Ae_itheta = sum( A(:) .* exp(1j* 2 * ori(:)) );
R = abs(Ae_itheta) / sum(A(:));

% pi periodic orientation needs to be converted back to 2pi periodic angle
meanAngle = 0.5 * rad2deg(angle(Ae_itheta));
stdAngle = rad2deg(sqrt(-2 * log(R))) / 2;  % This is the approxiamtion for von-misses distribution

stimMetrics.meanAngle = meanAngle;
stimMetrics.stdAngle = stdAngle;

% DEBUG
fprintf("Stim metrics: sampleAngle = %.2f, sampleStdDev = %.2f \n", meanAngle, stdAngle)
end

%% Check fixation 
function [fixationHeld, eyeUsed] = checkFixation(psychToolBoxConfig, eyeUsedCurr, fixThreshold, tStart, durFix, xCenter, yCenter, quitKey, checkFixationMaxDur)
fixationHeld = true;
tStartCheckFixation = tStart;
eyeUsed = eyeUsedCurr;

% Loop until the subject maintain fixation
while 1
    % For emergency exit
    [~,~,keyCode] = KbCheck(-1);
    if keyCode(quitKey)
        quitExperiment();
        return
        
    % Check fixation
    elseif Eyelink('NewFloatSampleAvailable') > 0
        evt = Eyelink('NewestFloatSample');
        gx = evt.gx(eyeUsed+1); % [left eye gaze x, right eye gaze x] +1 because this is matlab
        gy = evt.gy(eyeUsed+1);
        if ~isnan(gx) && ~isnan(gy)
            % If fixation was maintained for durFix, break while loop to show stimulus
            if norm([xCenter,yCenter]-[gx,gy]) <= fixThreshold
                if GetSecs - tStartCheckFixation >= durFix % In seconds
                    break
                end
            % Reset timer if the subject break fixation
            else
                tStartCheckFixation = GetSecs;
                
                % Eyelink sucks! Just calibrate again to avoid eyelink
                % problems.
                if keyCode(psychToolBoxConfig.calibrateKey)
                    % Eye used needs to be updated
                    eyeUsed = setupEyeLink(psychToolBoxConfig.w);
                end
                
                % If no stable fixation is achieved within 2 seconds then
                % abort
                % if (tStartCheckFixation - tStart) > checkFixationMaxDur
                %   % Show warning that some problem with fixation
                %   
                % end
            end
        end
    end
    % Avoid hogging CPU
    WaitSecs(0.001); % ~1000 Hz polling
end
            
end


%% Show response screen
function respData = showResponseScreen(psychToolBoxConfig, fixationWinCfg, trlCfg, eyeUsed, respScreenGazeHoldDur, beeperDur, respSuccessWaitDur)

arcRadi1 = 3.5;
arcRadi2 = 5.5;

% if rand > 0.5
%     arcRadi1 = 3.5;
%     arcRadi2 = 5.5;
% else
%     arcRadi1 = 5.5;
%     arcRadi2 = 3.5;
% end

% Two arcs
arcRadiusRed   = arcRadi1* psychToolBoxConfig.ppd;
arcRadiusGreen = arcRadi2 * psychToolBoxConfig.ppd;
arcTolerance   = 0.4 * psychToolBoxConfig.ppd; % +/- tolerance in pixels to match arc

% Initialize
responseGiven = false;
fixStartTime = NaN;
sampleDuration = 0.01;
reportedAngle = NaN;
reportedArc = NaN;
responseTime = NaN;

tRspWinStart = GetSecs;
% tPrint = GetSecs;

% Draw response arc
Screen('FillOval', psychToolBoxConfig.w, fixationWinCfg.fixColor, fixationWinCfg.fixRect);
arcRectRed = CenterRectOnPointd([-1 -1 1 1] * arcRadiusRed, psychToolBoxConfig.xCenter, psychToolBoxConfig.yCenter);
arcRectGreen = CenterRectOnPointd([-1 -1 1 1] * arcRadiusGreen, psychToolBoxConfig.xCenter, psychToolBoxConfig.yCenter);
Screen('FrameArc', psychToolBoxConfig.w, [255 0 0], arcRectRed, 90, -180, 4);   % red arc
Screen('FrameArc', psychToolBoxConfig.w, [0 255 0], arcRectGreen, 90, -180, 4); % green arc
[~, tStartOfRespScreen] = Screen('Flip', psychToolBoxConfig.w);

    
while ~responseGiven
    if Eyelink('NewFloatSampleAvailable') > 0
        evt = Eyelink('NewestFloatSample');
        gx = evt.gx(eyeUsed+1); % [left eye gaze x, right eye gaze x] +1 because this is matlab
        gy = evt.gy(eyeUsed+1);

        % if sample.gx(1) ~= el.MISSING_DATA && sample.gy(1) ~= el.MISSING_DATA
        if ~isnan(gx) && ~isnan(gy)

            % Compute relative to center
            dx = gx - psychToolBoxConfig.xCenter;
            dy = psychToolBoxConfig.yCenter - gy;  % screen coords: y is inverted

            % Polar coordinates
            % angleDeg = mod(atan2d(dy, dx), 180);
            angleDeg = mod(atan2d(dy, dx), 180);
            gazeRadius = sqrt(dx^2 + dy^2);

            % Determine arc (or none)
            if abs(gazeRadius - arcRadiusRed) < arcTolerance
                currentArc = "red";
                selectedRadius = arcRadiusRed;
            elseif abs(gazeRadius - arcRadiusGreen) < arcTolerance
                currentArc = "green";
                selectedRadius = arcRadiusGreen;
            else
                currentArc = "none";
            end

            % If gaze on one of the arcs
            if currentArc ~= "none"
                % Show dot
                % dotX = psychToolBoxConfig.xCenter + selectedRadius * cosd(angleDeg);
                % dotY = psychToolBoxConfig.yCenter - selectedRadius * sind(angleDeg);
                % Screen('DrawDots', psychToolBoxConfig.w, [dotX; dotY], 20, [255 0 0], [], 2);

                % fprintf("%.2f, %.2f, %.2f, %.2f", dotX, dotY, gazeRadius, selectedRadius)
                % fprintf("\n");
                
                % Confirm if fixated
                if isnan(fixStartTime)
                    fixStartTime = GetSecs;
                elseif GetSecs - fixStartTime > respScreenGazeHoldDur % record response when the gaze is stable for 0.5s
                    responseTime = GetSecs;
                    responseGiven = true;
                    reportedAngle = angleDeg;
                    reportedArc   = currentArc;

                    Eyelink('Message','RESPONSE');
                end
            else
                fixStartTime = NaN; % not on arc
            end
        end
    end
    
    % Max timeout occured - exit this
    if ~responseGiven && ( (GetSecs - tRspWinStart) > trlCfg.respMaxDur )
        responseGiven = false;
        break;
    end
    
    % Redraw arc ith current cursor
    Screen('FillOval', psychToolBoxConfig.w, fixationWinCfg.fixColor, fixationWinCfg.fixRect);
    arcRectRed = CenterRectOnPointd([-1 -1 1 1] * arcRadiusRed, psychToolBoxConfig.xCenter, psychToolBoxConfig.yCenter);
    arcRectGreen = CenterRectOnPointd([-1 -1 1 1] * arcRadiusGreen, psychToolBoxConfig.xCenter, psychToolBoxConfig.yCenter);
    Screen('FrameArc', psychToolBoxConfig.w, [255 0 0], arcRectRed, 90, -180, 4);   % red arc
    Screen('FrameArc', psychToolBoxConfig.w, [0 255 0], arcRectGreen, 90, -180, 4); % green arc
    
    if responseGiven
        % One response is given, draw the cursor on the arc and
        % make beep noise
        dotX = psychToolBoxConfig.xCenter + selectedRadius * cosd(reportedAngle);
        dotY = psychToolBoxConfig.yCenter - selectedRadius * sind(reportedAngle);
        Screen('DrawDots', psychToolBoxConfig.w, [dotX; dotY], 20, [0 0 0], [], 2);
    end
    
    Screen('Flip', psychToolBoxConfig.w);
    
    if responseGiven
        Beeper(400, .4, beeperDur);
    end
    
    WaitSecs(sampleDuration);
end

confReport = NaN;

if responseGiven
    fprintf("Reported ARC %s \n", reportedArc)
    if reportedArc == "green"
        confReport = 1;
    elseif reportedArc == "red"
        confReport = 0;
    end
end

respData.responseGiven = responseGiven;
respData.reportedAngle = reportedAngle;
respData.reportedArc   = reportedArc;
respData.reportedConf  = confReport;
respData.reactionTime  = responseTime - tStartOfRespScreen;

if responseGiven
    WaitSecs(respSuccessWaitDur - sampleDuration - beeperDur);
end

end

%% Calculate reward
function reward = calcReward(trueOri, reportedOri, confReport)

% Note: reported orientation is already pi-periodic, true orientation as
% well
maxTolerableError = 20; % In degrees
sigmaHC = sqrt(30);     % HC reward function std deviation
valLC   = 0.3;          % LC constant reward
absPerceptualError = abs(trueOri - reportedOri);
absPerceptualError = min(absPerceptualError, 180 - absPerceptualError);

if absPerceptualError > maxTolerableError
    reward = 0;
    return
end

if confReport == 1 % High confidence
    g = exp( - (absPerceptualError).^2 / (2 * sigmaHC^2) );
    reward = g;
elseif confReport == 0  % Low confidence
    reward = valLC;
else
    error('Unknown confidence report: %s. Must be "HC" or "LC".', confReport);
end

end

%% Abort trial
function abortTrial(w, durFeedbackFixBreak, beeperDur)

% Give feedback
DrawFormattedText(w, 'Poor fixation!', 'center', 'center');
Screen('Flip', w);
Beeper(400, .4, beeperDur);
WaitSecs( durFeedbackFixBreak - beeperDur );

% TODO: edit the data file
end

%% envelopeSpatial - Control aperture
function env = envelopeSpatial(N, name, par)

[X,Y] = getGrids( N, N, 1 );
switch name
    case 'Gaussian'
        sigma = par;
        env = exp(-(X.^2+Y.^2)/(2*sigma^2));
    case 'raised-cosine'
        beta = par;
        T = .5/(N/2)*(1+beta);
        env = .5*(1+cos(pi*T/beta*(abs(sqrt(X.^2+Y.^2))-(1-beta)/(2*T))));
        env( sqrt(X.^2+Y.^2) < (1-beta)/(2*T) ) = 1;
    otherwise
        error('Invalid input');
end
env( (X.^2 + Y.^2) > (N/2)^2 ) = 0; % Anything outside the spatial redius is set to zero
    
end

%% getGrids
function [fx, fy, ft] = getGrids(N_X, N_Y, N_frame)

% Even and odd
df = 1;
x = (-N_X/df/2:df:N_X/df/2-1) + mod(N_X/df,2) * df/2;
y = (-N_Y/df/2:df:N_X/df/2-1) + mod(N_Y/df,2) * df/2;
t = (-N_frame/1/2:1:N_frame/1/2-1) + mod(N_frame/1,2) * 1/2;

[fx, fy, ft] = meshgrid(x,y,t);

end

%% frequencyRadius
function fr = frequencyRadius(fx, fy, ft, sftScale, clean_division)

if sftScale == inf
    fr2 = fx.^2 + fy.^2;
else
    fr2 = fx.^2 + fy.^2 + (ft/sftScale).^2;
end

if clean_division
    fr2(fr2==0) = inf;
end

fr = sqrt(fr2);

end

%% envelopeColor
function env = envelopeColor(fx, fy, ft, alpha, sftScale)

if alpha == 0 % white noise
    env = ones(size(fx));
else
    fr = frequencyRadius(fx, fy, ft, sftScale, true);
    env = fr.^(-alpha);
end

end

%% envelopeRadial - SF envelope
function env = envelopeRadial(fx, fy, ft, sf_0, B_sf, sftScale, loggabor)

if sf_0 == 0 || B_sf == inf
    if loggabor
        env = envelopeColor(fx, fy, ft, 1, sftScale); % should I remove this?
    else
        env = ones(size(fx));
    end
elseif loggabor
    fr = frequencyRadius(fx, fy, ft, sftScale, true);
    env = exp(-.5*log(fr/sf_0).^2/log(1+B_sf/sf_0)^2);
    % env = exp(-.5*(log(fr)-log(sf_0^2/sqrt(sf_0^2+B_sf^2))).^2/log(1+B_sf^2/sf_0^2));
else
    fr = frequencyRadius(fx, fy, ft, sftScale, true);
    % env = double(abs(fr - sf_0) < 1e-2);  % allow tiny tolerance
    env = exp(-.5*(fr-sf_0).^2/B_sf^2);
end

end

%% envelopeSpeed
function env = envelopeSpeed(fx, fy, ft, V_X, V_Y, B_V, sftScale)

if size(ft,3) == 1
    env = ones(size(fx));
elseif B_V == 0
    env = zeros(size(fx));
    env(ft == 0) = 1;
else
    fr = frequencyRadius(fx, fy, ft, sftScale, true);
    env = exp(-.5*((ft + fx*V_X + fy*V_Y)).^2./(B_V*fr).^2);
    % env = double(abs((ft + fx*V_X + fy*V_Y)) < 1e-1);
end

end

%% envelopeTemporalFreq
function env = envelopeTemporalFreq(ft, ft0, Bft, counterphase)

if counterphase
    % should be dine for each x, y and t values
    env = exp(-.5*(ft - ft0).^2/Bft^2) + exp(-.5*(ft + ft0).^2/Bft^2);
end

end

%% envelopeOrientation
function env = envelopeOrientation(fx, fy, theta, B_theta)

if B_theta == inf
    env = ones(size(fx));
elseif theta == 0 && B_theta == 0
    env = zeros(size(fx));
    env(fy == 0) = 1;
else
    angl = atan2(fy, fx); % Orientation of the stimuli depends upon x and y spatial frequency component
    env = exp(cos(2*(angl-theta))/(2*B_theta)^2);
    % env = exp(-0.5*(angl-theta).^2/(2*B_theta)^2);
end

end

%% envelopeSpectral
function envelope = envelopeSpectral(fx, fy, ft, V_X, V_Y, B_V, sf_0, B_sf, ft0, Bft, sftScale, loggabor, theta, B_theta, alpha, counterphase)

if counterphase
    envelope = envelopeColor(fx, fy, ft, alpha, sftScale) ...
        .* envelopeOrientation(fx, fy, theta, B_theta) ...
        .* envelopeRadial(fx, fy, ft, sf_0, B_sf, sftScale, loggabor) ...
        .* envelopeTemporalFreq(ft, ft0, Bft, counterphase);
else
    envelope = envelopeColor(fx, fy, ft, alpha, sftScale) ...
        .* envelopeOrientation(fx, fy, theta, B_theta) ...
        .* envelopeRadial(fx, fy, ft, sf_0, B_sf, sftScale, loggabor) ...
        .* envelopeSpeed(fx, fy, ft, V_X, V_Y, B_V, sftScale);
end
    
end

%% randomCloud
function z = randomCloud(envelope, impulse, events, do_amp, symmetric)

[N_X, N_Y, N_frame] = size(envelope);
if symmetric % For counterphase stimuli we need the random phase to be symmetric
    phase = zeros(N_X, N_Y, N_frame);
    F_events = exp(1i * phase);
elseif impulse
    [fx, fy, ft] = getGrids(N_X, N_Y, N_frame);
    phase = -2*pi*(N_X/2*fx + N_Y/2*fy + N_frame/2*ft);
    F_events = exp(1i * phase);
elseif isempty(events) % This is the culprit... it assigns random phase to the points and causes asymmetry - Why the phases have to be random?
    phase = 2*pi*rand(N_X, N_Y, N_frame);
    F_events = exp(1i * phase);
    if do_amp
        F_events = F_events .* randn(N_X, N_Y, N_frame);
    end
else
    F_events = fftn( events );
    F_events = fftshift(F_events);
end
Fz = F_events .* envelope;

% De-center the spectrum
Fz = ifftshift(Fz);
Fz(1,1,1) = 0; % remove the DC component
z = real(ifftn(Fz));

end

%% rectify
function z = rectify(z, method, contrast, verbose)

if verbose
    disp('Before rectification of the frames')
    fprintf('mean = %.2g, std = %.2g, min = %.2g, max = %.2g, max(abs) = %.2g\n',...
        mean(z,'all'), std(z,[],'all'), min(z,[],'all'), max(z,[],'all'), max(abs(z),[],'all'))
end

z = z - mean(z(:));

switch method
    case 'Michelson'
        z = .5 * z/max(abs(z),[],'all') * contrast + .5;
    case 'RMS'
        z = .5 * z/std(z,[],'all') * contrast + .5;
    otherwise
        error('Invalid input.');
end

if verbose
    disp('After rectification of the frames')
    fprintf('mean = %.2f, std = %.2f, min = %.2f, max = %.2f, max(abs) = %.2f\n',...
        mean(z,'all'), std(z,[],'all'), min(z,[],'all'), max(z,[],'all'), max(abs(z),[],'all'))
    fprintf('Percentage pixels clipped = %.1f%%\n', sum(abs(z)>1,'all')*100/numel(z))
end
    
end

