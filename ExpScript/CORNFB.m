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

% TODO: set the stimulus parameter by looking into paper and then decide on
% the bounds

%% Start clean
close all;

clearvars;

% Random seed
rng shuffle

%% Use EyeLink?
dummymode = 0; % variable for eyelink usage. (1 for not using eye-tracking)

%% Initialize the script
% Experiment information

answer = str2double(inputdlg({'Subject ID:'}, 'Experiment Information')); % subject number = 0 for pilot
if sum(isnan(answer))
    error('Invalid input.')
end
subjectNum = answer(1);
% sessionNum = 1;          % Default value of session number
% newSession = true;

% Name the .mat file
matFile = ['CORNFB' num2str(subjectNum,'%.2d') '.mat'];
% Mat file for train block
tMatFile = ['CORNFB' num2str(subjectNum,'%.2d') '_train.mat'];

% dataFilePath = '/Users/gorislab/Desktop/psychophysics/experiments/COR/Data';
% addpath(dataFilePath)

% 1. If file exist 
%   - Get last run session i.e. 2 
%   - if there are remaining trials in last run session then ask user and run that
%     session again for the remaining trials
%   - if all the trials are completed then ask the user if they want to run
%     new session for this subject
% 2. If file does not exist then just take the subject number and initialize
%    session to 1
if ~isempty(dir(matFile)) % File exist
    % Load table from the previous session
    load(matFile)

    lastRunSessionNum = max(dat.session(:));
    remainingTrilsInLastRun = numel( find( (dat.session == lastRunSessionNum) & (dat.trialStatus == 0) ) );

    if isnan(lastRunSessionNum) && lastRunSessionNum > 0
        error('Something weird! Last run session num is nan.')
    end
    
    % There are still some remaining trials in last run. Ask user if they
    % want to continue with the last session or create a new one
    if remainingTrilsInLastRun > 0
        choice = questdlg(sprintf(['%d trials still needs to be completed from the last session! Do you want to' ...
            ' continue with the previous session or create a new one?'], remainingTrilsInLastRun), ...
            'Confirm', 'Continue', 'New', 'Continue');

        if strcmp(choice, 'Continue')
            sessionNum = lastRunSessionNum; % Continue with previous session
            newSession = false;
        elseif strcmp(choice, 'New')
            sessionNum = lastRunSessionNum + 1; % Create new session
            newSession = true;
        else
            error('Make a selection!')
        end
        
    else % No remaining trials. Continue with new sessions for this user
        choice = questdlg(['A session already exist for this subject. ' ...
            'Do you want to continue with a new session?'], ...
            'Confirm', 'Yes', 'No', 'No');

        if strcmp(choice, 'Yes')
            sessionNum = lastRunSessionNum + 1; % Start new session
            newSession = true;
        elseif strcmp(choice, 'No')
            error('Nothing to do here! Exiting!')
        else
            error('Make a selection!')
        end
    end

else % File does not exist
    sessionNum = 1; % Running new session
    newSession = true;
end

% Add instruction screen

% TODO: before every run, check
% 1. orientation
% 2. numblocks
% 3. trials per block
% 4. stim params
% 5. Reward params
% 6. Fixation tolerance
% 7. Eyelink setting - centroid, illumination, ( pupil: ~90, cornea: ~220), check right eye tracking
% 8. Response screen - switch arcs

%% Table
% Extreme values
% Contrast 0.015, duration 0.06, loc 4.5, spread 5
% Contrast 0.018, duration 0.06, loc 4.5, spread 45 
interTrlInterval = 1;                                                      % Inter trial interval in seconds
fixationDur = 0.5;                                                         % Fixation duration in seconds
stimOrientations = 0:15:179; %linspace(0, 179, 10);                                   % 
stimLoc_x = 0;                                                             % Stimulus location in visual field degrees
stimLoc_y = 5;                                                             % Stimulus location in visual field degrees
stimDur = [0.05, 0.3];  %0.05 0.15 0.5. 0.15 is  the minimum               % Stimulus duration in seconds
stimSpread = [5, 25]; % 5, 25, 45                                          % Stimulus spread in degrees
% stimContrast = [0.015, 0.05];                                            % Stimulus contrast levels
stimContrast = [0.015, 0.05];                                               % 0.022 - seems good contrast level to get 50:50 conf report (HC:Lc)
respMaxDur = 5;                                                            % 0.010 Maximum allowed time for user to respond (2 seconds)
respSuccessWaitDur = 0.5;
numBlocks   = 7; %6;                                                       % Number of blocks 
% nTrialsPerBlock = numel(stimOrientations)*numel(stimSpread)*numel(stimContrast)*numel(stimDur);    % Assuming each trial takes max of 5 second, a block should take ~8 minutes
% nTrialsPerBlock = numel(stimOrientations)*8;                             % TODO: delete
% nTrialsPerBlock = numel(stimOrientations)*2; % Train
nTrialsPerBlock = numel(stimOrientations)*6; % Actual ep, 6 corresponds to no of conditions and not the block

nTrials = numBlocks*nTrialsPerBlock;                                       % Total number of trials to run in this session
durFeedback = 1;
beeperDur = 0.05;
durFeedbackFixBreak = .5;
checkFixationMaxDur = 2;
respScreenGazeHoldDur = 0.2;
stimToRespScreenTransitionDelay = 0.5;
trlStatus = 0;                                                             % 0 for not completed yet, 1 for successful completion

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
% combinations = [d(:), s(:), c(:)];  % 8x3

% combinations = [combinations(2, :); combinations(4:end, :)]; %TODO: delete
% combinations = [
%     0.3000    10.0000    0.0200
%     0.3000    10.0000    0.0800
%     0.3000    30.0000    0.0800
%     0.3000    30.0000    0.0200
%     
%     0.1000    10.0000    0.0200 % can go up to 20ms, 40ms
%     0.1000    30.0000    0.0200
% ];

% combinations = [
%     % 0.0800    20.0000    0.0220
%     % 0.0800    30.0000    0.0250
%     % 0.3000    20.0000    0.0500  %%%%% Easiest
%     % 0.3000    35.0000    0.0500
%     0.0800    35.0000    0.0250    %%%%% Hardest
% ];

% combinations = [
%     % 0.0800    20.0000    0.0220
%     % 0.0800    30.0000    0.0250
%     %0.3000    20.0000    0.0500  %%%%% Easiest
%     %0.3000    30.0000    0.0800  %%%%% Easiest
%     0.0800    30.0000    0.0220
%     % 0.3000    35.0000    0.0500
% %     0.0800    35.0000    0.0250    %%%%% Hardest
% ];

% Yichao
% combinations = [
%     0.3000    10.0000     0.0600
%     0.3000    30.0000     0.0600
%     0.3000    10.0000     0.0220
%     0.3000    30.0000     0.0220
%     
%     0.0800    10.0000     0.0220
%     0.0800    30.0000     0.0220
% ];

% Jonathan
% Dispersion 30 or 35?
% combinations = [
%     0.3000    10.0000     0.0600
%     0.3000    30.0000     0.0600
% ];

% 0.3000    30.0000     0.0600


% Jonathan
combinations = [
    0.3000    10.0000     0.0600
    0.3000    30.0000     0.0600
    0.3000    10.0000     0.0220
    0.3000    30.0000     0.0220
    
    0.0800    10.0000     0.0220
    0.0800    30.0000     0.0220
];

% % Reduce contrast (0.018) and increase duration (0.1200) - for the next
% % subject
% % David
% % Stim value changed
% combinations = [
%     0.3000    10.0000     0.0400
%     0.3000    30.0000     0.0400
%     0.3000    10.0000     0.0180
%     0.3000    30.0000     0.0180
%     
%     0.1000    10.0000     0.0180
%     0.1000    30.0000     0.0180
% ];


    
%  % 0.0800    20.0000    0.0220
%     % 0.0800    30.0000    0.0250
%     %0.3000    20.0000    0.0500  %%%%% Easiest
%     %0.3000    30.0000    0.0800  %%%%% Easiest
%     0.0800    30.0000    0.0220
%     % 0.3000    35.0000    0.0500
% %     0.0800    35.0000    0.0250    %%%%% Hardest

% Dur: 0.3, 0.08 spread: 30, contrast: 0.022
% 

% Dispersion
% 0.3000    10.0000    0.0800
% 0.3000    30.0000    0.0800
% 0.3000    45.0000    0.0800

% Contrast and time
% 0.3000    10.0000    0.0250
% 0.3000    30.0000    0.0250
% 0.3000    45.0000    0.0250

% Contrast and time
% 0.0800    10.0000    0.0250
% 0.0800    30.0000    0.0250
% 0.0800    45.0000    0.0250

% Workflow - find easiest and  hardest and then set rest of the parameters

% % Jiaming & Tien
% rewardConfig.maxTolerableError = 36; %36;
% rewardConfig.y1                = 18; %18;
% rewardConfig.x1                = 16.1; %16.1;
% rewardConfig.c1                = 0.2;
% rewardConfig.g1                = 0.3; % 0.3
% rewardConfig.g2                = 0.3; % 0.3

rewardConfig.maxTolerableError = 30; %36;
rewardConfig.y1                = 17; %18;
rewardConfig.x1                = 15; %16;
rewardConfig.c1                = 0.15;
rewardConfig.g1                = 1.2; % 0.3
rewardConfig.g2                = 1.5; % 0.3

% Values set accrding to the values obtained from calibration
fixationScreenConfig.fixationWinTolerance = 0.75; % Radius: Based on the mean err in degree obtained from calibration
responseScreenConfig.responseArcTolerance = 1;    % Radius: Based on the mean err in degree obtained from calibration = diameter = 2*maxErr
responseScreenConfig.arcRadius1           = 4.5;
responseScreenConfig.arcRadius2           = 4.5 + 2*responseScreenConfig.responseArcTolerance + 0.3; % set this according to the calibrtion error
responseScreenConfig.switchRespAcs        = true; %true

comboRepeated = repmat(combinations, length(stimOrientations), 1);  % Repeat the combinations for each orientation (10 times)
orientationsRepeated = repelem(stimOrientations(:), size(combinations, 1));  % Repeat orientations to match combination matrix
stimMatrix = [orientationsRepeated(:), comboRepeated]; % Final matrix: each row = [orientation, duration, spread, contrast]
finalStimMatrix = repmat(stimMatrix, numBlocks, 1);

% Shuffle the stim matrix across all the trials (note: shuffling is not within a single block)
shuffledIdx = randperm(size(finalStimMatrix, 1));
finalStimMatrix = finalStimMatrix(shuffledIdx, :);

% stimOriVec        = rand(nTrials, 1)*179; % Replace fixed orientations with random orientations sampled from uniform distribution (0 - 179)
stimOriVec        = finalStimMatrix(:, 1); % Stim orientation vector for all the trials  in this block
stimDurVec        = finalStimMatrix(:, 2); % Stim duration vector for all the trials  in this block
stimDispersionVec = finalStimMatrix(:, 3); % Stim dispersion vector for all the trials  in this block
stimContrastVec   = finalStimMatrix(:, 4); % Stim contrast dispersion vector for all the trials  in this block

stimSampleMeanOriVec   = nan(nTrials,1);
stimSampleSpreadVec    = nan(nTrials,1);
reportedOriVec         = nan(nTrials,1);
reportedConfVec        = nan(nTrials,1);
rewardVec              = nan(nTrials,1);
respTimeVec            = nan(nTrials,1);
rawOriErrorVec         = nan(nTrials,1);
trlErrorVec            = nan(nTrials,1);
trlAbortedVec          = nan(nTrials,1);

if ~exist('dat', 'var') % Data file does not exist, most likely first session
    % Finalize the table
    variableNames = {...
        'subject','session','block','trial', 'ITI',... % TODO: experiment information
        'fixationDur', 'stimLocX', 'stimLocY', 'stimDur', 'stimOri', 'stimSpread', 'stimContrast',... % TODO: independent and confounding variables
        'stimSampleMeanOri', 'stimSampleSpread', ...
        'respMaxDur', 'reportedOri','reportedConf','reward','reactionTime', ... % TODO: dependent variables
        'trialStatus', 'rawOriError', 'trlError', 'trlAborted'
        };
    
    dat = table( ...
            subjectVec, sessionVec, blockVec, trialVec, itiVec, ...
            fixDurVec, stimLocVec_x, stimLocVec_y, stimDurVec, stimOriVec, stimDispersionVec, stimContrastVec, ...
            stimSampleMeanOriVec, stimSampleSpreadVec, ...
            respMaxDurVec, reportedOriVec, reportedConfVec, rewardVec, respTimeVec, ...
            trlStatusVec, rawOriErrorVec, trlErrorVec, trlAbortedVec, ...
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
        'Raworientation error - Wrapped to -90 to 90 degrees'; ...
        'Trial Error - subject performance {1: Error (negative reward), 0: otherwise}'; ...
        'Trial aborted - set it ot 1 if trial was aborted (note these are tried again and eventually completed)'; ...
        }),...
        'VariableNames',{'description'},'RowNames',variableNames');

elseif newSession % dat file already exists and this is a new session
    newData = table( ...
        subjectVec, sessionVec, blockVec, trialVec, itiVec, ...
        fixDurVec, stimLocVec_x, stimLocVec_y, stimDurVec, stimOriVec, stimDispersionVec, stimContrastVec, ...
        stimSampleMeanOriVec, stimSampleSpreadVec, ...
        respMaxDurVec, reportedOriVec, reportedConfVec, rewardVec, respTimeVec, ...
        trlStatusVec, rawOriErrorVec, trlErrorVec, trlAbortedVec, ...
        'VariableNames', dat.Properties.VariableNames ...
    );
    
    dat = [dat; newData];  % Append new rows to the existing table
end

trlBlockData = dat(randperm(numel(dat(:, {'subject'})), 30), :);

%% Reward
baselineReward           = 0; % If done over multiple session, then this should probably be reloaded from the previous session
currentTotalReward       = 0; % This will not be baseline if a session file for this subject already exist
currentTotalRewardPoints = 0; %convertRewardValToPoints(baselineReward);

if exist('metaData', 'var')
    currentTotalReward = metaData.currentTotalReward;
    currentTotalRewardPoints = metaData.currentTotalRewardPoints;
end

% Other stats
abortTrlCnt = 0;
%%
tSessionStart = GetSecs;
try
    %% Initialize EyeLink
    edfFile = initializeEyeLink(subjectNum, dummymode);
    
    % Start timer
    tStartSession = GetSecs;
    
    % Initialize Psychtoolbox
    psychToolBoxConfig = initPsychToolBox();   
    
    % Train block goes here
    showTrainBlock(trlBlockData, psychToolBoxConfig, tMatFile, fixationScreenConfig, responseScreenConfig)
    % Maybe showing feedback is not such a good idea.
    
    
    %% Main experiment block goes here
    % Show trial and instruction before the experiment starts
    % Spectral env should be generated before
    % Before the experiment block starts - do following
    % Obtain fixation window (No need to do this in every trial). Just
    % obtain once and be done with it
    fixationWinCfg = initializeFixationWin(psychToolBoxConfig.ppd, ...
        psychToolBoxConfig.xCenter, psychToolBoxConfig.yCenter, ...
        psychToolBoxConfig.black, psychToolBoxConfig.white, fixationScreenConfig);

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
    % TODO: don't go to completed block. Just start with block yet to
    % be completed.
    completedTrialsData = dat((dat.session == sessionNum) & (dat.trialStatus == 1), :);
    if numel(completedTrialsData) > 0
        completedBlockIdxes = unique(completedTrialsData.block);
        currBlockIdx = max(completedBlockIdxes);
        
        thisBlockData = find( (dat.session == sessionNum) & ...
            (dat.block == currBlockIdx) & (dat.trialStatus == 1) );
        nCompletedTrialsCurrBlock = numel( thisBlockData ); 
        
        if nCompletedTrialsCurrBlock >= nTrialsPerBlock
            startBlockIdx = currBlockIdx + 1;
        else
            startBlockIdx = currBlockIdx;
        end
    else
        startBlockIdx = 1;
    end
    
    % Get total count of trials completed till last session
    if (sessionNum > 1)
        trialsDataPrevSession = find((dat.session <= (sessionNum - 1)));
        totalTrialsTillPrevSession = numel(trialsDataPrevSession);
    else
        totalTrialsTillPrevSession = 0;
    end
    
    fprintf("Session %d, block %d, Trials till prev session %d", ...
        sessionNum, startBlockIdx, totalTrialsTillPrevSession);
    fprintf("\n")
    
    % Run trials
    for blockIDx=startBlockIdx:numBlocks
        
        tStartBlock = GetSecs;

        % Calibrate eyelink at start of each block
        eyeUsed = setupEyeLink(psychToolBoxConfig.w);
        
        % Wait for user to start the experiment
        beginExpScreen(psychToolBoxConfig)
        Eyelink('Message','START');

        % Get total number of completed trials in current block and session
        thisBlockData = find( (dat.session == sessionNum) & ...
            (dat.block == blockIDx) & (dat.trialStatus == 1) );
        nCompletedTrialsCurrBlock = numel( thisBlockData ); 
        fprintf("Completed trial in this block %d / %d \n", nCompletedTrialsCurrBlock, nTrialsPerBlock)
        
        % Set current cursor at first not completed trial
        currTrlCursor = find( (dat.session == sessionNum) & ...
            (dat.block == blockIDx) & (dat.trialStatus == 0) , 1 );
        
        % Dummy value - for ITI calculations
        tEndPrevTrl = GetSecs;
        
        % Print trial info at the end of each trial
        colsToPrint = {'block', 'trial', 'stimOri', 'stimDur', 'stimSpread', 'stimContrast', 'reportedOri', 'reportedConf', 'rawOriError', 'trlError'};
        % disp(dat(1, colsToPrint)) % Print dummy row
        fprintf("Block \t Trial \t stimOri \t stimDur \t stimSpread \t stimContrast \t reportedOri \t reportedConf \t rawOriError \t trlError")
        fprintf("\n")
        
        while nCompletedTrialsCurrBlock < nTrialsPerBlock
        % while nCompletedTrialsCurrBlock < 10
            
            % === Fetch the trial configuration from the dat file ===
            trlCfgIdx = currTrlCursor;
            trlCfg = dat(trlCfgIdx, :);
            
            % Just do a sanity check to make sure this is a correct row.
            % Ideally this should always be false (i.e. no error)
            if (trlCfg.session ~= sessionNum) || (trlCfg.block ~= blockIDx)
                fprintf("Session: %d, %d; Block IDx: %d, %d", ...
                    trlCfg.session, sessionNum, trlCfg.block, blockIDx);
                error("Mismatched session num or block idx!")
            end
            
            % % === Fetch the trial configuration from the dat file ===
            % trlCfgIdx = find( (dat.session == sessionNum) & (dat.block == blockIDx) & (dat.trial == trialIDx) );
            % trlCfg = dat( (dat.session == sessionNum) & (dat.block == blockIDx) & (dat.trial == trialIDx), :);
            
            % Just to confirm, make sure that trial status is not completed
            if trlCfg.trialStatus == 1
                error("!!!! Trial already done!!!! " + ...
                    "Ideally the control should never go here! Bug somewhere.")
            end
            
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
            % fprintf("Actual Stim: angle = %.2f, stdev = %.2f \n", trlCfg.stimOri, trlCfg.stimSpread);
            
            % Inter trial interval - Do whateven preperation needs to be
            % done before this ITI
            WaitSecs( trlCfg.ITI - ( GetSecs - tEndPrevTrl ) );
            
            % === Fixation point ====
            fixnData = showFixation(psychToolBoxConfig, fixationWinCfg);
            
            % Wait indefinately untill subject maintains fixation
            [fixationHeld, eyeUsed] = checkFixation(eyeUsed, fixationWinCfg.fixThreshold, ...
                fixnData.tStart, trlCfg.fixationDur, ...
                psychToolBoxConfig.xCenter, psychToolBoxConfig.yCenter, psychToolBoxConfig.quitKey);
            
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
                % Append the current trial to the end of the block. Swap 
                % the current trial at the "end of the block" for current
                % session only
                dat.trlAborted(trlCfgIdx) = 1; % Set abort status to 1
                blockIdxes = find( (dat.session == sessionNum) & (dat.block == blockIDx) );
                blockEndIdx = blockIdxes(end);
                currTrlIdx = trlCfgIdx;
                dat = [ dat( [1:currTrlIdx - 1, currTrlIdx + 1:blockEndIdx], :); 
                    dat(currTrlIdx, :); dat(blockEndIdx+1:end, :) ]; % Test edge cases
                
                % Sanity check
                % TODO: fix hack
                if numel(dat(:,1)) > totalTrialsTillPrevSession + numBlocks*nTrialsPerBlock
                    error("Block length greater than what it should be!")
                end
                
                abortTrlCnt = abortTrlCnt + 1;
                fprintf("Trial aborted (Poor fixation)")
                fprintf('\n')

                abortTrial(psychToolBoxConfig.w, durFeedbackFixBreak, beeperDur);
                % Everytime trial is aborted, update the timer for ITI calculation
                [~, tEndPrevTrl] = Screen('Flip', psychToolBoxConfig.w);
                continue;
            end
            
            % Wait for some time before showing response scrren
            Screen('FillOval', psychToolBoxConfig.w, fixationWinCfg.fixColor, fixationWinCfg.fixRect);
            Screen('Flip', psychToolBoxConfig.w);
            WaitSecs(stimToRespScreenTransitionDelay); % Wait 0.5 secs
            
            % === Response screen ====
            % Screen('Flip', psychToolBoxConfig.w);
            respData = showResponseScreen(psychToolBoxConfig, fixationWinCfg, ...
                trlCfg, eyeUsed, respScreenGazeHoldDur, beeperDur, respSuccessWaitDur, responseScreenConfig);
            
            if ~respData.responseGiven
                % Append the current trial to the end of the block. Swap 
                % the current trial at the "end of the block" for current
                % session only
                dat.trlAborted(trlCfgIdx) = 1; % Set abort status to 1
                blockIdxes = find( (dat.session == sessionNum) & (dat.block == blockIDx) );
                blockEndIdx = blockIdxes(end);
                currTrlIdx = trlCfgIdx;
                dat = [ dat( [1:currTrlIdx - 1, currTrlIdx + 1:blockEndIdx], :); 
                    dat(currTrlIdx, :); dat(blockEndIdx+1:end, :) ];
                
                % Sanity check
                % TODO: fix this - this isnot right - no of blocks might
                % var for different session
                if numel(dat(:,1)) > totalTrialsTillPrevSession + numBlocks*nTrialsPerBlock
                    error("Block length greater than what it should be!")
                end
                
                % Show error text
                DrawFormattedText(psychToolBoxConfig.w, 'Timeout occured while waiting for response', ...
                    'center', 'center',[255 0 0]);
                Screen('Flip', psychToolBoxConfig.w);
                
                abortTrlCnt = abortTrlCnt + 1;
                fprintf("Trial aborted (Response screen timeout)")
                fprintf('\n')
                
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
            reward = calcReward(trlCfg.stimOri, reportedOri, confReport, rewardConfig);
            % currentTotalReward = currentTotalReward + reward;
            
            % Conver reward value to points - 100 points = 2 dollar (but maybe make it 1 dollar)
            rewardPoints = convertRewardValToPoints(reward);
            currentTotalRewardPoints = currentTotalRewardPoints + rewardPoints;
            currentTotalReward = baselineReward + convertPointsToUSD(currentTotalRewardPoints);
            
            % Show beeper and wait for the response screen
            Beeper(800, .4, beeperDur);
            
            % WaitSecs(respSuccessWaitDur - beeperDur);
            WaitSecs(respSuccessWaitDur);
            
            % Start timer for inter-trial interval
            [~, tEndPrevTrl] = Screen('Flip', psychToolBoxConfig.w);
            
            % Record everything - Like reaction time, choice and other data
            % in dat file.
            
            rawError = reportedOri - trlCfg.stimOri;
            rawOriError = mod(rawError + 90, 180) - 90;
        
            % Save info in file
            dat.stimSampleMeanOri(trlCfgIdx)  = trlCfg.stimOri + (stimMetrics.meanAngle - oriRefStim);
            dat.stimSampleSpread(trlCfgIdx)   = stimMetrics.stdAngle;
            dat.reportedOri(trlCfgIdx)        = reportedOri;
            dat.reportedConf(trlCfgIdx)       = confReport;
            dat.reward(trlCfgIdx)             = reward;
            dat.reactionTime(trlCfgIdx)       = respData.reactionTime;
            dat.trialStatus(trlCfgIdx)        = 1; % This trial was completed
            dat.rawOriError(trlCfgIdx)        = rawOriError;
            dat.trlError(trlCfgIdx)           = reward < 0; % Error if negative reward
            
            % Update number of completed trials in this block
            nCompletedTrialsCurrBlock = nCompletedTrialsCurrBlock + 1;
            currTrlCursor = currTrlCursor + 1; % Move to next to be completed trial
            
            % Print End of trial stats
            fprintf("%02d \t %03d \t %06.2f \t %.4f \t %04.2f \t\t %05.4f \t %06.2f \t\t %1d \t %06.2f \t %1d", ...
                dat(trlCfgIdx, colsToPrint).block, ...
                dat(trlCfgIdx, colsToPrint).trial, ...
                round( dat(trlCfgIdx, colsToPrint).stimOri, 2), ...
                round( dat(trlCfgIdx, colsToPrint).stimDur, 4), ...
                round( dat(trlCfgIdx, colsToPrint).stimSpread, 2), ...
                round( dat(trlCfgIdx, colsToPrint).stimContrast, 4), ...
                round( dat(trlCfgIdx, colsToPrint).reportedOri, 2), ...
                dat(trlCfgIdx, colsToPrint).reportedConf, ...
                round( dat(trlCfgIdx, colsToPrint).rawOriError, 2), ...
                dat(trlCfgIdx, colsToPrint).trlError);
            fprintf("\n")
        
        end
        
        
        % End of block related operations
        Eyelink('Message', 'END'); % Send event to eyelink
        Eyelink('StopRecording'); % Stop recording
        
        % Close the .edf file
        WaitSecs(1);
        Eyelink('CloseFile');

        
        % =====================================================================
        % Print stats 
        % =====================================================================
        thisBlockData = dat( (dat.session == sessionNum) & (dat.block == blockIDx), :);
        
        % Calculate maximum reward subject could have earned in this block
        trueOriVec_        = thisBlockData.stimOri;
        reportedOriVec_    = thisBlockData.reportedOri;
        reportedConfigVec_ = thisBlockData.reportedConf;
        
        thisBlockOptRewardData = calcMaxEarnedReward(trueOriVec_, ...
        reportedOriVec_, reportedConfigVec_, rewardConfig);
        
        summaryTable = groupsummary(thisBlockData, {'reportedConf'}, ...
                            {@(x) sqrt(nanmean(x.^2)), 'mean', 'std', 'numel'}, 'rawOriError');

        blockDuration = GetSecs - tStartBlock;
        minutes = floor(blockDuration / 60);
        seconds = mod(blockDuration, 60);
         
        fprintf('\n\n')
        disp(summaryTable)
        
        fprintf('\n\n')                        
        fprintf("Block duration: %d min %.2f s\n", minutes, seconds);
        fprintf("Completed trials %d/%d \n", nCompletedTrialsCurrBlock, nTrialsPerBlock);
        fprintf('Current total reward: $%.2f \n', currentTotalReward);
        fprintf('Percent Wrong (negative reward): %.2f \n', 100*sum(thisBlockData.trlError)/numel(thisBlockData.trlError)); % this should be less than 50%
        fprintf('Aborted trials: %d/%d', sum(thisBlockData.trlAborted, 'omitnan'), numel(thisBlockData.trlAborted))
        fprintf('\n\n')
        
        metaData.currentTotalReward         = currentTotalReward;
        metaData.currentTotalRewardPoints   = currentTotalRewardPoints;
        
        save(matFile,'dat', 'description', 'metaData'); % Save info
        endOfBlockScreen(psychToolBoxConfig, currentTotalReward, thisBlockOptRewardData) % Show waiting screen
        
    end
    
    % =====================================================================
    % Print stats - probably don't print session stats!
    % A session is not necessarily run continuously in this design.
    % =====================================================================
    thisSessionData = dat( (dat.session == sessionNum), :);
    summaryTable = groupsummary(thisSessionData, {'reportedConf'}, ...
                            {@(x) sqrt(nanmean(x.^2)), 'mean', 'std', 'numel'}, 'rawOriError');

    fprintf('\n\n')
    disp(summaryTable)
        
    sessionDuration = GetSecs - tSessionStart;
    minutes = floor(sessionDuration / 60);
    seconds = mod(sessionDuration, 60);
        
    fprintf('\n\n')                        
    fprintf("Session duration: %d min %.2f s\n", minutes, seconds);
    fprintf("Completed trials %d/%d \n", sum(thisSessionData.trialStatus), numel(thisSessionData.trialStatus));
    fprintf('Total earned reward: $%.2f \n', currentTotalReward);
    fprintf('Percent Wrong (negative reward): %.2f\n', 100*sum(thisSessionData.trlError)/numel(thisSessionData.trlError)); % this should be less than 50%
    fprintf('Aborted trials: %d/%d', sum(thisSessionData.trlAborted, 'omitnan'), numel(thisSessionData.trlAborted))
    fprintf('\n\n')
    
    % Announce the end of the session
    DrawFormattedText(psychToolBoxConfig.w,'The end', 'center', 'center', psychToolBoxConfig.black);
    Screen('Flip', psychToolBoxConfig.w);
    WaitSecs(3);
    
    % Save the .edf file
    fprintf('\nReceiving data file ''%s''\n', edfFile);
    Eyelink('ReceiveFile', edfFile, pwd, 1);
    fprintf('Data file ''%s'' received\n', edfFile);
        
    % Save the .mat file
    save(matFile,'dat', 'description', 'metaData');
    
    % Clean up
    quitExperiment();
    
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

%% ----------------------------- FUNCTION ----------------------------- %%

%% Show trial block
function showTrainBlock(trlBlockData, psychToolBoxConfig, tMatFile, fixationScreenConfig, responseScreenConfig)

respSuccessWaitDur = 0.5;
durFeedback = 1;
beeperDur = 0.05;
durFeedbackFixBreak = .5;
respScreenGazeHoldDur = 0.2; 

% --- ASK USER WHETHER TO RUN BLOCK ---
DrawFormattedText(psychToolBoxConfig.w, ...
    'Hi! \n\nBefore starting, do you want to run a train block to get familiar with the experiment?\n\nPress Y for YES, N for NO.', ...
    'center', 'center', 0);
Screen('Flip', psychToolBoxConfig.w);

% Wait for a key press
runBlock = false; % default
while true
    [keyIsDown, ~, keyCode] = KbCheck;
    if keyIsDown
        if keyCode(KbName('Y'))  % yes
            runBlock = true;
            break;
        elseif keyCode(KbName('N'))  % no
            runBlock = false;
            break;
        end
    end
end

% Clear screen after response
Screen('Flip', psychToolBoxConfig.w);

if runBlock
       
    fixationWinCfg = initializeFixationWin(psychToolBoxConfig.ppd, ...
            psychToolBoxConfig.xCenter, psychToolBoxConfig.yCenter, ...
            psychToolBoxConfig.black, psychToolBoxConfig.white, fixationScreenConfig);
    
    oriRefStim = 45; 
    
    % Calibrate eyelink at start of each block
    eyeUsed = setupEyeLink(psychToolBoxConfig.w);
    
    % Wait for user to start the experiment
    beginExpScreen(psychToolBoxConfig)
    Eyelink('Message','TRIAL BLOCK START');
    
    nTrialInTrialBlock = numel( trlBlockData{:, 1} );
    currTrlCursor = 1;
    nCompletedTrialsCurrBlock = 0;
    tEndPrevTrl = GetSecs;
    
    % Print trial info at the end of each trial
    colsToPrint = {'block', 'trial', 'stimOri', 'stimDur', 'stimSpread', 'stimContrast', 'reportedOri', 'reportedConf', 'rawOriError', 'trlError'};
    fprintf("Block \t Trial \t stimOri \t stimDur \t stimSpread \t stimContrast \t reportedOri \t reportedConf \t rawOriError \t trlError")
    fprintf("\n")
    
    while nCompletedTrialsCurrBlock < nTrialInTrialBlock
    
        trlCfgIdx = currTrlCursor;
        trlCfg = trlBlockData(trlCfgIdx, :);
        
        stimCfg = generateStimuli(psychToolBoxConfig, trlCfg.stimDur, trlCfg.stimSpread, oriRefStim);      
        stimParam = drawStimuli(stimCfg, trlCfg, psychToolBoxConfig);
        stimMetrics = calcStimMetrics(stimParam.movie);
        
        WaitSecs( trlCfg.ITI - ( GetSecs - tEndPrevTrl ) );
    
        % === Fixation point ====
        fixnData = showFixation(psychToolBoxConfig, fixationWinCfg);
    
        % Wait indefinately untill subject maintains fixation
        [~, eyeUsed] = checkFixation(eyeUsed, fixationWinCfg.fixThreshold, ...
            fixnData.tStart, trlCfg.fixationDur, ...
            psychToolBoxConfig.xCenter, psychToolBoxConfig.yCenter, psychToolBoxConfig.quitKey );
                
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
            
            nCompletedTrialsCurrBlock = nCompletedTrialsCurrBlock + 1;
            currTrlCursor = currTrlCursor + 1; % Move to next to be completed trial

            abortTrial(psychToolBoxConfig.w, durFeedbackFixBreak, beeperDur);
            % Everytime trial is aborted, update the timer for ITI calculation
            [~, tEndPrevTrl] = Screen('Flip', psychToolBoxConfig.w);
            continue;
        end
    
        % === Response screen ====
        respData = showResponseScreen(psychToolBoxConfig, fixationWinCfg, ...
            trlCfg, eyeUsed, respScreenGazeHoldDur, beeperDur, respSuccessWaitDur, responseScreenConfig);
                
        if ~respData.responseGiven
            
            nCompletedTrialsCurrBlock = nCompletedTrialsCurrBlock + 1;
            currTrlCursor = currTrlCursor + 1; % Move to next to be completed trial

            % All aboted trial will be tried again towards the end
            abortTrial(psychToolBoxConfig.w, durFeedbackFixBreak, beeperDur);
            % Everytime trial is aborted, update the timer for ITI calculation
            [~, tEndPrevTrl] = Screen('Flip', psychToolBoxConfig.w);
            continue;
        end
        
        stimOri       = trlCfg.stimOri;
        reportedOri   = respData.reportedAngle;
        confReport    = respData.reportedConf;
          
        % Show beeper and wait for the response screen
        Beeper(800, .4, beeperDur);
        
        WaitSecs(respSuccessWaitDur);
        
        Screen('Flip', psychToolBoxConfig.w);
        
        % WaitSecs(0.5);
        % Show feedback screen
        % showTrainTrialFeedbackScreen(psychToolBoxConfig, fixationWinCfg, stimOri, reportedOri, confReport);
        
        % WaitSecs(durFeedback);
        
        % Start timer for inter-trial interval
        [~, tEndPrevTrl] = Screen('Flip', psychToolBoxConfig.w);     
        
        rawError = reportedOri - trlCfg.stimOri;
        rawOriError = mod(rawError + 90, 180) - 90;
            
        % Save info in file
        trlBlockData.stimSampleMeanOri(trlCfgIdx)  = trlCfg.stimOri + (stimMetrics.meanAngle - oriRefStim);
        trlBlockData.stimSampleSpread(trlCfgIdx)   = stimMetrics.stdAngle;
        trlBlockData.reportedOri(trlCfgIdx)        = reportedOri;
        trlBlockData.reportedConf(trlCfgIdx)       = confReport;
        trlBlockData.reward(trlCfgIdx)             = nan;
        trlBlockData.reactionTime(trlCfgIdx)       = respData.reactionTime;
        trlBlockData.trialStatus(trlCfgIdx)        = 1; % This trial was completed
        trlBlockData.rawOriError(trlCfgIdx)        = rawOriError;
        trlBlockData.trlError(trlCfgIdx)           = nan; % Error if negative reward
            
        % Print End of trial stats
        fprintf("%02d \t %03d \t %06.2f \t %.4f \t %04.2f \t\t %05.4f \t %06.2f \t\t %1d \t %06.2f \t %1d", ...
            trlBlockData(trlCfgIdx, colsToPrint).block, ...
            trlBlockData(trlCfgIdx, colsToPrint).trial, ...
            round( trlBlockData(trlCfgIdx, colsToPrint).stimOri, 2), ...
            round( trlBlockData(trlCfgIdx, colsToPrint).stimDur, 4), ...
            round( trlBlockData(trlCfgIdx, colsToPrint).stimSpread, 2), ...
            round( trlBlockData(trlCfgIdx, colsToPrint).stimContrast, 4), ...
            round( trlBlockData(trlCfgIdx, colsToPrint).reportedOri, 2), ...
            trlBlockData(trlCfgIdx, colsToPrint).reportedConf, ...
            round( trlBlockData(trlCfgIdx, colsToPrint).rawOriError, 2), ...
            trlBlockData(trlCfgIdx, colsToPrint).trlError);
        fprintf("\n")
            
        nCompletedTrialsCurrBlock = nCompletedTrialsCurrBlock + 1;
        currTrlCursor = currTrlCursor + 1; % Move to next to be completed trial
        
    end

    Eyelink('Message','TRIAL BLOCK END');
    
    % =====================================================================
    % Print stats 
    % =====================================================================
    thisBlockData = trlBlockData;
    summaryTable = groupsummary(thisBlockData, {'reportedConf'}, ...
                        {@(x) sqrt(nanmean(x.^2)), 'mean', 'std', 'numel'}, 'rawOriError');
                    
    fprintf('\n\n')
    disp(summaryTable)
        
    save(tMatFile,'trlBlockData');
    
    WaitSecs(0.5);
    
    % --- ASK USER TO CONTINUE TO MAIN EXP ---
    DrawFormattedText(psychToolBoxConfig.w, ...
        'You are now done with the training. \n\nPress SPACE to continue to the main experiment ...', ...
        'center', 'center', 0);
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

end

%% Experiment begin screen
function beginExpScreen(psychToolBoxConfig)

% Before starting the experiment wait for user to press space key
% to begin the experiment.
DrawFormattedText(psychToolBoxConfig.w, 'Press SPACE whenever you are ready to begin...', ...
    'center', 'center', 0);
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
function endOfBlockScreen(psychToolBoxConfig, rewardUSD, thisBlockOptRewardData)

timeoutSecs = 60;  % maximum wait time
tStart = GetSecs;

% TODO: before printing check overconfident or underconfident
% Display end-of-block message
if thisBlockOptRewardData.optimalReward - thisBlockOptRewardData.reward > 2 % if the difference is geater than 2 USD

    if thisBlockOptRewardData.flagUnderconfident
        DrawFormattedText(psychToolBoxConfig.w, ...
            ['End of block.', ...
            '\n\nCurrent reward in USD: $',  num2str(rewardUSD, '%.2f'), ...
            '\n\nYou seem slightly conservative with your confidence report!', ...
            '\nIn this block, you earned $', num2str(thisBlockOptRewardData.reward, '%.2f'), '.', ...
            '\nYou could have earned $', num2str(thisBlockOptRewardData.optimalReward, '%.2f'), ' had you been slightly more confident.', ...
            '\n\n\n\n\nPress SPACE to continue...\n\n(Automatically continuing in 60 seconds)'], ...
            'center', 'center', [0 0 0]);
    elseif thisBlockOptRewardData.flagOverconfident
        DrawFormattedText(psychToolBoxConfig.w, ...
            ['End of block.', ...
            '\n\nCurrent reward in USD: $',  num2str(rewardUSD, '%.2f'), ...
            '\n\nYou seem a bit more confident with your confidence reports!', ...
            '\nIn this block, you earned $', num2str(thisBlockOptRewardData.reward, '%.2f'), '.', ... 
            '\nYou could have earned $', num2str(thisBlockOptRewardData.optimalReward, '%.2f'), ' had you been slightly less overconfident.', ...
            '\n\n\n\n\nPress SPACE to continue...\n\n(Automatically continuing in 60 seconds)'], ...
            'center', 'center', [0 0 0]);
    else
        DrawFormattedText(psychToolBoxConfig.w, ...
            ['End of block.', ...
            '\n\nCurrent reward in USD: $',  num2str(rewardUSD, '%.2f'), ...
            '\nIn this block, you earned $', num2str(thisBlockOptRewardData.reward, '%.2f'), ... 
            '\n\n\n\n\nPress SPACE to continue...\n\n(Automatically continuing in 60 seconds)'], ...
            'center', 'center', [0 0 0]);
    end
        
else
    DrawFormattedText(psychToolBoxConfig.w, ...
        ['End of block.', ...
        '\n\nCurrent reward in USD: $',  num2str(rewardUSD, '%.2f'), ...
        '\n\n\n\n\nPress SPACE to continue...\n\n(Automatically continuing in 60 seconds)'], ...
        'center', 'center', [0 0 0]);
end

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
edfFile = ['CORNFB' num2str(subjectNum,'%.2d')];

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

fprintf("\nScreen info %.2f, %.2f\n", screenWidth, screenHeight);

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
yKey = KbName('y');
nKey = KbName('n');
RestrictKeysForKbCheck([yKey nKey spaceKey quitKey calibrateKey]);

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
function fixationWinCfg = initializeFixationWin(ppd, xCenter, yCenter, black, white, fixationScreenConfig)

fixThreshold = fixationScreenConfig.fixationWinTolerance * ppd;
% fixThreshold = .75 * ppd; % radius of the fixation window
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

% Note: the significant difference between sample mean and actual angle is 
% becuase of asymmetry in the power specturm which is due to the addition 
% of random phase. Moreover, the power spectra of the movie needs to be
% asymmetric otherwise the component will just cancel out leaving zero
% orientation (or noiseless oriented grating).
stimMetrics.meanAngle = meanAngle;
stimMetrics.stdAngle = stdAngle;

% DEBUG
% fprintf("Stim metrics: sampleAngle = %.2f, sampleStdDev = %.2f \n", meanAngle, stdAngle)
end

%% Check fixation 
function [fixationHeld, eyeUsed] = checkFixation(eyeUsedCurr, fixThreshold, tStart, durFix, xCenter, yCenter, quitKey)
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
                % if keyCode(psychToolBoxConfig.calibrateKey)
                %     % Eye used needs to be updated
                %     eyeUsed = setupEyeLink(psychToolBoxConfig.w);
                % 
                %     % Do i need to draw fixation screen again?
                % end
                
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
function respData = showResponseScreen(psychToolBoxConfig, fixationWinCfg, trlCfg, eyeUsed, respScreenGazeHoldDur, beeperDur, respSuccessWaitDur, responseScreenConfig)

% arcRadi1 = 5;
% % arcRadi2 = 6.2; 
% arcRadi2 = 6.8; 

% responseScreenConfig.switchArcs

arcRadi1 = responseScreenConfig.arcRadius1;
arcRadi2 = responseScreenConfig.arcRadius2; 

% If switching is enabled, then switch for every even blocks
if ( responseScreenConfig.switchRespAcs ) && ( mod(trlCfg.block, 2) == 0 )
    arcRadi1 = responseScreenConfig.arcRadius2;
    arcRadi2 = responseScreenConfig.arcRadius1;
end

% Two arcs
redRGBLevel = 255;
greenRGBLevel = 50;
arcRadiusRed   = arcRadi1* psychToolBoxConfig.ppd;
arcRadiusGreen = arcRadi2 * psychToolBoxConfig.ppd;
% arcTolerance   = 0.4 * psychToolBoxConfig.ppd; % +/- tolerance in pixels to match arc
% arcTolerance   = 0.6 * psychToolBoxConfig.ppd;
arcTolerance   = responseScreenConfig.responseArcTolerance * psychToolBoxConfig.ppd;

% Initialize
responseGiven = false;
fixStartTime = NaN;
reportedAngle = NaN;
reportedArc = NaN;
responseTime = NaN;

tRspWinStart = GetSecs;
% tPrint = GetSecs;

% Draw response arc
Screen('FillOval', psychToolBoxConfig.w, fixationWinCfg.fixColor, fixationWinCfg.fixRect);
arcRectRed = CenterRectOnPointd([-1 -1 1 1] * arcRadiusRed, psychToolBoxConfig.xCenter, psychToolBoxConfig.yCenter);
arcRectGreen = CenterRectOnPointd([-1 -1 1 1] * arcRadiusGreen, psychToolBoxConfig.xCenter, psychToolBoxConfig.yCenter);
Screen('FrameArc', psychToolBoxConfig.w, [redRGBLevel 0 0], arcRectRed, 90, -180, 4);   % red arc
Screen('FrameArc', psychToolBoxConfig.w, [0 greenRGBLevel 0], arcRectGreen, 90, -180, 4); % green arc
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
            % angleDeg = mod(atan2d(dy, dx), 180); % Restrict this between 0 and 90 - no need to take mod
            % gazeRadius = sqrt(dx^2 + dy^2);
    
            currentArc = "none";
            angleDeg = atan2d(dy, dx);        
            if angleDeg >= 0 && angleDeg <= 180
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
    Screen('FrameArc', psychToolBoxConfig.w, [redRGBLevel 0 0], arcRectRed, 90, -180, 4);   % red arc
    Screen('FrameArc', psychToolBoxConfig.w, [0 greenRGBLevel 0], arcRectGreen, 90, -180, 4); % green arc
    
    if responseGiven
        % One response is given, draw the cursor on the arc and
        % make beep noise
        dotX = psychToolBoxConfig.xCenter + selectedRadius * cosd(reportedAngle);
        dotY = psychToolBoxConfig.yCenter - selectedRadius * sind(reportedAngle);
        Screen('DrawDots', psychToolBoxConfig.w, [dotX; dotY], 20, [0 0 0], [], 2);
    end
    
    Screen('Flip', psychToolBoxConfig.w);
    
end

confReport = NaN;

if responseGiven
    % fprintf("Reported ARC %s \n", reportedArc)
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

% if responseGiven
%     WaitSecs(respSuccessWaitDur - sampleDuration - beeperDur);
% end

end

%% Show feedback screen for train block
function showTrainTrialFeedbackScreen(psychToolBoxConfig, fixationWinCfg, stimOri, reportedOri, confReport)

arcRadi1 = 5; 
arcRadi2 = 6.2;

% Two arcs
redRGBLevel    = 255;
greenRGBLevel  = 50;
arcRadiusRed   = arcRadi1* psychToolBoxConfig.ppd;
arcRadiusGreen = arcRadi2 * psychToolBoxConfig.ppd;

% Draw response arc
Screen('FillOval', psychToolBoxConfig.w, fixationWinCfg.fixColor, fixationWinCfg.fixRect);
arcRectRed   = CenterRectOnPointd([-1 -1 1 1] * arcRadiusRed, psychToolBoxConfig.xCenter, psychToolBoxConfig.yCenter);
arcRectGreen = CenterRectOnPointd([-1 -1 1 1] * arcRadiusGreen, psychToolBoxConfig.xCenter, psychToolBoxConfig.yCenter);
Screen('FrameArc', psychToolBoxConfig.w, [redRGBLevel 0 0], arcRectRed, 90, -180, 4);   % red arc
Screen('FrameArc', psychToolBoxConfig.w, [0 greenRGBLevel 0], arcRectGreen, 90, -180, 4); % green arc

% Draw response and actual stimuli
if confReport == 0
    selectedRadius = arcRadiusRed;
else
    selectedRadius = arcRadiusGreen;
end

dotX = psychToolBoxConfig.xCenter + selectedRadius * cosd(reportedOri);
dotY = psychToolBoxConfig.yCenter - selectedRadius * sind(reportedOri);
Screen('DrawDots', psychToolBoxConfig.w, [dotX; dotY], 15, [0 0 0], [], 2);

xCenter   = psychToolBoxConfig.xCenter + (selectedRadius - 0.25) * cosd(stimOri);
yCenter   = psychToolBoxConfig.yCenter - (selectedRadius - 0.25) * sind(stimOri);

tickLength = 18;
tickWidth = 5;

rect = CenterRectOnPoint([0 0 tickWidth tickLength], xCenter, yCenter );
     
img = uint8(zeros(tickLength, tickWidth, 3));
img(:,:,1) = 0;
img(:,:,2) = 0;
img(:,:,3) = 200;

tex = Screen('MakeTexture', psychToolBoxConfig.w, img);
Screen('DrawTexture', psychToolBoxConfig.w, tex, [], rect, 90 - stimOri)

% Show response
rawError = reportedOri - stimOri;
rawOriError = mod(rawError + 90, 180) - 90;
rawOriError = floor(rawOriError);

DrawFormattedText(psychToolBoxConfig.w,['\n\n\n\nError: ' num2str(rawOriError,'%d')],'center', ...
            psychToolBoxConfig.yCenter - 30,[0 0 0]);
           
[~, ~] = Screen('Flip', psychToolBoxConfig.w);

end

%% Calculate reward
function reward = calcReward(trueOri, reportedOri, confReport, rewardConfig)

% Note: reported orientation is already pi-periodic, true orientation as
% well
maxTolerableError = rewardConfig.maxTolerableError; % 25; % In degrees
absPerceptualError = abs(reportedOri - trueOri);
absPerceptualError = min(absPerceptualError, 180 - absPerceptualError);

% y1 = 12 + (14 - 12)*abs(sind(2*trueOri)); %12; fprintf("%.2f, %.2f \n", trueOri, y1);
% x1 = 11 + (13 - 11)*abs(sind(2*trueOri));
y1 = rewardConfig.y1; %15; %12; fprintf("%.2f, %.2f \n", trueOri, y1);
x1 = rewardConfig.x1; %14
c1 = rewardConfig.c1; %0.2 % Set to 0.1 if the limit is 100 USD %0.5; % 0.2
m1 = - c1./y1;
m2 = ( m1.*x1 + c1 ) ./ ( x1 - maxTolerableError );
c2 = -  m2.*maxTolerableError;

if numel(confReport) > 1 % it is a vector and not a single value
    % High confidence
    gHC = m1 .* absPerceptualError + c1;
    reward(:) = gHC;
    idx = (absPerceptualError <= y1) & (confReport == 1);
    reward(idx) = gHC(idx);
    idx = (absPerceptualError > y1) & (confReport == 1);
    reward(idx) = rewardConfig.g1*gHC(idx);
    
    % Low confidence
    gLC = m2 .* absPerceptualError + c2;
    reward(confReport == 0) = gLC(confReport == 0);
    
    % Greater than max tolerance
    idx = (absPerceptualError > maxTolerableError) & (confReport == 1);
    reward(idx) = -rewardConfig.g2*c1;
    
    idx = (absPerceptualError > maxTolerableError) & (confReport == 0);
    reward(idx) = 0;
else
    
    if absPerceptualError > maxTolerableError
        if confReport == 1
            reward = -rewardConfig.g2*c1; %-1.6*c1;
        else
            reward = 0;
        end
        return
    end

    if confReport == 1 % High confidence
        gHC = m1 .* absPerceptualError + c1;
        if absPerceptualError > y1
            gHC = rewardConfig.g1*gHC;
        end
        reward = gHC;
    elseif confReport == 0  % Low confidence
        gLC = m2 .* absPerceptualError + c2;
        reward = gLC;
    else
        error('Unknown confidence report: %s. Must be "HC" or "LC".', confReport);
    end

end

end

function points = convertRewardValToPoints(rewardVal)

% 1000 points = 1 USD
points = ceil(1000*rewardVal / 10) * 10;

end

function val = convertPointsToUSD(points)
val = points / 1000;
end

function retData = calcMaxEarnedReward(trueOri, reportedOri, confReport, rewardConfig)

absPerceptualError = abs(reportedOri - trueOri);
absPerceptualError = min(absPerceptualError, 180 - absPerceptualError);

optimalConfReport = zeros(size(absPerceptualError));

optimalConfReport( absPerceptualError >= rewardConfig.x1 ) = 0; % HC
optimalConfReport( absPerceptualError < rewardConfig.x1 )  = 1; % LC

% Reward based on optimal confidence choices
optRewardVec    = calcReward(trueOri, reportedOri, optimalConfReport, rewardConfig);
optPointsVec    = convertRewardValToPoints(optRewardVec);
optimalPoints   = sum(optPointsVec);
optimalReward   = convertPointsToUSD(optimalPoints);
optHCLCRatio    = sum(optimalConfReport == 1) / sum(optimalConfReport == 0);

% Reward based on actual data
rewardVec  = calcReward(trueOri, reportedOri, confReport, rewardConfig);
pointsVec  = convertRewardValToPoints(rewardVec);
points     = sum(pointsVec);
reward     = convertPointsToUSD(points);
HCLCRatio  = sum(confReport == 1) / sum(confReport == 0);

% Overconfident
if ( ( HCLCRatio - optHCLCRatio )/optHCLCRatio ) > 0.1 % Greater than 10% range
    flagOverconfident  = true;
    flagUnderconfident = false;
elseif  ( ( HCLCRatio - optHCLCRatio )/optHCLCRatio ) < 0.1 
    flagOverconfident  = false;
    flagUnderconfident = true;
else
    flagOverconfident  = false;
    flagUnderconfident = false;
end
         
retData.optimalReward      = optimalReward; % Maximum reward that the subject could have earned
retData.reward             = reward;        % Reward that the subject earned
retData.flagOverconfident  = flagOverconfident;
retData.flagUnderconfident = flagUnderconfident;

end

%% Abort trial
function abortTrial(w, durFeedbackFixBreak, beeperDur)

% Give feedback
DrawFormattedText(w, 'Poor fixation!', 'center', 'center');
Screen('Flip', w);
Beeper(200, .4, beeperDur);
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
    % Note: this addition of small spread at all angles maybe reflected in
    % huge differences in actual spread vs sample spread.
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
    env = env + max(env(:))*0.4;
    % env = env + max(env(:))*0.4; % Add small power to all orientations
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

