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

%% Start clean
close all;
clearvars;

% Random seed
rng shuffle

%% Use EyeLink?
dummymode = 0; % variable for eyelink usage. (1 for not using eye-tracking)

%% Initialize the script
% Experiment information
answer = str2double(inputdlg({'Subject number:','Session number:'},'Experiment Information')); % subject number = 0 for pilot
if sum(isnan(answer))
    error('Invalid input.')
end
subjectNum = answer(1);
sessionNum = answer(2);

% Name the .mat file
matFile = ['Test' num2str(subjectNum,'%.2d') '.mat'];

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

%% Table
% Experiment parameter
prior = [1 2]; % 1 = delta prior, 2 = uniform prior
% contrast = [.025 .05 .1]; % Gaussian (Akash 2025-05-15)
% contrast = [.025 .075 .225]; % Gaussian
% contrast = [.015 .045 .135]; % Raised-cosine
% contrast = [.02 .04 .08]; % Tien 2025-05-21
% contrast = [.01 .03 .09]; % Avinash 2025-05-23
% contrast = [.015 .025 .12]; % Akash 2025-06-09
contrast = [.0125 .025 .05]; % Jongmin 06-13
decisionBoundary = 45;
L = 13.5; % [deg]
nAbsTh = 5;
absTh = decisionBoundary + (L/(2*nAbsTh-1):L*2/(2*nAbsTh-1):L);
nRepsInUniformBlock = 2; % determines nTrialsInBlock
nRepsOfDeltaBlock = 4; % determines nTotalReps. Should be a even number to randomize stimulus location across blocks.
nSessions = 4; % should be a factor of nRepsOfDeltaBlock to have full conditions within a session

% Calculate some numbers
nPrior = length(prior);
nContrast = length(contrast);
nTrialsInBlock = nContrast * 2*nAbsTh * nRepsInUniformBlock; % set for uniform block. 2 for CW vs CCW.
nRepsInDeltaBlock = nTrialsInBlock / ( nContrast * 2 ); % 2 for CW vs CCW
nTotalReps = nRepsInDeltaBlock * nRepsOfDeltaBlock; 
nUniformBlocks = nTotalReps / nRepsInUniformBlock;
nBlocks = nPrior * nUniformBlocks;
nBlocksInSession = nBlocks / nSessions;
nTrials = nTrialsInBlock * nBlocks;

% Create table if this is the first session
if sessionNum == 1

    % Simple vectors
    subjectVec = subjectNum*ones(nTrials,1);
    sessionVec = nan(nTrials,1); % fill in sesion number when each session starts
    blockVec = reshape(repmat(1:nBlocksInSession,nTrialsInBlock,nSessions),nTrials,1);
    trialVec = repmat((1:nTrialsInBlock)',nBlocks,1);
    priorVec = reshape(repmat(prior,nTrialsInBlock,nUniformBlocks),nTrials,1); % alternates across blocks

    % Stimulus orientation
    th = decisionBoundary + sort([-(absTh-decisionBoundary) absTh-decisionBoundary]);

    % Confounding variables
    location = [-1 1]; % -1 = left visual field, 1 = right visual field
    
    % Preallocation
    thVec = nan(nTrials,1);
    contrastVec = nan(nTrials,1);
    locationVec = nan(nTrials,1);
    thOrder = nan(nAbsTh,nRepsOfDeltaBlock);

    % Loop over repetitions of full conditions
    for iR = 1:nRepsOfDeltaBlock

        % Shuffle delta block order
        thOrder(:,iR) = randperm(nAbsTh);

        % Loop over absolute values of orientation
        for iTh = 1:nAbsTh

            % Delta block
            thDelta = reshape( repmat( decisionBoundary + (absTh(thOrder(iTh,iR))-decisionBoundary) * [-1 1], nContrast * nRepsInDeltaBlock, 1 ), nTrialsInBlock, 1 );
            contrastDelta = reshape( repmat( contrast, nRepsInDeltaBlock, 2 ), nTrialsInBlock, 1 ); % 2 for CW vs CCW

            % Uniform block
            thUniform = reshape( repmat( th, nContrast * nRepsInUniformBlock, 1 ), nTrialsInBlock, 1 );
            contrastUniform = reshape( repmat( contrast, nRepsInUniformBlock, 2*nAbsTh ), nTrialsInBlock, 1 );

            % Shuffle trial order
            trialOrderDelta = randperm(nTrialsInBlock);
            trialOrderUniform = randperm(nTrialsInBlock);

            % Trial indices for these nPrior blocks
            iT = (iR-1)*nAbsTh*nPrior*nTrialsInBlock+(iTh-1)*nPrior*nTrialsInBlock+1:(iR-1)*nAbsTh*nPrior*nTrialsInBlock+iTh*nPrior*nTrialsInBlock;

            % Finalize vectors
            thVec(iT) = [thDelta(trialOrderDelta); thUniform(trialOrderUniform)];
            contrastVec(iT) = [contrastDelta(trialOrderDelta); contrastUniform(trialOrderUniform)];

        end

    end

    % Loop over repetitions of full conditions
    for iR = 1:nRepsOfDeltaBlock

        % Trial indices for these nBlocks/nRepsOfDeltaBlock blocks
        iT = (iR-1)*nTrials/nRepsOfDeltaBlock+1:iR*nTrials/nRepsOfDeltaBlock;

        % Stimulus location
        if mod(iR,2) == 1

            % Counter-balance within repetition
            locationVecOdd = repmat(location,nTrialsInBlock,nAbsTh);
            locationVecOdd = reshape(locationVecOdd(:,randperm(nPrior*nAbsTh)),[],1);
            locationVec(iT) = locationVecOdd;

        else

            % Preallocation
            locationVecEven = nan(nTrials/nRepsOfDeltaBlock,1);

            % Loop over absolute values of orientation
            for iTh = 1:nAbsTh

                % Delta block: counter-balance across repetitions with the same orientation
                locationVecEven( priorVec(iT) == 1 & abs(thVec(iT)) == absTh( thOrder(:,iR-1) == thOrder(iTh,iR) ) ) = ...
                    -1 * locationVecOdd( priorVec(iT) == 1 & abs(thVec(iT)) == absTh(thOrder(iTh,iR-1)) );

            end

            % Uniform block: counter-balance within repetition
            distLocation = nAbsTh * histcounts( locationVecOdd( priorVec(iT) == 2 ), 2, 'normalization','probability' );
            locationVecEvenUniform = repmat([repmat(location(1),1,distLocation(2)) repmat(location(2),1,distLocation(1))],nTrialsInBlock,1);
            locationVecEven( priorVec(iT) == 2 ) = reshape(locationVecEvenUniform(:,randperm(nAbsTh)),[],1);

            % Finalize the vector
            locationVec(iT) = locationVecEven;

        end
       
    end

    % Finalize the table
    variableNames = {...
        'subject','session','block','trial',... % experiment information
        'prior','location','contrast','th',... % independent and confounding variables
        'response','correct','confidence','reward','rt','rd'... % dependent variables
        };
    dat = array2table(...
        [ ...
        subjectVec sessionVec blockVec trialVec ...
        priorVec locationVec contrastVec thVec ...
        nan(nTrials,6) ...
        ],...
        'VariableNames',variableNames);
    description = table(categorical(...
        {...
        'subject number';...
        'session number (within each subject)';...
        'block number (within each session)';...
        'trial number (within each block)';...
        'prior context: 1 = delta prior, 2 = uniform prior';...
        'stimulus location: -1 = left visual field, 1 = right visual field';...
        'stimulus RMS contrast [sigma]';...
        'stimulus orientation [deg]';...
        'orientation judgment: 1 = clockwise, 0 = counterclockwise';...
        'correctness of orientation judgment: 1 = correct, 0 = incorrect';...
        'confidence judgment: 1 = confident, 0 = not confident';...
        'reward on this trial [USD]';...
        'response time (time interval between stimulus offset and key press) [s]';...
        'response duration (time interval between key press and key release) [s]';...
        }),...
        'VariableNames',{'description'},'RowNames',variableNames');

    % Preallocation
    durBlock = nan(nSessions,nBlocksInSession);
    durSession = nan(nSessions,1);
    
    % Temp
    dat( dat.session == 1, : ) = [dat( dat.session == 1 & dat.prior == 1, : ); dat( dat.session == 1 & dat.prior == 2, : )];

end

% Initial trial for this session
trl = max([1, 1 + find( dat.session == sessionNum-1 & dat.trial == nTrialsInBlock, 1, 'last' )]);

% Fill in session number
dat.session( trl:trl+nTrialsInBlock*nBlocksInSession-1 ) = sessionNum;

%% Try & catch
try
    %% Eyelink setting
    % Name the .edf file
    edfFile = ['Test' num2str(subjectNum,'%.2d')];

    % Exit program if this fails
    if ~EyelinkInit(dummymode, 1) 
        error('Eyelink Init failed.');
    end

    % Open file to record data
    res = Eyelink('Openfile',edfFile);
    if res ~= 0
        error('Cannot create EDF file.');
    end

    % Select which events are saved in the EDF file or available online. Include everything just in case
    Eyelink('Command','file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
    Eyelink('Command','link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');

    % Select which sample data is saved in EDF file or available online. Include everything just in case
    Eyelink('Command','file_sample_data = LEFT,RIGHT,GAZE,HREF,RAW,AREA,GAZERES,BUTTON,STATUS,INPUT');
    Eyelink('Command','link_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT');

    %% Psychtoolbox setting
    % Start timer
    tStartSession = GetSecs;

    % Some preferences
    Screen('Preference','SkipSyncTests', 1);
    oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
    oldSupressAllWarning = Screen('Preference', 'SuppressAllWarnings', 1);
    screenNumber = max(Screen('Screens'));
    doublebuffer = 1;

    % Luminance
    white = 255;
    black = 1;
    gray = ( white - black ) / 2;

    % For color correction
    load('/Users/gorislab/Desktop/psychophysics/Calib/2018/humanRig20180612.mat');
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask','FinalFormatting','DisplayColorCorrection','SimpleGamma');

    % Open window (use PsychImaging('OpenWindow',...)?)
    [w, rect] = PsychImaging('OpenWindow', screenNumber, ...
        gray, [], 32, doublebuffer+1, [], 6);
    
    % Color correction
    PsychColorCorrection('SetEncodingGamma',w,gam.power,'FinalFormatting');
    
    % Get the width and height of the window in pixels
    [screenWidth, screenHeight] = Screen('WindowSize',w);

    % Get the center coordinate of the window
    [xCenter, yCenter] = RectCenter(rect);

    % Get the refresh rate of the screen
    frameRate = 75; % [Hz]

    % Darkroom spec [cm]
    monitorWidth = 36;
    distance = 81;
    ppd = screenWidth/(2*atand(monitorWidth/2/distance));
    
    % Intial flip
    vbl = Screen('Flip',w);

    % Set up blending for drawing smooth points and lines
    Screen('BlendFunction',w,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    
    % Set up audio
    InitializePsychSound();
    PsychPortAudio('Open');

    % Set text size
    Screen('TextSize',w,25);

    % Keyboard setting 
    KbName('UnifyKeyNames');
    spaceKey = KbName('space');
    quitKey = KbName('q');
    zKey = KbName('z');
    xKey = KbName('x');
    periodKey = KbName('.>');
    slashKey = KbName('/?');
    RestrictKeysForKbCheck([spaceKey quitKey zKey xKey periodKey slashKey]);

    % Provide against interference
    HideCursor;
    ListenChar(2);
    PL = MaxPriority(w);
    Priority(PL);
    
    %% Stimuli
    % Timing [s]
    durOrientationDemo = 3;
    durLocationDemo = 1;
    durFix = .5;
    durStim = .5;
    durKeyReleaseTimeout = 1;
    durFeedback = 1;
    durFeedbackFixBreak = .5;
    durITI = 1;
    durIBI = 60;

    % Stimulus location
    xStim = 8 * ppd;
    yStim = -6 * ppd;

    % Spatial envelope
%     sigma = 2/3 * ppd;
%     radius = round( 3 * sigma );
%     diameter = 2*radius + 1;
%     envSpatial = white * envelopeSpatial( diameter, 'Gaussian', sigma );
    diameter = round( 3 * ppd );
    diameter = diameter + ~mod(diameter,2);
    envSpatial = white * envelopeSpatial( diameter, 'raised-cosine', .75 );

    % Stimulus parameters
    nFrames = ceil( durStim * frameRate );
    sf = 1.5 * diameter/ppd; % cpd * p/ppd = cycle
    Bsf = sf * ( 2^(.5*1) - 1 ); % cycle * ( 2^(.5*octave) - 1 ) = cycle. Bsf / (diameter/ppd) is in cpd
    Bv = 1.5 * nFrames/frameRate / (diameter/ppd); % dps * f/fps / (p/ppd) = unitless
    Bth = deg2rad(15);
    ft0 = 1; % spatiotemporal scaling factor
    logGabor = true; % (logical) if true it uses a log-Gabor kernel
    alpha = 1; % exponent for the color envelope

    % Spectral envelope
    [fx,fy,ft] = getGrids( diameter, diameter, nFrames );
    envSpectral = envelopeSpectral( fx, fy, ft, 0, 0, Bv, sf, Bsf, ft0, logGabor, 0, Bth, alpha );
    
    % Fixation
    fixThreshold = .75 * ppd; % radius of the fixation window
    fixDiameter = .2 * ppd; % diameter of the fixation point
    fixOuterRect = CenterRectOnPointd(2*fixDiameter*[0 0 1 1],xCenter,yCenter);
    fixInnerRect = CenterRectOnPointd(fixDiameter*[0 0 1 1],xCenter,yCenter);
    fixRect = [fixOuterRect; fixInnerRect]';
    fixColor = repmat([black white],3,1);

    % Monetary reward [USD]
    baselineReward = 10;
    payoff = [.05 -.05; .03 -.01]; % optimal criterion: (.05-.01)/((.05-.01)+(.05-.03)) = 2/3

    % For progress report
    responseString = {'CW','CCW'};
    correctString = {'O','X'};
    confidenceString = {'Yes','No'};
    
    %%
    % Spatial envelope
%     sigma = 2/3 * ppd;
%     radius = round( 3 * sigma );
%     diameter = 2*radius + 1;
%     envSpatialGaussian = white * envelopeSpatial( diameter, 'Gaussian', sigma );
%     
%     % Generate the motion cloud
%     nFramesDemo = ceil( durStim * frameRate );
%     BvDemo = 3 * nFramesDemo/frameRate / (diameter/ppd); % dps * f/fps / (p/ppd) = unitless
%     [fx,fy,ft] = getGrids( diameter, diameter, nFramesDemo );
%     envSpectralDemo = envelopeSpectral( fx, fy, ft, 0, 0, BvDemo, sf, Bsf, ft0, logGabor, 0, Bth, alpha );
%     movie = randomCloud( envSpectralDemo, false, [], false );
%  
%     % Finalize
%     movieRectified = rectify( movie, 'RMS', .025, false );
%     movieGaussian =( ( movieRectified(:,:,1) -.5 ) .* envSpatialGaussian/white + .5) * (white-black) + black;
%     
%     % Spatial envelope
%     diameter = round( 3 * ppd );
%     diameter = diameter + ~mod(diameter,2);
%     envSpatialRaisedCosine = white * envelopeSpatial( diameter, 'raised-cosine', .25 );
% 
%     % Generate the motion cloud
%     nFramesDemo = ceil( durStim * frameRate );
%     BvDemo = 3 * nFramesDemo/frameRate / (diameter/ppd); % dps * f/fps / (p/ppd) = unitless
%     [fx,fy,ft] = getGrids( diameter, diameter, nFramesDemo );
%     envSpectralDemo = envelopeSpectral( fx, fy, ft, 0, 0, BvDemo, sf, Bsf, ft0, logGabor, 0, Bth, alpha );
%     movie = randomCloud( envSpectralDemo, false, [], false );
%     
%     % Finalize
%     movieRectified = rectify( movie, 'RMS', .015, false );
%     movieRaisedCosine = ( ( movieRectified(:,:,1) -.5 ) .* envSpatialRaisedCosine/white + .5) * (white-black) + black;
% 
%     % Demo
%     figure
%     subplot(121); mesh(movieGaussian,[0 255]); axis square; 
%     subplot(122); mesh(movieRaisedCosine,[0 255]); axis square; 
%     mesh(round(movie(:,:,1)),[0 255]); colormap gray; axis square; 
%     imagesc(round(movie(:,:,1)),[0 255]); colormap gray; axis square; axis off
%     image(round(movie(:,:,1))/255); colormap gray; axis square; axis off
    % Demo
%     figure('color',[.5 .5 .5])
%     for iF = 1:5:nFramesDemo
%         image(movie(:,:,iF)+2), colormap gray, hold on, box off, axis off, axis square
%         pause(5/frameRate)
%     end
    
    %% Experiment loop
    % Loop over blocks
%     for blockNum = 1:nBlocksInSession
%     % Check the uniform block
% %     trl = trl + nTrialsInBlock;
% %     for blockNum = 2:nBlocksInSession
% 
%         % -------
%         % Eyelink
%         % -------
%     
%         % Eyelink setting
%         el = EyelinkInitDefaults(w);
%         
%         % Calibration & validation
%         EyelinkDoTrackerSetup(el);
%     
%         % Start recording
%         WaitSecs(0.05);
%         Eyelink('StartRecording');
%     
%         % Check which eye is available
%         eyeUsed = Eyelink('EyeAvailable'); % 0 = left, 1 = right, 2 = binocular
%         % Get samples from right eye if binocular
%         if eyeUsed == 2
%             eyeUsed = 1;
%         end
%     
%         % ----
%         % Demo
%         % ----
%     
%         % Stimulus orientation in this block
%         th = unique(dat.th( dat.session == sessionNum & dat.block == blockNum ));
%         
%         % Generate stimulus
%         X = meshgrid((-.5*diameter:.5*diameter-1) + mod(diameter,2) * .5);
%         sineWave = sin( sf/diameter * 2*pi * X );
%         squareWave = (ceil(abs(sineWave)) .* sign(sineWave) - .5) * 2*max(dat.contrast) * gray + gray;
%         squareWaveTex = Screen('MakeTexture',w,cat(3,squareWave,envSpatial));
% 
%         % Destination rectangle for orientation demo
%         distDemo = screenWidth/(2*nAbsTh);
%         xDemo = linspace(.5*distDemo,(2*nAbsTh-.5)*distDemo,2*nAbsTh) - .5*screenWidth;
%         if dat.prior(trl) == 1 % delta prior
%             thisTheta = find(absTh == max(th));
%             demoRect = nan(4,2);
%             demoRect(:,1) = CenterRectOnPointd(distDemo*[0 0 1 1],xCenter+xDemo(nAbsTh+1-thisTheta),yCenter);
%             demoRect(:,2) = CenterRectOnPointd(distDemo*[0 0 1 1],xCenter+xDemo(nAbsTh+thisTheta),yCenter);
%         elseif dat.prior(trl) == 2 % uniform prior
%             demoRect = nan(4,2*nAbsTh);
%             for iTh = 1:2*nAbsTh
%                 demoRect(:,iTh) = CenterRectOnPointd(distDemo*[0 0 1 1],xCenter+xDemo(iTh),yCenter);
%             end
%         end
% 
%         % Get ready for orientation demo
%         DrawFormattedText(w,'Orientation','center','center',white);
%         [~,tStartBlock] = Screen('Flip',w);
%         KbStrokeWait(-1);
% 
%         % Orientation demo
%         Screen('DrawTextures',w,squareWaveTex,[],demoRect,th);
%         Screen('Flip',w);
%         WaitSecs(durOrientationDemo);
% 
%         % Get ready for location demo
%         DrawFormattedText(w,'Location','center','center',white);
%         Screen('Flip',w);
%         KbStrokeWait(-1);
%     
%         % Stimulus location in this block
%         location = dat.location(trl); % -1 = left, 1 = right
%         movieRect = CenterRectOnPointd(diameter*[0 0 1 1],xCenter+location*xStim,yCenter+yStim);
%         
%         % Location Demo
%         Screen('FillOval',w,fixColor,fixRect);
%         Screen('FrameOval',w,[],movieRect,.1*ppd);
%         Screen('Flip',w);
%         WaitSecs(durLocationDemo);
% 
%         % ----------
%         % Trial loop
%         % ----------
% 
%         % Get ready for the main experiment
%         [~,~,textRect] = DrawFormattedText(w,'Ready','center','center',white);
%         textHeight = textRect(4) - textRect(2); % for reward feedback
%         Screen('Flip',w);
%         KbStrokeWait(-1);
% 
%         % Send event to eyelink
%         Eyelink('Message','START');
%     
%         % For progress report
%         fprintf('\n\nSubject\tSession\tBlock\tTrial\tCont\tStim\tResp\tO/X\tConf\tReward\tTotal\tRT\tRD\n')
%     
%         % Initial value
%         nCompletedTrials = 0;
%     
%         % Start timer for initial ITI
%         tStart = GetSecs;
%     
%         % Loop over trials
%         while nCompletedTrials < nTrialsInBlock
%     
%             % Parameters on this trial
%             contrast = dat.contrast(trl); % [0-1]
%             th = dat.th(trl); % [deg]
%     
%             % Generate stimulus
%             movie = randomCloud( envSpectral, false, [], false );
%             movie = rectify( movie, 'RMS', contrast, false ) * (white-black) + black;
%             movietex = nan(1,nFrames);
%             for iF = 1:length(movietex)
%                 movietex(iF) = Screen('MakeTexture',w,cat(3,movie(:,:,iF),envSpatial));
%             end
%     
%             % Default value
%             fixated = true;
%     
%             % Inter-trial interval
%             WaitSecs( durITI - ( GetSecs - tStart ) );
%     
%             % Fixation point
%             Screen('FillOval',w,fixColor,fixRect);
%             [~,tStart] = Screen('Flip',w);
%             Eyelink('Message','FIXATION');
%             
%             % Loop until the subject maintain fixation
%             while 1
%                 % For emergency exit
%                 [~,~,keyCode] = KbCheck(-1);
%                 if keyCode(quitKey)
%                     Eyelink('ShutDown'); 
%                     Screen('CloseAll');
%                     Priority(0);
%                     ShowCursor;
%                     ListenChar(0);
%                     return
%                 % Check fixation
%                 elseif Eyelink('NewFloatSampleAvailable') > 0
%                     evt = Eyelink('NewestFloatSample');
%                     gx = evt.gx(eyeUsed+1); % [left eye gaze x, right eye gaze x] +1 because this is matlab
%                     gy = evt.gy(eyeUsed+1);
%                     if ~isnan(gx) && ~isnan(gy)
%                         % If fixation was maintained for durFix, break while loop to show stimulus
%                         if norm([xCenter,yCenter]-[gx,gy]) <= fixThreshold
%                             if GetSecs - tStart >= durFix
%                                 break
%                             end
%                         % Reset timer if the subject break fixation
%                         else
%                             tStart = GetSecs;
%                         end
%                     end
%                 end
%                 % Avoid hogging CPU
%                 WaitSecs(0.001); % ~1000 Hz polling
%             end
%             
%             % Loop over frames
%             for iF = 1:length(movietex)
%                 Screen('FillOval',w,fixColor,fixRect);
%                 Screen('DrawTextures',w,movietex(iF),[],movieRect,th);
%                 Screen('Flip',w);
%                 % Send event to eyelink
%                 if iF == 1
%                     Eyelink('Message','STIMULUS');
%                 end
%                 % Check fixation
%                 if Eyelink('NewFloatSampleAvailable') > 0
%                     evt = Eyelink('NewestFloatSample');
%                     gx = evt.gx(eyeUsed+1); % [left eye gaze x, right eye gaze x] +1 because this is matlab
%                     gy = evt.gy(eyeUsed+1);
%                     % Abort trial if subject break fixation
%                     if norm([xCenter,yCenter]-[gx,gy]) > fixThreshold
%                         fixated = false;
%                         break
%                     end
%                 end
%             end
%             
%             % If fixation was maintained during stimulus presentation
%             if fixated
%   
%                 % Orientation and confidence judgment
%                 [~,tStart] = Screen('Flip',w);
%     
%                 % Wait until all keys are released
%                 while KbCheck(-1); end
%     
%                 % Wait for key press
%                 while 1
%                     [keyIsDown,tKeyDown,keyCode] = KbCheck(-1);
%                     if keyIsDown && sum(keyCode([zKey xKey periodKey slashKey])) == 1 % ignore if subject presses more than one key at the same time
%                         % Send event to eyelink
%                         Eyelink('Message','RESPONSE');
%                         % Compute response time
%                         rt = tKeyDown - tStart;
%                         break
%                     end
%                 end
%     
%                 % Wait for key release
%                 while 1
%                     [keyIsDown,tKeyUp] = KbCheck(-1);
%                     if ~keyIsDown
%                         % Compute response duration
%                         rd = tKeyUp - tKeyDown;
%                         break
%                     % Break while loop if subject never release key within timeout
%                     elseif GetSecs - tKeyDown > durKeyReleaseTimeout
%                         rd = nan;
%                         break
%                     end
%                 end
%     
%                 % Map orientation and confidence judgments
%                 if keyCode(zKey) || keyCode(xKey) 
%                     response = 0; % CCW
%                     confidence = keyCode(zKey); % 'z' = confident, 'x' = not confidnet
%                 elseif keyCode(periodKey) || keyCode(slashKey) 
%                     response = 1; % CW
%                     confidence = keyCode(slashKey); % '/' = confident, '.' = not confident
%                 end
%                 correct = response == (th > 45);
%     
%                 % Compute reward
%                 reward = payoff(2-confidence,2-correct);
%                 totalReward = baselineReward + sum(dat.reward(1:max(1,trl-1)),'omitnan') + reward;
%     
%                 % Give feedback
%                 if reward >= 0
%                     DrawFormattedText(w,['+' num2str(reward,'%.2f')],'center',yCenter-textHeight,[0 white 0]);
%                 else
%                     DrawFormattedText(w,['-' num2str(abs(reward),'%.2f')],'center',yCenter-textHeight,[white 0 0]);
%                 end
%                 DrawFormattedText(w,['$' num2str(totalReward,'%.2f')],'center',yCenter+textHeight,white);
%                 Screen('Flip',w);
%                 Eyelink('Message','FEEDBACK');
%                 WaitSecs(durFeedback);
% 
%                 % Start timer for inter-trial interval
%                 [~,tStart] = Screen('Flip',w);
%     
%                 % Record everything
%                 dat.rt(trl) = rt;
%                 dat.rd(trl) = rd;
%                 dat.response(trl) = response;
%                 dat.correct(trl) = correct;
%                 dat.confidence(trl) = confidence;
%                 dat.reward(trl) = reward;
%                 
%                 % Report progress
%                 fprintf('%d\t%d\t%d\t%d\t%.3f\t%.1f°\t%s\t%s\t%s\t$%.2f\t$%.2f\t%.1f s\t%.1f s\n', ...
%                     subjectNum, sessionNum, blockNum, dat.trial(trl), ...
%                     contrast, th, responseString{2-response}, correctString{2-correct}, confidenceString{2-confidence}, ...
%                     reward, totalReward, rt, rd);
%     
%                 % Count the number of completed trials
%                 nCompletedTrials = nCompletedTrials + 1;
%             
%             % If the subject broke fixation
%             else
%     
%                 % Give feedback
%                 DrawFormattedText(w,'Poor fixation!','center','center');
%                 Screen('Flip',w);
%                 Beeper(400,.4,.15);
%                 WaitSecs( durFeedbackFixBreak - .15 );
%     
%                 % Start timer for inter-trial interval
%                 [~,tStart] = Screen('Flip',w);
%     
%                 % Adjust trial number
%                 remainingTrials = dat.session == sessionNum & dat.block == blockNum & dat.trial > dat.trial(trl);
%                 dat.trial( remainingTrials ) = dat.trial( remainingTrials ) - 1;
%                 dat.trial(trl) = nTrialsInBlock;
%     
%                 % Append this trial to the end of the block
%                 dat = [dat(1:trl+nTrialsInBlock-nCompletedTrials-1,:); dat(trl,:); dat(trl+nTrialsInBlock-nCompletedTrials:end,:)];
%     
%                 % Mark this trial as aborted
%                 dat.trial(trl) = nan;
%     
%                 % Report progress
%                 disp('Poor fixation!')
%     
%             end
%     
%             % Trial count
%             trl = trl + 1;
%             
%         end
% 
%         % -----------------
%         % Wrap up the block
%         % -----------------
%     
%         % Report statistics
%         thisBlock = dat.session == sessionNum & dat.block == blockNum;
%         nAbortedTrials = sum(isnan(dat.trial(thisBlock)));
%         pCorrect = mean(dat.correct(thisBlock),'omitnan');
%         pConfident = mean(dat.confidence(thisBlock),'omitnan');
%         totalRewardInBlock = sum(dat.reward(thisBlock),'omitnan');
%         durBlock(sessionNum,blockNum) = GetSecs - tStartBlock;
%         fprintf('\nThe number of aborted trials: %d\n', nAbortedTrials);
%         fprintf('Proportion correct: %.2f\n', pCorrect);
%         fprintf('Proportion confident: %.2f\n', pConfident);
%         fprintf('Reward earned in this block: $%.2f\n', totalRewardInBlock);
%         fprintf('Total reward: $%.2f\n', totalReward);
%         fprintf('Elapsed time: %d min %d sec\n\n', floor(durBlock(sessionNum,blockNum)/60), round(mod(durBlock(sessionNum,blockNum),60)))
% 
%         % Send event to eyelink
%         Eyelink('Message', 'END');
%     
%         % Stop recording
%         Eyelink('StopRecording');
%     
%         % Close the .edf file
%         WaitSecs(1);
%         Eyelink('CloseFile');
%         
%         % Inter-block interval
%         if blockNum ~= nBlocksInSession
% 
%             % Start timer
%             tStart = GetSecs;
% 
%             % Loop until durIBI has passed or subject wants to start
%             while 1
%                 
%                 % Compute remaining time
%                 remainingTime = ceil( durIBI - ( GetSecs - tStart ) );
%                 if remainingTime <= 0
%                     break
%                 end
% 
%                 % Display countdown
%                 countdown = sprintf('The next block will start in %d seconds.\nPress any key to start now.', remainingTime);
%                 DrawFormattedText(w,countdown,'center','center',white);
%                 Screen('Flip',w);
% 
%                 % Check if subject presses a key to start early
%                 keyIsDown = KbCheck(-1);
%                 if keyIsDown
%                     break
%                 end
%                 
%                 % Avoid hogging CPU
%                 WaitSecs(0.1);
% 
%             end
% 
%         end
% 
%     end

    % -------------------
    % Wrap up the session
    % -------------------

%     % Report statistics
%     thisSession = dat.session == sessionNum;
%     nAbortedTrials = sum(isnan(dat.trial(thisSession)));
%     pCorrect = mean(dat.correct(thisSession),'omitnan');
%     pConfident = mean(dat.confidence(thisSession),'omitnan');
%     totalRewardInSession = sum(dat.reward(thisSession),'omitnan');
%     durSession(sessionNum) = GetSecs - tStartSession;
%     fprintf('\nProportion aborted: %.2f\n', nAbortedTrials/(nAbortedTrials+sum(thisSession)));
%     fprintf('Proportion correct: %.2f\n', pCorrect);
%     fprintf('Proportion confident: %.2f\n', pConfident);
%     fprintf('Reward earned in this block: $%.2f\n', totalRewardInSession);
%     fprintf('Total reward: $%.2f\n', totalReward);
%     fprintf('Elapsed time: %d min %d sec\n\n', floor(durSession(sessionNum)/60), round(mod(durSession(sessionNum),60)))

    % Announce the end of the session
    DrawFormattedText(w,'The end','center','center',white);
    Screen('Flip',w);
    WaitSecs(3);

    % Save the .edf file
%     fprintf('\nReceiving data file ''%s''\n', edfFile);
%     Eyelink('ReceiveFile', edfFile, pwd, 1);
%     fprintf('Data file ''%s'' received\n', edfFile);
        
    % Save the .mat file
    save(matFile,'dat','description','durBlock','durSession');

    % Clean up
    Eyelink('ShutDown'); 
    Screen('CloseAll');
    Priority(0); % Shutdown realtime mode
    ShowCursor; % Show cursor again, if it has been disabled.
    ListenChar(0);

catch ME

    % Notify error
    disp([newline...
        '====================' newline...
        '====== Caught ======' newline...
        '====================' newline])

    % Save everything
    save('Caught.mat');
    
    % Clean up
    Eyelink('ShutDown'); 
    Screen('CloseAll');
    Priority(0); % Shutdown realtime mode
    ShowCursor; % Show cursor again, if it has been disabled.
    ListenChar(0);
    
end


%% ----------------------------- FUNCTION ----------------------------- %%

%% envelopeSpatial
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
env( (X.^2 + Y.^2) > (N/2)^2 ) = 0;
    
end

%% getGrids
function [fx, fy, ft] = getGrids(N_X, N_Y, N_frame)

% [fx, fy, ft] = meshgrid(-N_X/2:(N_X-1)/2+1, -N_Y/2:(N_Y-1)/2+1, -N_frame/2:(N_frame-1)/2+1);

% Even and odd
df = 1;
x = (-N_X/df/2:df:N_X/df/2-1) + mod(N_X/df,2) * df/2;
y = (-N_Y/df/2:df:N_X/df/2-1) + mod(N_Y/df,2) * df/2;
t = (-N_frame/1/2:1:N_frame/1/2-1) + mod(N_frame/1,2) * 1/2;
[fx, fy, ft] = meshgrid(x,y,t);

end

%% frequencyRadius
function fr = frequencyRadius(fx, fy, ft, ft_0, clean_division)

if ft_0 == inf
    fr2 = fx.^2 + fy.^2;
else
    fr2 = fx.^2 + fy.^2 + (ft/ft_0).^2;
end

if clean_division
    fr2(fr2==0) = inf;
end

fr = sqrt(fr2);

end

%% envelopeColor
function env = envelopeColor(fx, fy, ft, alpha, ft_0)

if alpha == 0 % white noise
    env = ones(size(fx));
else
    fr = frequencyRadius(fx, fy, ft, ft_0, true);
    env = fr.^(-alpha);
end

end

%% envelopeRadial
function env = envelopeRadial(fx, fy, ft, sf_0, B_sf, ft_0, loggabor)

if sf_0 == 0 || B_sf == inf
    if loggabor
        env = envelopeColor(fx, fy, ft, 1, ft_0); % should I remove this?
    else
        env = ones(size(fx));
    end
elseif loggabor
    fr = frequencyRadius(fx, fy, ft, ft_0, true);
    env = exp(-.5*log(fr/sf_0).^2/log(1+B_sf/sf_0)^2);
    % env = exp(-.5*(log(fr)-log(sf_0^2/sqrt(sf_0^2+B_sf^2))).^2/log(1+B_sf^2/sf_0^2));
else
    fr = frequencyRadius(fx, fy, ft, ft_0, true);
    env = exp(-.5*(fr-sf_0).^2/B_sf^2);
end

end

%% envelopeSpeed
function env = envelopeSpeed(fx, fy, ft, V_X, V_Y, B_V, ft_0)

if size(ft,3) == 1
    env = ones(size(fx));
elseif B_V == 0
    env = zeros(size(fx));
    env(ft == 0) = 1;
else
    fr = frequencyRadius(fx, fy, ft, ft_0, true);
    env = exp(-.5*((ft+fx*V_X+fy*V_Y)).^2./(B_V*fr).^2);
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
    angl = atan2(fy, fx);
    env = exp(cos(2*(angl-theta))/(2*B_theta)^2);
end

end

%% envelopeSpectral
function envelope = envelopeSpectral(fx, fy, ft, V_X, V_Y, B_V, sf_0, B_sf, ft_0, loggabor, theta, B_theta, alpha)

% envelope = envelopeColor(fx, fy, ft, alpha, inf) ...
envelope = envelopeColor(fx, fy, ft, alpha, ft_0) ...
    .* envelopeOrientation(fx, fy, theta, B_theta) ...
    .* envelopeRadial(fx, fy, ft, sf_0, B_sf, ft_0, loggabor) ...
    .* envelopeSpeed(fx, fy, ft, V_X, V_Y, B_V, ft_0);

end

%% randomCloud
function z = randomCloud(envelope, impulse, events, do_amp)

[N_X, N_Y, N_frame] = size(envelope);
if impulse
    [fx, fy, ft] = getGrids(N_X, N_Y, N_frame);
    phase = -2*pi*(N_X/2*fx + N_Y/2*fy + N_frame/2*ft);
    F_events = exp(1i * phase);
elseif isempty(events)
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
