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


%% Table
interTrlInterval = 1;                                                      % Inter trial interval in seconds
fixationDur = 0.5;                                                         % Fixation duration in seconds
stimOrientations = linspace(0, 180, 10);                                   % 
stimLoc_x = 3;                                                            % Stimulus location in visual field degrees
stimLoc_y = 3;                                                            % Stimulus location in visual field degrees
stimDur = [0.5, 0.5];                                                          % Stimulus duration in seconds
stimSpread = [3, 5];                                                      % Stimulus spread in degrees
stimContrast = [0.5, .9];                                                 % Stimulus contrast levels
respMaxDur = 2;                                                            % Maximum allowed time for user to respond (2 seconds)
numBlocks   = 2;                                                           % Number of blocks 
nTrialsPerBlock = numel(stimOrientations)*numel(stimSpread)*numel(stimContrast)*numel(stimDur);    % Assuming each trial takes max of 5 second, a block should take ~8 minutes
nTrials = numBlocks*nTrialsPerBlock;                                       % Total number of trials to run in this session

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



% Finalize the table
variableNames = {...
    'subject','session','block','trial', 'ITI',... % TODO: experiment information
    'fixationDur', 'stimLocX', 'stimLocY', 'stimDur', 'stimOri', 'stimSpread', 'stimContrast',... % TODO: independent and confounding variables
    'stimSampleMeanOri', 'stimSampleSpread', ...
    'respMaxDur', 'reportedOri','reportedConf','reward','respTime'... % TODO: dependent variables
    };

dat = table( ...
        subjectVec, sessionVec, blockVec, trialVec, itiVec, ...
        fixDurVec, stimLocVec_x, stimLocVec_y, stimDurVec, stimOriVec, stimDispersionVec, stimContrastVec, ...
        stimSampleMeanOriVec, stimSampleSpreadVec, ...
        respMaxDurVec, reportedOriVec, reportedConfVec, rewardVec, respTimeVec, ...
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
    'Response time (time interval between stimulus offset and key press) (in s)';...
    }),...
    'VariableNames',{'description'},'RowNames',variableNames');


%% Define Stimulus constants
try
    %% Eyelink setting
    % Name the .edf file
    edfFile = ['COR' num2str(subjectNum,'%.2d')];

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
    load('/Users/gorislab/Desktop/psychophysics/Calib/2018/humanRig20180612.mat');  % TODO: set according to the system - change this to a common config variable
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
    monitorWidth = 36;  % TODO: set according to the system - change this to a common config variable
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
    leftArrowKey = KbName('LeftArrow');
    rightArrowKey = KbName('RightArrow');
    RestrictKeysForKbCheck([spaceKey quitKey zKey xKey periodKey slashKey leftArrowKey rightArrowKey]);

    % Provide against interference
    HideCursor;
    ListenChar(2);
    PL = MaxPriority(w);
    Priority(PL);
    
    %% For each block run the eyelinksetup again
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
        
    %% Main experiment block goes here
    % TODO: ITI start times goes here
    blockNum = 1; % Dummy variable
    trlNum = 1;   % Dummy variable - later goes into a loop

    row = dat( (dat.session == sessionNum) & (dat.block == blockNum) & (dat.trial == trlNum), :);
    
    %% Fixation
    fixThreshold = .75 * ppd; % radius of the fixation window
    fixDiameter = .2 * ppd; % diameter of the fixation point
    fixOuterRect = CenterRectOnPointd(2*fixDiameter*[0 0 1 1],xCenter,yCenter);
    fixInnerRect = CenterRectOnPointd(fixDiameter*[0 0 1 1],xCenter,yCenter);
    fixRect = [fixOuterRect; fixInnerRect]';
    fixColor = repmat([black white],3,1);
    

    %% Stimuli
    % Either generate stimuli of each combination (in factorial design) or
    % generate it on the fly during each iteration - check if this has any
    % implication for timing

    % Stimulus location
    xStim = row.stimLocX * ppd;
    yStim = row.stimLocY * ppd;

    % Stimulus diameter
    stimDiameter = round( 3 * ppd );
    stimDiameter = stimDiameter + ~mod(stimDiameter,2); % Why does stim diameter have to be even?
    envSpatial = white * envelopeSpatial( stimDiameter, 'raised-cosine', 0.75 ); % Is this like aperture?
    
    % Stimulus parameters
    nFrames = ceil( row.stimDur * frameRate );
    sf = 1.5 * stimDiameter/ppd;                       % cpd * p/ppd = cycle
    Bsf = sf * ( 2^(.5*1) - 1 );                       % 0.1; 1 % cycle * ( 2^(.5*octave) - 1 ) = cycle. Bsf / (diameter/ppd) is in cpd
    Bv = 1.5 * nFrames/frameRate / (stimDiameter/ppd); % dps * f/fps / (p/ppd) = unitless - What is this velocity? Is it bandwidth? More standard would be to define temporal frequency (Hz) or velocity in pixels/frame, not deg/s.
    Vx = 0;                                            % 1 * ppd / frameRate; % dps * ppd / (1/s) = pixels/frame (velocity)
    Vy = 0;                                            % 0 * ppd / frameRate;
    Bth = deg2rad(row.stimSpread);                     % Orientation bandwidth - dispersion
    theta = deg2rad(row.stimOri);                      % Orietation
    sftScale = 1;                                      % spatiotemporal scaling factor
    logGabor = true;                                   % (logical) if true it uses a log-Gabor kernel
    alpha = 1;                                         % exponent for the color envelope
    ft0 = 2;                                           % 2 Hz contrast modulation
    Bft = 0.1;
    counterphase = false;
    
    
    % Spectral envelope
    % Spatial and temporal frequency grids (Fourier space)
    [fx,fy,ft] = getGrids( stimDiameter, stimDiameter, nFrames );
    envSpectral = envelopeSpectral( fx, fy, ft, Vx, Vy, Bv, sf, Bsf, ft0, Bft, sftScale, logGabor, theta, Bth, alpha, counterphase );
    
    
    %% After setting the parameters run the experiment now
    % Draw stimulus
    % Generate stimulus
    movie = randomCloud( envSpectral, false, [], false, counterphase );
    movie = rectify( movie, 'RMS', row.stimContrast, false ) * (white-black) + black;
    movietex = nan(1,nFrames);
    for iF = 1:length(movietex)
        movietex(iF) = Screen('MakeTexture', w, cat(3,movie(:,:,iF), envSpatial));
    end
    
    stimX = row.stimLocX * ppd;
    stimY = row.stimLocY * ppd;
    movieRect = CenterRectOnPointd(stimDiameter*[0 0 1 1], xCenter + stimX, yCenter + stimY);

    % TODO: Wait for the remaining ITI

    % === Fixation point ====
    Screen('FillOval',w,fixColor,fixRect);
    [~,tStart] = Screen('Flip',w);
    % TODO: check if fixation is maintained for given no of secs using
    % eyelink
    WaitSecs(row.fixationDur);  % waits for 0.5 seconds before moving on
    
    
    % TODO: Loop until subject maintains fixation

    % === Show stimulus ====
    % Loop over frames
    for iF = 1:length(movietex)
        Screen('FillOval', w, fixColor, fixRect);
        Screen('DrawTextures', w, movietex(iF), [], movieRect);
        Screen('Flip',w);

        % For each frame make sure fixation is maintained
        % Abort trial if subject breaks fixation
    end
    
    Screen('Flip',w);
    
    % === Response window ====
    % Initial cursor settings
    cursorAngle = 90; % degrees
    
    % Define two arcs: red and green
    arcRadii = [3.5, 5]*ppd;             % in pixels
    arcColors = {[255 0 0], [0 255 0]}; % red and green
    currentArcIndex = 1;               % start on first arc (red)
    
    % Initialize mouse
    % [mouseX0, mouseY0, ~] = GetMouse(w);  % Initial mouse position
    lastSwitchTime = GetSecs;            % for debounce
    
    % Recenter mouse to avoid hitting screen bounds
    SetMouse(xCenter, yCenter, w);

    responseGiven = false;
    while ~responseGiven
    
        % === Clear and redraw arcs ===
        Screen('FillOval', w, fixColor, fixRect);
    
        % Draw both arcs (semi-circles)
        for iArc = 1:2
            radius = arcRadii(iArc);
            arcRect = [ ...
                xCenter - radius, ...
                yCenter - radius, ...
                xCenter + radius, ...
                yCenter + radius ...
            ];
            Screen('FrameArc', w, arcColors{iArc}, arcRect, 90, -180, 2);
        end
    
        % === Compute cursor position on selected arc ===
        radius = arcRadii(currentArcIndex);
        angleRad = deg2rad(cursorAngle);
        cursorX = xCenter + radius * cos(angleRad);
        cursorY = yCenter - radius * sin(angleRad);
    
        % Draw cursor
        cursorSize = 10;
        Screen('FillOval', w, [255 255 255], ...
            [cursorX cursorY cursorX cursorY] + [-1 -1 1 1]*cursorSize);
    
        % Flip screen
        Screen('Flip',w);
    
        % === Mouse input ===
        [mouseX, mouseY, buttons] = GetMouse(w);
        deltaX = mouseX - xCenter;
        deltaY = mouseY - yCenter;
        
        % Recenter mouse to avoid hitting screen bounds
        SetMouse(xCenter, yCenter, w);

        % Controls movement along the arc
        % aVecror = a: cursorX - xCenter, b: cursorY - yCenter
        % bVector = cursorY - yCenter, - (cursorX - xCenter)
        a = cursorX - xCenter;
        b = cursorY - yCenter;
        deltaArc = (b*deltaX - a*deltaY) / sqrt(a^2 + b^2);

        % Controls movement between the arcs
        deltaBetweenArc = - (a*deltaX + b*deltaY) / sqrt(a^2 + b^2);

        % Move cursor along arc using horizontal movement
        sensitivity = 0.3;
        cursorAngle = cursorAngle + sensitivity * deltaArc;
        cursorAngle = max(0, min(180, cursorAngle));
    
        % Switch arcs based on vertical movement (if large enough & cooldown)
        if deltaBetweenArc > 20 && (GetSecs - lastSwitchTime) > 0.2 % Positive deltaY => set to the upper arc
            currentArcIndex = 1;  % toggles between 1 and 2
            lastSwitchTime = GetSecs;
        elseif deltaBetweenArc < -20 && (GetSecs - lastSwitchTime) > 0.2 % Negative deltaY => set to the lower arc
            currentArcIndex = 2;  % toggles between 1 and 2
            lastSwitchTime = GetSecs;
        end
    
        % Confirm selection on mouse click
        if any(buttons)
            responseGiven = true;
        end
    
    end

    % WaitSecs(row.respMaxDur);

    % Stop recording
    Eyelink('StopRecording');
    
    % Close the .edf file
    WaitSecs(1);
    Eyelink('CloseFile');
        
    % Test code
    while 1
        % For emergency exit
        [~,~,keyCode] = KbCheck(-1);
        if keyCode(quitKey)
            Eyelink('ShutDown'); 
            Screen('CloseAll');
            Priority(0);
            ShowCursor;
            ListenChar(0);
            return
        end
    end
    
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
    % phase = zeros(N_X, N_Y, N_frame);
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

