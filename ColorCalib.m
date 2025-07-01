try 
    
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


    % Stimulus parameters
    flickerHz = 30;              % Flicker frequency
    framesPerCycle = round(frameRate / flickerHz);
    baseRect = [0 0 round( 2 * ppd) round( 2 * ppd)];
    centeredRect = CenterRectOnPointd(baseRect, rect(3)/2, rect(4)/2);
    
    % Color settings
    redColor = [1 0 0];           % Full red
    greenLevel = 0.7;             % Start point for green intensity
    step = 0.01;                  % Step size for adjustment

    % Instructions
    Screen('TextSize', w, 24);
    DrawFormattedText(w, 'Adjust green until flicker disappears.\n\n← / → to adjust, Enter to confirm.', 'center', 'center', [255 255 255]);
    Screen('Flip', w);
    KbWait;

    HideCursor;
    ListenChar(2);

    keepGoing = true;

    while keepGoing

        [~,~,keyCode] = KbCheck(-1);
        if keyCode(KbName('q'))
            Screen('CloseAll');
            ShowCursor;
            ListenChar(0);
            return
        end

        for f = 1:framesPerCycle
            % Alternate between red and green
            if mod(f, 2) == 1
                color = redColor;
            else
                color = [0 greenLevel 0];
            end

            Screen('FillRect', w, color * 255, centeredRect);
            Screen('Flip', w);
        end

        % Check for key input
        [keyIsDown, ~, keyCode] = KbCheck;
        if keyIsDown
            if keyCode(KbName('LeftArrow'))
                greenLevel = max(0, greenLevel - step);
            elseif keyCode(KbName('RightArrow'))
                greenLevel = min(1, greenLevel + step);
            elseif keyCode(KbName('Return'))
                keepGoing = false;
            end
            WaitSecs(0.2); % debounce
        end
    end

    ShowCursor;
    ListenChar(0);
    sca;

    fprintf('Estimated equiluminant green intensity: %.3f\n', greenLevel);

catch ME

    % Notify error
    disp([newline...
        '====================' newline...
        '====== Caught ======' newline...
        '====================' newline])
    
    % Save everything
    save('Caught.mat');

    % Clean up
    Screen('CloseAll');
    ShowCursor;
    ListenChar(0);
    
end
