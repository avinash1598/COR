clear all
close all

% Luminance
white = 255;
black = 1;
gray = ( white - black ) / 2;

% Get the refresh rate of the screen
frameRate = 75; % [Hz] % This affects the dispersion computation
dt = 1/frameRate;

% Darkroom spec [cm]
screenWidth  = 1280;
screenHeight = 960;
monitorWidth = 36;  % TODO: set according to the system - change this to a common config variable
distance = 81;
ppd = screenWidth/(2*atand(monitorWidth/2/distance)); % Get this from monitor

config.ppd        = ppd;
config.frameRate  = frameRate;
config.white      = white;

% % Plot stimuli in freq domain
% % stimDur, stimSpread, stimContrast
% combinations = [
%     0.1500    5.0000    0.0250
%     0.3000   40.0000    0.0500
%     0.3000    5.0000    0.0150
%     0.3000   25.0000    0.0150
%     0.0500    5.0000    0.0500
%     0.3000    5.0000    0.0500
%     0.0500   25.0000    0.0500
%     0.3000   25.0000    0.0500
% ];

stimOri      = 45;
stimDuration = 0.0500;
stimSpread   = 25.0000;
stimContrast = 0.0500;

stimCfg = generateStimuli(config, stimDuration, stimSpread, stimOri);

% Adjust contrast
movie = randomCloud( stimCfg.envSpectral, false, [], false, stimCfg.counterphase );
movie = rectify( movie, 'RMS', stimContrast, true ) * (white-black) + black;

% % Plot
% figure
% for i=1:size(movie, 3)
%     imagesc(movie(:,:,i), 'AlphaData', stimCfg.envSpatial); 
%     colormap gray
%     set(gca, 'clim', [min(movie(:)), max(movie(:))]); 
%     axis off
%     axis square
%     pause(dt);
% end

% Flatten all variables into vectors
Z = fftshift(fftn(movie));
A = abs(Z); % Amplitude spectrum
A(floor(size(A, 1)/2) + 1, floor(size(A, 2)/2) + 1, floor(size(A, 3)/2) + 1) = 0; % Remove DC

[N_X, N_Y, N_frame] = size(Z);
[fx,fy,ft] = getGrids(N_X, N_Y, N_frame);

envVec = A(:);
fxVec  = fx(:);
fyVec  = fy(:);
ftVec  = ft(:);

figure
subplot(2, 3, 1)
scatter(ftVec, envVec)
xlabel("Temporal freq (a.u.)")
ylabel("Amplitude")

subplot(2, 3, 2)
fr = frequencyRadius(fx, fy, ft, 1, true);
scatter(fr(:), envVec)
xlabel("Radial frequency")
ylabel("Amplitude")

subplot(2, 3, 3)
angleVec = atan2d(fyVec, fxVec);
angleVec = mod(angleVec, 180);
scatter(angleVec, envVec);
xlabel("Angle (deg)")
ylabel("Amplitude")


% % Plot stimuli in freq domain
% % stimDur, stimSpread, stimContrast
% combinations = [
%     0.1500    5.0000    0.0250
%     0.3000   40.0000    0.0500
%     0.3000    5.0000    0.0150
%     0.3000   25.0000    0.0150
%     0.0500    5.0000    0.0500
%     0.3000    5.0000    0.0500
%     0.0500   25.0000    0.0500
%     0.3000   25.0000    0.0500
% ];

%% Stim configuration
function stimCfg = generateStimuli(config, stimDuration, stimSpread, stimOri)

% !!!!!!!!!!
% Beware: this function is completely deterministic. No randomness!
% Changing this will have implications!
% !!!!!!!!!!

% Either generate stimuli of each combination (in factorial design) or
% generate it on the fly during each iteration - check if this has any
% implication for timing

% Stimulus location
ppd = config.ppd;
frameRate = config.frameRate;
white = config.white;

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
    env = env + max(env(:))*0.4; % Add small power to all orientations
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