% close all
clear all

% data = load('COR01_Ranjan random ori.mat');
% data = load('COR01_600_trials.mat');
% data = load('COR03.mat');
% data = load('COR04_ArcRad 4.2G 5.4R.mat');
% data = load('COR05 ArcRadi 5G 6.2R.mat');
% data = load('COR06 ArcRadi 6.2G 5R.mat');
% data = load('COR06.mat');
% data = load('COR07.mat');
data = load('COR06.mat');

stimOri = data.dat.stimOri;
reportedOri = data.dat.reportedOri;

absErr = abs(stimOri - reportedOri);
absOriError = min(absErr, 180 - absErr);  % ensures error is in [0, 90]
rawError = stimOri - reportedOri;
rawOriError = mod(rawError + 90, 180) - 90;

data.dat.absOriError = absOriError;
data.dat.rawOriError = rawOriError;

%% Histogram
theoreticalPropHC = zeros(1, numel(unique(data.dat.stimContrast)));
actualPropHC = zeros(1, numel(unique(data.dat.stimContrast)));

tableSummary = groupsummary(data.dat, {'stimContrast', 'reportedConf'}, ...
    {'mean', 'std', 'numel'}, 'rawOriError');
tableSummary = tableSummary(tableSummary.reportedConf == 1, : );
actualPropHC(:) = tableSummary.GroupCount / 120;

[G, contrastLevels, spreadLevels] = findgroups(data.dat.stimContrast, data.dat.stimSpread);
grpIdxes = unique(G);

% Sort conditions by uncertainty difference
[B_, sidx] = sort(spreadLevels + contrastLevels);

% sSpread   = spreadLevels(sidx);
% sContrast = contrastLevels(sidx);

figure
for i=1:numel(sidx)
    gidx_ = sidx(i);

    subplot(2, 2, i)
    hold on

    grpOriErr = rawOriError(G == gidx_);
    histogram(grpOriErr, -90:3:90, "Normalization", "pdf")
    
    % fltOriErr = grpOriErr( zscore(grpOriErr) < 3 ); % Data within 3 SD
    fltOriErr = grpOriErr( abs(grpOriErr) < 50 );
    pd = fitdist(fltOriErr, 'Normal');
    x_values = -90:1:90;
    y = pdf(pd, x_values);
    % Compute from actual data as well - HC and LC
    
    mu = pd.mu;
    sigma = pd.sigma;

    plot(x_values, y, LineWidth=1.5)
    
    xlabel("Orientation (deg)")
    ylabel("count")
    title(sprintf("Contrast: %.3f, Spread: %.3f, \n sigma = %.3f", ...
        contrastLevels(gidx_),  spreadLevels(gidx_), sigma))
    hold off
    
end