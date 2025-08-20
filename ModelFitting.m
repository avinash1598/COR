close all
clear all

%% Load data
data = load('COR06 ArcRadi 6.2G 5R.mat');
stimOri = data.dat.stimOri;
reportedOri = data.dat.reportedOri;

absErr = abs(stimOri - reportedOri);
absOriError = min(absErr, 180 - absErr);  % ensures error is in [0, 90]
rawError = stimOri - reportedOri;
rawOriError = mod(rawError + 90, 180) - 90;

data.dat.absOriError = absOriError;
data.dat.rawOriError = rawOriError;

%%
sigma_e = 12.24;
sigma_i = 1;
varGain = 0.5;
sampleCnt = 100;
% x = linspace(-90, 90, 1001);
% 
% shapeParam = 1./varGain;
% scaleParam = varGain;
% gammaSamples = gaminv(linspace(1/sampleCnt, 1 - 1/sampleCnt, sampleCnt), ...
%     shapeParam, scaleParam);
% sigma_theta = sqrt( sigma_e^2 + (sigma_i*gammaSamples).^2 );
% 
% y = normpdf(x, zeros(size(sigma_theta')), sigma_theta');
% mixturePDF = sum(y, 1);
% 
% 
% figure
% plot(x, mixturePDF)
% xlabel("X")
% ylabel("PDF")

%% MLE estimation

% grpIDx = 2;
% [G, contrastLevels] = findgroups(data.dat.stimContrast);
% grpIdxes = unique(G);
% grpOriErr = rawOriError(G == grpIdxes(grpIDx));
% dataPts = -90:3:90;
% 
% confReport = data.dat.reportedConf;
% grpConfReport = confReport(G == grpIdxes(grpIDx));
% grpOriErrHC = grpOriErr(grpConfReport == 1);
% grpOriErrLC = grpOriErr(grpConfReport == 0);

modelData = load('modelContOriData.mat');
grpOriErr = modelData.data.err;
dataPts = -10:0.1:10;

% confReport = modelData.data.confReport;
% grpOriErrHC = grpOriErr(confReport == 1);
% grpOriErrLC = grpOriErr(confReport == 0);


% options = statset( ...
%     'MaxIter', 1e4, ...
%     'MaxFunEvals', 1e5, ...
%     'TolFun', 1e-6, ...
%     'TolX', 1e-12, ...
%     'Display', 'iter');  % show progress

options = optimset('Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-12);
x0 = [rand(1, size(grpOriErr, 1)), rand, rand]; 
lb = [zeros(1, size(grpOriErr, 1)), 0, 0];       % lower bounds
ub = [];                                         % no upper bounds unless needed

optimalValues = fmincon(@(x) totalNLL(x, grpOriErr), ...
    x0, [], [], [], [], lb, ub, [], options);

% % True s_e = 10, s_i = 1, v_G = 0.5
% paramEsts = mle(grpOriErr(:), ...
%     'pdf', @myCustomPDF, ...
%     'start', [rand(1, size(grpOriErr, 1)), rand, rand], ... % Later include mean as well
%     'lowerbound', [zeros(1, size(grpOriErr, 1)), 0, 0], ...
%     'Options', options); % 'NTrials', 1000

% paramEsts = mle(grpOriErr(10,:), ...
%     'pdf', @myCustomPDF, ...
%     'start', [rand, rand, rand], ... % Later include mean as well
%     'lowerbound', [0, 0, 0], ...
%     'Options', options); % 'NTrials', 1000


% % Plot results
% s_e = paramEsts(1);
% s_i = paramEsts(2);
% v_G = paramEsts(3);
% % s_e = 2;
% % s_i = 1;
% % v_G = 0.5;
% % x = linspace(-90, 90, 1001);
% x = linspace(-10, 10, 1001);

% shapeParam = 1./v_G;
% scaleParam = v_G;
% gammaSamples = gaminv(linspace(1/100, 1 - 1/100, 100), ...
%     shapeParam, scaleParam);
% sigma_theta = sqrt( s_e^2 + (s_i*gammaSamples).^2 );
% 
% y = normpdf(x, zeros(size(sigma_theta')), sigma_theta');
% mixturePDF = sum(y, 1);
% mixturePDF = mixturePDF./trapz(x, mixturePDF);

% % Find Cc and sigma_m
% options = optimset('Display', 'iter', 'TolFun', 1e-6, 'TolX', 1e-12);
% targetResp.errHC = grpOriErrHC;
% targetResp.errLC = grpOriErrLC;
% x0 = [rand, rand]; % sigma_m, Cc
% 
% optimalValues = fminsearch(@(x) minimizeError(x, targetResp, paramEsts), x0, options);
% 
% % sigma_m          = 0.2;
% % Cc               = 0.5;
% sigma_m = optimalValues(1);
% Cc = optimalValues(2);

% Plot
% figure
% histogram(grpOriErr, dataPts, "Normalization", "pdf")
% hold on
% plot(x, mixturePDF, LineWidth=2, DisplayName="Fit")
% xlabel("theta error")
% ylabel("PDF")
% hold off


%% Custom PDF
% function p = myCustomPDF(x, varargin)
% 
% % Determine how many sigma_e components there are
% nSigmaE = numel(varargin) - 2; % last two are sigma_i and varGain
% sigma_e = cell2mat(varargin(1:nSigmaE));
% sigma_i = varargin{end-1};
% varGain = varargin{end};
% 
% % dataPts = -90:3:90;
% dataPts = -10:0.1:10;
% sampleCnt = 100;
% shapeParam = 1./varGain;
% scaleParam = varGain;
% gammaSamples = gaminv(linspace(1/sampleCnt, 1 - 1/sampleCnt, sampleCnt), ...
%     shapeParam, scaleParam);
% sigma_theta = sqrt( sigma_e'.^2 + (sigma_i*gammaSamples).^2 );
% data = reshape(x, nSigmaE, []);
% 
% totalProb = [];
% 
% for i=1:nSigmaE
%     sigma_theta_i = sigma_theta(i, :);
%     y = normpdf(dataPts, zeros(size(sigma_theta_i')), sigma_theta_i');
%     mixturePDF = sum(y, 1);
%     normVal = trapz(dataPts, mixturePDF);
% 
%     p_X = exp(- (data(i, :)').^2 ./ (2*sigma_theta_i.^2 ) ) ./ sqrt(2*pi*sigma_theta_i.^2);
%     p_ = sum(p_X, 2)./normVal;
% 
%     totalProb = [totalProb; p_'];
% end
% 
% p = totalProb(:);
% 
% end 

function p=myCustomPDF(x, sigma_e, sigma_i, varGain)

% dataPts = -90:3:90;
dataPts = -10:0.1:10;
sampleCnt = 100;
shapeParam = 1./varGain;
scaleParam = varGain;
gammaSamples = gaminv(linspace(1/sampleCnt, 1 - 1/sampleCnt, sampleCnt), ...
    shapeParam, scaleParam);
sigma_theta = sqrt( sigma_e'.^2 + (sigma_i*gammaSamples).^2 );

y = normpdf(dataPts, zeros(size(sigma_i')), sigma_i');
mixturePDF = sum(y, 1);
normVal = trapz(dataPts, mixturePDF);

p_X = exp(- (x').^2 ./ (2*sigma_theta.^2 ) ) ./ sqrt(2*pi*sigma_theta.^2);
p = sum(p_X, 2)./normVal;

end

function nll = totalNLL(params, data)

nSigmaE = numel(params) - 2; 
sigma_e = params(1:nSigmaE);
sigma_i = params(end-1);
varGain = params(end);

% nSigmaE = size(data, 1);
% sigma_e = params.sigma_e;
% sigma_i = params.sigma_i;
% varGain = params.varGain;

nll = 0;

for k = 1:nSigmaE
    x = data(k, :);
    p = myCustomPDF(x, sigma_e(k), sigma_i, varGain);
    p(p <= 0 | isnan(p)) = eps;
    nll = nll - sum(log(p));
end

end

function errorVal = minimizeError(x, targetResp, params)

errHC = targetResp.errHC;
errLC = targetResp.errLC;
stdHC = sqrt( sum( errHC.^2 ) / numel(errHC) );
stdLC = sqrt( sum( errLC.^2 ) / numel(errLC) );

sigma_m = x(1);
Cc = x(2);

sigma_e = params(1);
sigma_i = params(2);
varGain = params(3);

sampleCnt = 100;
shapeParam = 1./varGain;
scaleParam = varGain;
gammaSamples = gaminv(linspace(1/sampleCnt, 1 - 1/sampleCnt, sampleCnt), ...
    shapeParam, scaleParam);
sigma_theta = sqrt( sigma_e^2 + (sigma_i*gammaSamples).^2 );

% For each value of sigma_theta, find the probability of high
% confidence and low confidence. These probabilitites are
% calculated by taking the inverse of the log normal distribution
% (which basically is a distribution of confidence variable).
muLogN    = - log((sigma_theta.^2) ./ sqrt(sigma_m.^2 + sigma_theta.^2));
sigmaLogN = sqrt(log((sigma_m.^2)./(sigma_theta.^2) + 1));

x1 = repmat(Cc, [1 sampleCnt]);
cdf_vals  = logncdf(x1, muLogN, sigmaLogN);
    
pHC = 1 - cdf_vals;
pLC = cdf_vals;

exp_sigma_theta_HC = sum( pHC.*sigma_theta ) / sum(pHC);
exp_sigma_theta_LC = sum( pLC.*sigma_theta ) / sum(pLC);

errorVal = (stdHC - exp_sigma_theta_HC)^2 + (stdLC - exp_sigma_theta_LC)^2;

% disp(errorVal)

end

