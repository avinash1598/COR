% close all
clear all

% %% Analytical reward function
% maxTolerableError = 20; % In degrees
% sigmaHC = 3.25;         % HC reward function std deviation
% valLC   = 0.3;          % LC constant reward
% x = linspace(-20, 20, 100);
% 
% mu = 0;
% sigma = sigmaHC;
% y1 = exp( - (abs(x)).^2 / (2 * sigma^2) );
% y2 = ones(1, numel(x))*valLC;
% 
% x1 = sqrt( - log(0.3)* 2 * sigma^2 );
% x2 = -x1;
% 
% figure
% hold on
% plot(x, y1, DisplayName="HC", LineWidth=2)
% plot(x, y2, DisplayName="LC", LineWidth=2)
% xlabel("Perceptual error")
% ylabel("Reward")
% legend
% 
% xline(x1, LineWidth=2, LineStyle="--", HandleVisibility="off")
% xline(x2, LineWidth=2, LineStyle="--", HandleVisibility="off")
% 
% hold off
% 
% x1, x2
% %% Analytical reward function
% maxTolerableError = 20; % In degrees
% x = linspace(-20, 20, 201);
% 
% y1 = 7; %20
% x1 = 5.97;
% sigmaErrSD_C0_015 = 8.85;
% prob = cdf('Normal', x1, 0, sigmaErrSD_C0_015) - cdf('Normal', -x1, 0, sigmaErrSD_C0_015);
% disp(prob)
% 
% c1 = 1;
% c2 = 0.21;
% 
% m1 = 1/y1;
% m2 = m1 - (c1 - c2)/x1;
% 
% yHC = - m1 * sign(x) .* x + c1; yHC(abs(x) > y1) = 0;
% yLC = - m2 * sign(x) .* x + c2; yLC(abs(x) > maxTolerableError) = 0;
% 
% x2 = -x1;
% 
% figure
% hold on
% plot(x, yHC, DisplayName="HC", LineWidth=2)
% plot(x, yLC, DisplayName="LC", LineWidth=2)
% xlabel("Perceptual error")
% ylabel("Reward")
% legend
% 
% xline(x1, LineWidth=2, LineStyle="--", HandleVisibility="off")
% xline(x2, LineWidth=2, LineStyle="--", HandleVisibility="off")
% 
% hold off



%% Analytical reward function
maxTolerableError = 25; % In degrees
x = linspace(-30, 30, 201);

y1 = 11; %20
x1 = 10;

c1 = 5;
c2 = 0.3;

m1 = c1/y1;
m2 = c2 / maxTolerableError; %m1 - (c1 - c2)/x1;

yHC = - m1 * sign(x) .* x + c1; yHC(abs(x) > y1) = 0.3*yHC(abs(x) > y1);
yLC = - m2 * sign(x) .* x + c2; yLC(abs(x) > maxTolerableError) = 0;


figure
hold on
plot(x, yHC, 'DisplayName', "HC", 'LineWidth', 2)
plot(x, yLC, 'DisplayName', "LC", 'LineWidth', 2)
xlabel("Perceptual error")
ylabel("Reward")
legend

yline(0, 'LineStyle', "--")
hold off


