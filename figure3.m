clear; close all; clc;

load("ADCP.mat");
load("BOUY.mat");

figure;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);

hold on;
plot(b_Date, movmean(b_SignificantWaveHeight, 9), 'LineWidth', 1);
plot(a_Date02, movmean(a_Hs02, 3), 'LineWidth', 1);
hold off;

title('Significant wave height time series');
legend('Bouy', 'ADCP');
xlim([min(a_Date02), max(a_Date02)]);