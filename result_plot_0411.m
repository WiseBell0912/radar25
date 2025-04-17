clear; close all; clc;

%% 불러오기
% ADCP
load("ADCP.mat");
a_Date = a_Date02;
a_Hs = a_Hs02;

clear a_Date01 a_Hs01 a_Pdir01 a_Tp01 a_Date02 a_Hs02 a_Pdir02 a_Tp02

% Bouy
load("BOUY.mat");
b_Hs = b_SignificantWaveHeight;
b_Wind = b_WindSpeed;
b_Pdir = b_WaveDirection;
b_Tp = b_MaximumWavePeriod;

b_Date = b_Date(~isnan(b_Hs));
b_Hs = b_Hs(~isnan(b_Hs));
b_Wind = b_Wind(~isnan(b_Hs));
b_Pdir = b_Pdir(~isnan(b_Hs));
b_Tp = b_Tp(~isnan(b_Hs));

b_Hs(b_Hs == 0) = NaN;
b_Hs = fillmissing(b_Hs, 'linear');

load("snr_y1910_0416.mat");
r_WaveE = WaveE;
r_wave_Uy = wave_Uy;
r_wave_Ux = wave_Ux;
r_wave_Tp = wave_Tp;
r_wave_SNR = wave_SNR;
r_wave_Pdir = wave_Pdir;
%r_wave_SpectrumMax = wave_SpectrumMax;
r_SurfE = SurfE;
r_surf_Uy = surf_Uy;
r_surf_Ux = surf_Ux;
r_surf_Tp = surf_Tp;
r_surf_SNR = surf_SNR;
r_surf_Pdir = surf_Pdir;
%r_surf_SpectrumMax = surf_SpectrumMax;
r_LandE = LandE;
r_Date = Date;

% load("snr_y1911_0412.mat");
% r_WaveE = [r_WaveE; WaveE];
% r_wave_Uy = [r_wave_Uy; wave_Uy];
% r_wave_Ux = [r_wave_Ux; wave_Ux];
% r_wave_Tp = [r_wave_Tp; wave_Tp];
% r_wave_SNR = [r_wave_SNR; wave_SNR];
% r_wave_Pdir = [r_wave_Pdir; wave_Pdir];
% r_wave_SpectrumMax = [r_wave_SpectrumMax; wave_SpectrumMax];
% r_SurfE = [r_SurfE; SurfE];
% r_surf_Uy = [r_surf_Uy; surf_Uy];
% r_surf_Ux = [r_surf_Ux; surf_Ux];
% r_surf_Tp = [r_surf_Tp; surf_Tp];
% r_surf_SNR = [r_surf_SNR; surf_SNR];
% r_surf_Pdir = [r_surf_Pdir; surf_Pdir];
% r_surf_SpectrumMax = [r_surf_SpectrumMax; surf_SpectrumMax];
% r_LandE = [r_LandE; LandE];
% r_Date = [r_Date; Date];
% 
% load("snr_y1912_0412.mat");
% r_WaveE = [r_WaveE; WaveE];
% r_wave_Uy = [r_wave_Uy; wave_Uy];
% r_wave_Ux = [r_wave_Ux; wave_Ux];
% r_wave_Tp = [r_wave_Tp; wave_Tp];
% r_wave_SNR = [r_wave_SNR; wave_SNR];
% r_wave_Pdir = [r_wave_Pdir; wave_Pdir];
% r_wave_SpectrumMax = [r_wave_SpectrumMax; wave_SpectrumMax];
% r_SurfE = [r_SurfE; SurfE];
% r_surf_Uy = [r_surf_Uy; surf_Uy];
% r_surf_Ux = [r_surf_Ux; surf_Ux];
% r_surf_Tp = [r_surf_Tp; surf_Tp];
% r_surf_SNR = [r_surf_SNR; surf_SNR];
% r_surf_Pdir = [r_surf_Pdir; surf_Pdir];
% r_surf_SpectrumMax = [r_surf_SpectrumMax; surf_SpectrumMax];
% r_LandE = [r_LandE; LandE];
% r_Date = [r_Date; Date];

r_wave_SNR = sqrt(r_wave_SNR);
r_surf_SNR = sqrt(r_surf_SNR);

r_wave_U = sqrt(r_wave_Ux.^2 + r_wave_Uy.^2);
r_surf_U = sqrt(r_surf_Ux.^2 + r_surf_Uy.^2);

r_wave_filter = (r_LandE < mean(r_LandE) + 2 * std(r_LandE));% & (r_wave_SpectrumMax > mean(r_wave_SpectrumMax) - 0.7 * std(r_wave_SpectrumMax));
r_surf_filter = (r_LandE < mean(r_LandE) + 2 * std(r_LandE));% & (r_surf_SpectrumMax > mean(r_surf_SpectrumMax) - 0.7 * std(r_surf_SpectrumMax));

% 예전 데이터
load("./0327/snr_y1910_wave_0327.mat");
r_wave_SNR = SNR;

load("./0327/snr_y1911_wave_0327.mat");
r_wave_SNR = [r_wave_SNR ; SNR];

load("./0327/snr_y1912_wave_0327.mat");
r_wave_SNR = [r_wave_SNR ; SNR];

r_wave_SNR = sqrt(r_wave_SNR);

load("./0327/snr_y1910_surf_0327.mat");
r_surf_SNR = SNR;

load("./0327/snr_y1911_surf_0327.mat");
r_surf_SNR = [r_surf_SNR ; SNR];

load("./0327/snr_y1912_surf_0327.mat");
r_surf_SNR = [r_surf_SNR ; SNR];

r_surf_SNR = sqrt(r_surf_SNR);

clear Date LandE surf_Pdir surf_SNR surf_SpectrumMax surf_Tp surf_Ux surf_Uy SurfE wave_Pdir wave_SNR wave_SpectrumMax wave_Tp wave_Ux wave_Uy WaveE

%% 처리
mask_r2b = ismember(r_Date, b_Date);

rr_Date = r_Date(mask_r2b);
rr_wave_SNR = r_wave_SNR(mask_r2b);
rr_surf_SNR = r_surf_SNR(mask_r2b);
rr_LandE = r_LandE((mask_r2b));
%rr_wave_SpectrumMax = r_wave_SpectrumMax(mask_r2b);
%rr_surf_SpectrumMax = r_surf_SpectrumMax(mask_r2b);
rr_wave_filter = r_wave_filter(mask_r2b);
rr_surf_filter = r_surf_filter(mask_r2b);

mask_b2r = ismember(b_Date, r_Date);

bb_Date = b_Date(mask_b2r);
bb_Hs = b_Hs(mask_b2r);
bb_Wind = b_Wind(mask_b2r);

%% 개발
modelfun = @(x, SNR) x(1) + x(2) * SNR;
initial_guess = [1, 1];
options = optimoptions('lsqcurvefit', 'Display', 'off');
x = lsqcurvefit(modelfun, initial_guess, rr_wave_SNR(rr_wave_filter), bb_Hs(rr_wave_filter), [], [], options);

rr_wave_Hs = x(1) + x(2) .* rr_wave_SNR(rr_wave_filter);

modelfun = @(y, SNR) y(1) + y(2) * SNR;
initial_guess = [1, 1];
options = optimoptions('lsqcurvefit', 'Display', 'off');
y = lsqcurvefit(modelfun, initial_guess, rr_surf_SNR(rr_surf_filter), bb_Hs(rr_surf_filter), [], [], options);

rr_surf_Hs = y(1) + y(2) .* rr_surf_SNR(rr_surf_filter);

%% 확인 Graph
for i = 1 : 1
    a = figure(1);
    hold on;
    % yyaxis left
    plot(bb_Date, movmean(bb_Hs, 1), 'Color', [0, 0, 0, 0.8], 'LineStyle', '-', 'LineWidth',0.3);
    plot(rr_Date(rr_wave_filter), movmean(rr_wave_Hs, 1), 'Color', [1, 0, 0, 0.8], 'LineStyle', '-', 'LineWidth',0.3);
    plot(rr_Date(rr_surf_filter), movmean(rr_surf_Hs, 1), 'Color', [0, 0, 1, 0.8], 'LineStyle', '-', 'LineWidth',0.3);
    ylabel("Hs [m]");
    % yyaxis right
    % plot(b2r_Date, movmean(b2r_Wind, 3), 'Color', [0.9, 0.2, 0.9, 0.7], 'LineStyle', '-.');
    % ylabel("Wind Velocity [m/s]");
    hold off;

    set(gcf, 'Position', [0, 0, 1820, 980]);
    xlim([min(rr_Date), max(rr_Date)]);
    title("Significant Wave Height", "FontSize", 15);
    xlabel("Date [mm dd]");
    legend('Bouy', 'Wave', 'Surf');

    % if i < 10
    %     saveas(gcf, ['./Image0324/w4_2019120', num2str(i), '.fig']);
    %     saveas(gcf, ['./Image0324/w4_2019120', num2str(i), '.png']);
    % else
    %     saveas(gcf, ['./Image0324/w4_201912', num2str(i), '.fig']);
    %     saveas(gcf, ['./Image0324/w4_201912', num2str(i), '.png']);
    % end

    % close(a);
end

%% 확인 R2
% wave
y_mean = mean(bb_Hs(rr_wave_filter));
SS_tot = sum((bb_Hs(rr_wave_filter) - y_mean).^2);
SS_res = sum((bb_Hs(rr_wave_filter) - rr_wave_Hs).^2);
R_squared = 1 - (SS_res / SS_tot);

figure(2);
hold on;
plot(bb_Hs(rr_wave_filter), rr_wave_Hs, '.', 'MarkerSize', 3);
plot(0:0.5:5, 0:0.5:5, '--', 'Color', 'r', 'LineWidth', 1);
plot(0:0.5:5, 1:0.5:6, '--', 'Color', 'k', 'LineWidth', 1);
plot(0:0.5:5, -1:0.5:4, '--', 'Color', 'k', 'LineWidth', 1);
hold off;

set(gcf, 'Position', [0, 0, 980, 980]);
xlim([0, max(bb_Hs)]);
ylim([0, max(bb_Hs)]);
title(['R^2 = ', num2str(R_squared)], "FontSize", 15);
subtitle('Wave Zone');
xlabel("Bouy Hs [m]"); ylabel("Rader Hs [m]");

% surf
y_mean = mean(bb_Hs(rr_surf_filter));
SS_tot = sum((bb_Hs(rr_surf_filter) - y_mean).^2);
SS_res = sum((bb_Hs(rr_surf_filter) - rr_surf_Hs).^2);
R_squared = 1 - (SS_res / SS_tot);

figure(3);
hold on;
plot(bb_Hs(rr_surf_filter), rr_surf_Hs, '.', 'MarkerSize', 3);
plot(0:0.5:5, 0:0.5:5, '--', 'Color', 'r', 'LineWidth', 1);
plot(0:0.5:5, 1:0.5:6, '--', 'Color', 'k', 'LineWidth', 1);
plot(0:0.5:5, -1:0.5:4, '--', 'Color', 'k', 'LineWidth', 1);
hold off;

set(gcf, 'Position', [0, 0, 980, 980]);
xlim([0, max(bb_Hs)]);
ylim([0, max(bb_Hs)]);
title(['R^2 = ', num2str(R_squared)], "FontSize", 15);
subtitle('Surf Zone');
xlabel("Bouy Hs [m]"); ylabel("Rader Hs [m]");

%% 확인 SNR-Hs
% wave
figure(4);
hold on;
plot(rr_wave_SNR, bb_Hs, '.', 'MarkerSize', 3);
plot(0:1:10, x(1)+x(2)*(0:1:10), 'r--');
plot(0:1:10, x(1)+x(2)*(0:1:10)+1, 'k--');
plot(0:1:10, x(1)+x(2)*(0:1:10)-1, 'k--');
hold off;

set(gcf, 'Position', [0, 0, 1820, 980]);
xlim([min([rr_wave_SNR ; rr_surf_SNR]), max([rr_wave_SNR ; rr_surf_SNR])]);
ylim([0, max(bb_Hs)]);
title("SNR-Hs", "FontSize", 15);
subtitle('Wave Zone');
xlabel("SNR"); ylabel("Hs [m]");

% surf
figure(5);
hold on;
plot(rr_surf_SNR, bb_Hs, '.', 'MarkerSize', 3);
plot(0:1:10, y(1)+y(2)*(0:1:10), 'r--');
plot(0:1:10, y(1)+y(2)*(0:1:10)+1, 'k--');
plot(0:1:10, y(1)+y(2)*(0:1:10)-1, 'k--');
hold off;

set(gcf, 'Position', [0, 0, 1820, 980]);
xlim([min([rr_wave_SNR ; rr_surf_SNR]), max([rr_wave_SNR ; rr_surf_SNR])]);
ylim([0, max(bb_Hs)]);
title("SNR-Hs", "FontSize", 15);
subtitle('Surf Zone');
xlabel("SNR"); ylabel("Hs [m]");