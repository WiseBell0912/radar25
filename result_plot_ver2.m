clear; close all; clc;

%% 불러오기
% ADCP
load("ADCP.mat");
a_Date = a_Date02;
a_Hs = a_Hs02;

clear a_Date01 a_Hs01 a_Pdir01 a_Tp01 a_Date02 a_Hs02 a_Pdir02 a_Tp02

% Bouy
load("new_Bouy_no_NAN_ZERO.mat");

% Radar wave
load("snr_y1910_wave_0327.mat");
r_Date = Date;
r_wave_SNR = SNR;
r_wave_Ux = Ux;
r_wave_Uy = Uy;

load("snr_y1911_wave_0327.mat");
r_Date = [r_Date ; Date];
r_wave_SNR = [r_wave_SNR ; SNR];
r_wave_Ux = [r_wave_Ux ; Ux];
r_wave_Uy = [r_wave_Uy ; Uy];

load("snr_y1912_wave_0327.mat");
r_Date = [r_Date ; Date];
r_wave_SNR = [r_wave_SNR ; SNR];
r_wave_Ux = [r_wave_Ux ; Ux];
r_wave_Uy = [r_wave_Uy ; Uy];

r_wave_SNR = sqrt(r_wave_SNR);

% Radar surf
load("snr_y1910_surf_0327.mat");
r_Date = Date;
r_surf_SNR = SNR;
r_surf_Ux = Ux;
r_surf_Uy = Uy;

load("snr_y1911_surf_0327.mat");
r_Date = [r_Date ; Date];
r_surf_SNR = [r_surf_SNR ; SNR];
r_surf_Ux = [r_surf_Ux ; Ux];
r_surf_Uy = [r_surf_Uy ; Uy];

load("snr_y1912_surf_0327.mat");
r_Date = [r_Date ; Date];
r_surf_SNR = [r_surf_SNR ; SNR];
r_surf_Ux = [r_surf_Ux ; Ux];
r_surf_Uy = [r_surf_Uy ; Uy];

r_surf_SNR = sqrt(r_surf_SNR);

clear Date SNR Ux Uy r_wave_Ux r_wave_Uy r_surf_Ux r_surf_Uy

%% 처리
mask_r2b = ismember(r_Date, b_Date);

r2b_Date = r_Date(mask_r2b);
r2b_wave_SNR = r_wave_SNR(mask_r2b);
r2b_surf_SNR = r_surf_SNR(mask_r2b);

mask_b2r = ismember(b_Date, r_Date);

b2r_Date = b_Date(mask_b2r);
b2r_Hs = b_Hs(mask_b2r);
b2r_Wind = b_Wind(mask_b2r);

%% 개발
modelfun = @(x, r2b_wave_SNR) x(1) + x(2) * r2b_wave_SNR;
initial_guess = [1, 1];
options = optimoptions('lsqcurvefit', 'Display', 'off');
x = lsqcurvefit(modelfun, initial_guess, r2b_wave_SNR, b2r_Hs, [], [], options);

r2b_wave_Hs = x(1) + x(2) .* r2b_wave_SNR;

modelfun = @(y, r2b_surf_SNR) y(1) + y(2) * r2b_surf_SNR;
initial_guess = [1, 1];
options = optimoptions('lsqcurvefit', 'Display', 'off');
y = lsqcurvefit(modelfun, initial_guess, r2b_surf_SNR, b2r_Hs, [], [], options);

r2b_surf_Hs = y(1) + y(2) .* r2b_surf_SNR;

%% 확인 Graph
for i = 1 : 1
    a = figure(1);
    hold on;
    % yyaxis left
    plot(b2r_Date, movmean(b2r_Hs, 12), 'Color', [0, 0, 0, 0.8], 'LineStyle', '-', 'LineWidth',0.3);
    plot(r2b_Date, movmean(r2b_wave_Hs, 12), 'Color', [1, 0, 0, 0.8], 'LineStyle', '-', 'LineWidth',0.3);
    plot(r2b_Date, movmean(r2b_surf_Hs, 12), 'Color', [0, 0, 1, 0.8], 'LineStyle', '-', 'LineWidth',0.3);
    ylabel("Hs [m]");
    % yyaxis right
    % plot(b2r_Date, movmean(b2r_Wind, 3), 'Color', [0.9, 0.2, 0.9, 0.7], 'LineStyle', '-.');
    % ylabel("Wind Velocity [m/s]");
    hold off;

    set(gcf, 'Position', [0, 0, 1820, 980]);
    xlim([datetime(2019, 10, 1), datetime(2019, 12, 30) + days(2)]);
    title("Significant Wave Height", "FontSize", 15);
    xlabel("Date [mm dd]");
    legend('Bouy', 'Wave', 'Surf', 'Wind Velocity');

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
y_mean = mean(b2r_Hs);
SS_tot = sum((b2r_Hs - y_mean).^2);
SS_res = sum((b2r_Hs - r2b_wave_Hs).^2);
R_squared = 1 - (SS_res / SS_tot);

figure(2);
hold on;
plot(b2r_Hs, r2b_wave_Hs, '.', 'MarkerSize', 3);
plot(0:0.5:5, 0:0.5:5, '--', 'Color', 'r', 'LineWidth', 1);
plot(0:0.5:5, 1:0.5:6, '--', 'Color', 'k', 'LineWidth', 1);
plot(0:0.5:5, -1:0.5:4, '--', 'Color', 'k', 'LineWidth', 1);
hold off;

set(gcf, 'Position', [0, 0, 980, 980]);
xlim([0, 5]);
ylim([0, 5]);
title(['R^2 = ', num2str(R_squared)], "FontSize", 15);
subtitle('Wave Zone');
xlabel("Bouy Hs [m]"); ylabel("Rader Hs [m]");

% surf
y_mean = mean(b2r_Hs);
SS_tot = sum((b2r_Hs - y_mean).^2);
SS_res = sum((b2r_Hs - r2b_surf_Hs).^2);
R_squared = 1 - (SS_res / SS_tot);

figure(3);
hold on;
plot(b2r_Hs, r2b_surf_Hs, '.', 'MarkerSize', 3);
plot(0:0.5:5, 0:0.5:5, '--', 'Color', 'r', 'LineWidth', 1);
plot(0:0.5:5, 1:0.5:6, '--', 'Color', 'k', 'LineWidth', 1);
plot(0:0.5:5, -1:0.5:4, '--', 'Color', 'k', 'LineWidth', 1);
hold off;

set(gcf, 'Position', [0, 0, 980, 980]);
xlim([0, 5]);
ylim([0, 5]);
title(['R^2 = ', num2str(R_squared)], "FontSize", 15);
subtitle('Surf Zone');
xlabel("Bouy Hs [m]"); ylabel("Rader Hs [m]");

%% 확인 SNR-Hs
% wave
figure(4);
hold on;
plot(r2b_wave_SNR, b2r_Hs, '.', 'MarkerSize', 3);
plot(0:1:10, x(1)+x(2)*(0:1:10), 'r--');
plot(0:1:10, x(1)+x(2)*(0:1:10)+1, 'k--');
plot(0:1:10, x(1)+x(2)*(0:1:10)-1, 'k--');
hold off;

set(gcf, 'Position', [0, 0, 1820, 980]);
xlim([0.4, 2]);
ylim([0, 5]);
title("SNR-Hs", "FontSize", 15);
subtitle('Wave Zone');
xlabel("SNR"); ylabel("Hs [m]");

% surf
figure(5);
hold on;
plot(r2b_surf_SNR, b2r_Hs, '.', 'MarkerSize', 3);
plot(0:1:10, y(1)+y(2)*(0:1:10), 'r--');
plot(0:1:10, y(1)+y(2)*(0:1:10)+1, 'k--');
plot(0:1:10, y(1)+y(2)*(0:1:10)-1, 'k--');
hold off;

set(gcf, 'Position', [0, 0, 1820, 980]);
xlim([0.4, 2]);
ylim([0, 5]);
title("SNR-Hs", "FontSize", 15);
subtitle('Surf Zone');
xlabel("SNR"); ylabel("Hs [m]");