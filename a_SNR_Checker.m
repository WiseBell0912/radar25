clear; close all; clc;

%% Load Data
load("new.mat");
s_Date = result_date;
s_SNR = double(result_SNR);
load("ADCP.mat");
load("BOUY.mat");
gangleng = readtable("gangleng.csv");
gangleng_date = gangleng.x___1;
gangleng_precipitation = gangleng.x______mm_;
gangleng_precipitation_yn = 1 * (gangleng_precipitation > 0);


%% Operation of Data
d_start = min(s_Date);
d_end = max(s_Date);

% cut rain
gangleng_date_cut = gangleng_date(d_start <= gangleng_date & gangleng_date <= d_end);
gangleng_precipitation_yn_cut = gangleng_precipitation_yn(d_start <= gangleng_date & gangleng_date <= d_end);
[iscommon, idx] = ismember(gangleng_date_cut, s_Date);
metch_date = gangleng_date_cut(iscommon);
metch_idx = idx(iscommon);
metch_yn = gangleng_precipitation_yn_cut(iscommon);
%

b_Date_cut = b_Date(d_start <= b_Date & b_Date <= d_end);
b_Hs_cut = b_Hs(d_start <= b_Date & b_Date <= d_end);

b_Hs_interpol = interp1(b_Date_cut, b_Hs_cut, s_Date, "linear");

s_Sqrt_SNR = sqrt(s_SNR);

valid_indices = isfinite(b_Hs_interpol) & isfinite(s_Sqrt_SNR);
s_Date = s_Date(valid_indices);
b_Hs_interpol = b_Hs_interpol(valid_indices);
s_Sqrt_SNR = s_Sqrt_SNR(valid_indices);

modelfun = @(x, sqrt_SNR) x(1) + x(2) * sqrt_SNR;
initial_guess = [1, 1];
options = optimoptions('lsqcurvefit', 'Display', 'off');
x = lsqcurvefit(modelfun, initial_guess, s_Sqrt_SNR, b_Hs_interpol, [], [], options);

s_Hs = x(1) + x(2) .* s_Sqrt_SNR;

%% 플롯(전체)
figure(1)
hold on;
plot(gangleng_date, gangleng_precipitation_yn, "DisplayName", "Precipitation", "Marker",".", "Color","k","LineStyle","none");
plot(s_Date, movmean(b_Hs_interpol, 20), "DisplayName", "Bouy", "Color", "b");
plot(s_Date, movmean(s_Hs, 20), "DisplayName", "Radar", "Color", "r");
%plot(a_Date02, a_Hs02, "DisplayName", "ADCP");
hold off;
% option
xlabel("Time"); ylabel("Hs[m]");
legend;
title("Significant Wave Height", "FontSize", 15);
xlim([min(s_Date), max(s_Date)]);
ylim([0, 5]);
set(gcf, 'Position', [0, 0, 1820, 980]);

%% Error
error = b_Hs_interpol - s_Hs;

% figure(2);
% hold on;
% plot(s_Date, error);
% plot([min(s_Date), max(s_Date)], [0, 0], '--');
% hold off;
% % option
% xlabel("Time"); ylabel("Error[m]");
% title("Error", "FontSize", 15);
% xlim([min(s_Date), max(s_Date)]);
% set(gcf, 'Position', [0, 0, 1820, 980]);

% Check best 10
abs_error = abs(error);
[sort_abs_error_val, sort_abs_error_idx] = sort(abs_error, 'descend');

bestidx = 1:200;
bestidx = bestidx';
best10_error_date = s_Date(sort_abs_error_idx(bestidx));
best10_error_sqrt_SNR = s_Sqrt_SNR(sort_abs_error_idx(bestidx));
best10_error_radar_Hs = s_Hs(sort_abs_error_idx(bestidx));
best10_error_bouy_Hs = b_Hs_interpol(sort_abs_error_idx(bestidx));

%% R2
y_mean = mean(b_Hs_interpol);
SS_tot = sum((b_Hs_interpol - y_mean).^2);
SS_res = sum((b_Hs_interpol - s_Hs).^2);
R_squared = 1 - (SS_res / SS_tot);

figure(3);
hold on;
plot(s_Hs, b_Hs_interpol, '.', 'Color', 'k');
plot(0:0.5:5, 0:0.5:5, '--', 'Color', 'b');
plot(0:0.5:5, 1:0.5:6, '--', 'Color', 'b');
plot(0:0.5:5, -1:0.5:4, '--', 'Color', 'b');
plot(best10_error_radar_Hs, best10_error_bouy_Hs, 'ro', 'MarkerSize', 6)
hold off;
% option
xlabel("Radar Hs[m]"); ylabel("Bouy Hs[m]");
title(['R^2 = ', num2str(R_squared)], "FontSize", 15);
xlim([0, 5]);
ylim([0, 5]);
set(gcf, 'Position', [0, 0, 980, 980]);

%% SNR-Hs
figure(4);
hold on;
plot(s_Sqrt_SNR, b_Hs_interpol, 'k.');
plot((-5:0.1:10), x(1)+x(2)*(-5:0.1:10), 'r--', LineWidth=1);
hold off;
% option
xlabel('sqrt(SNR)'); ylabel("Radar Hs[m]");
title('SNR-Hs', "FontSize", 15);
xlim([min(s_Sqrt_SNR), max(s_Sqrt_SNR)]);
ylim([min(b_Hs_interpol), max(b_Hs_interpol)]);
set(gcf, 'Position', [0, 0, 1820, 980]);