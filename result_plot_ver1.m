clear; close all; clc;

%% 수집
load("BOUY.mat");

load("snr_y1910_wave.mat");
r_Date = Date;
r_SNR = SNR;
r_Ux = Ux;
r_Uy = Uy;

% load("snr_y1911_wave.mat");
% r_Date = [r_Date ; Date];
% r_SNR = [r_SNR ; SNR];
% r_Ux = [r_Ux ; Ux];
% r_Uy = [r_Uy ; Uy];

r_SNR = sqrt(r_SNR);

clear Date SNR Ux Uy r_Ux r_Uy;

%% 가공
r_start = min(r_Date);
r_finish = max(r_Date);

cut_idx = (r_start <= b_Date) & (b_Date <= r_finish) & (0 < b_Hs);
b_Date_cut = b_Date(cut_idx);
b_Hs_cut = b_Hs(cut_idx);

b_Hs_interp = interp1(b_Date_cut, b_Hs_cut, r_Date, "linear");

clear r_start r_finish cut_idx b_Date b_Hs b_Date_cut b_Hs_cut;

%% 처리
valid_idx = ~isnan(b_Hs_interp) & ~isnan(r_SNR);
b_Hs_interp = b_Hs_interp(valid_idx);
r_SNR = r_SNR(valid_idx);
r_Date = r_Date(valid_idx);

clear valid_idx

%% 개발
modelfun = @(x, r_SNR) x(1) + x(2) * r_SNR;
initial_guess = [1, 1];
options = optimoptions('lsqcurvefit', 'Display', 'off');
x = lsqcurvefit(modelfun, initial_guess, r_SNR, b_Hs_interp, [], [], options);

r_Hs = x(1) + x(2) .* r_SNR;

clear modelfun options initial_guess

%% 확인 graph
figure(1);
hold on;
plot(r_Date, movmean(b_Hs_interp, 1), 'Color', [0, 0, 1, 0.3]);
plot(r_Date, movmean(r_Hs, 1), 'Color', [1, 0, 0, 0.3]);
hold off;

set(gcf, 'Position', [0, 0, 1820, 980]);
title("Significant Wave Height", "FontSize", 15);
xlabel("Date [mm dd]"); ylabel("Hs [m]");
legend('Bouy', 'Rader');

%% 확인 r2
y_mean = mean(b_Hs_interp);
SS_tot = sum((b_Hs_interp - y_mean).^2);
SS_res = sum((b_Hs_interp - r_Hs).^2);
R_squared = 1 - (SS_res / SS_tot);

figure(2);
hold on;
plot(b_Hs_interp, r_Hs, '.', 'MarkerSize', 3);
plot(0:0.5:5, 0:0.5:5, '--', 'Color', 'r', 'LineWidth', 1);
plot(0:0.5:5, 1:0.5:6, '--', 'Color', 'k', 'LineWidth', 1);
plot(0:0.5:5, -1:0.5:4, '--', 'Color', 'k', 'LineWidth', 1);
hold off;

set(gcf, 'Position', [0, 0, 980, 980]);
xlim([0, 5]);
ylim([0, 5]);
title(['R^2 = ', num2str(R_squared)], "FontSize", 15);
xlabel("Bouy Hs [m]"); ylabel("Rader Hs [m]");

clear R_squared SS_res SS_tot y_mean

%% 확인 SNR-Hs
figure(3);
hold on;
plot(r_SNR, b_Hs_interp, '.', 'MarkerSize', 3);
plot(0:1:10, x(1)+x(2)*(0:1:10), 'r--');
plot(0:1:10, x(1)+x(2)*(0:1:10)+1, 'k--');
plot(0:1:10, x(1)+x(2)*(0:1:10)-1, 'k--');
hold off;

set(gcf, 'Position', [0, 0, 1820, 980]);
xlim([0.4, 2]);
ylim([0, 5]);
title("SNR-Hs", "FontSize", 15);
xlabel("SNR"); ylabel("Hs [m]");

%% 아웃라이어
E = b_Hs_interp - r_Hs;
abs_E = abs(E);
%idx_E_over_1m = abs_E > 1;
idx_E_over_1m = (abs_E > 0.5);
date_E_over_1m = r_Date(idx_E_over_1m);

y_mean = mean(b_Hs_interp);
SS_tot = sum((b_Hs_interp - y_mean).^2);
SS_res = sum((b_Hs_interp - r_Hs).^2);
R_squared = 1 - (SS_res / SS_tot);

figure(5);
hold on;
plot(b_Hs_interp, r_Hs, '.', 'MarkerSize', 3);
plot(0:0.5:5, 0:0.5:5, '--', 'Color', 'r', 'LineWidth', 1);
plot(0:0.5:5, 1:0.5:6, '--', 'Color', 'k', 'LineWidth', 1);
plot(0:0.5:5, -1:0.5:4, '--', 'Color', 'k', 'LineWidth', 1);
for i = 1 : length(idx_E_over_1m)
    if idx_E_over_1m(i) == 1
        plot(b_Hs_interp(i), r_Hs(i), 'o', 'Color', 'm')
    end
end
hold off;

set(gcf, 'Position', [0, 0, 980, 980]);
xlim([0, 5]);
ylim([0, 5]);
title(['R^2 = ', num2str(R_squared)], "FontSize", 15);
xlabel("Bouy Hs [m]"); ylabel("Rader Hs [m]");

figure(6);
hold on;
plot(r_SNR, b_Hs_interp, '.', 'MarkerSize', 3);
plot(0:1:10, x(1)+x(2)*(0:1:10), 'r--');
plot(0:1:10, x(1)+x(2)*(0:1:10)+1, 'k--');
plot(0:1:10, x(1)+x(2)*(0:1:10)-1, 'k--');
for i = 1 : length(idx_E_over_1m)
    if idx_E_over_1m(i) == 1
        plot(r_SNR(i), b_Hs_interp(i), 'o', 'Color', 'm')
    end
end
hold off;

set(gcf, 'Position', [0, 0, 1820, 980]);
xlim([0.4, 2]);
ylim([0, 5]);
title("SNR-Hs", "FontSize", 15);
xlabel("SNR"); ylabel("Hs [m]");

%% 이미지 분류
close all;

d = dir('./img/*.png');

for i = 1 : length(d)
    if idx_E_over_1m(i) == 1


        figure(7);
        set(gcf, 'Position', [0, 0, 1820, 500]);
        t = tiledlayout(2,2);
        maintitle = [datestr(r_Date(i), 'yyyy/mm/dd HH:MM'), '  &  Bouy = ', num2str(round(b_Hs_interp(i), 2)), 'm  &  ', 'Radar = ', num2str(round(r_Hs(i), 2)), 'm'];
        title(t, maintitle, "FontSize", 15);


        nexttile;
        hold on;
        plot(r_Date, movmean(r_Hs, 1));
        plot(r_Date, movmean(b_Hs_interp, 1));
        highlight_start = r_Date(i-1); % 시작 시간
        highlight_end = r_Date(i+1); % 끝나는 시간
        y_limits = ylim; % y축 범위 가져오기
        fill([highlight_start highlight_end highlight_end highlight_start], ...
            [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
            [1 1 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', 'm'); % 노란색 형광펜 효과
        hold off;

        if i < 60
            xlim([r_Date(1), r_Date(i+60)]);
        elseif i > length(d) - 60
            xlim([r_Date(i-60), r_Date(length(d))]);
        else
            xlim([r_Date(i-60), r_Date(i+60)]);
        end
        title("Significant Wave Height");
        xlabel("Date [mm dd]"); ylabel("Hs [m]");
        legend('Rader', 'Bouy');


        nexttile(2, [2, 1]);
        hold on;
        plot(b_Hs_interp, r_Hs, '.', 'MarkerSize', 3);
        plot(0:0.5:5, 0:0.5:5, '--', 'Color', 'r', 'LineWidth', 1);
        plot(0:0.5:5, 1:0.5:6, '--', 'Color', 'k', 'LineWidth', 1);
        plot(0:0.5:5, -1:0.5:4, '--', 'Color', 'k', 'LineWidth', 1);
        plot(b_Hs_interp(i), r_Hs(i), 'o', 'Color', 'm', 'LineWidth', 2)
        hold off;

        xlim([0, 5]);
        ylim([0, 5]);
        title(['R^2 = ', num2str(R_squared)]);
        xlabel("Bouy Hs [m]"); ylabel("Rader Hs [m]");


        nexttile;
        hold on;
        plot(r_SNR, b_Hs_interp, '.', 'MarkerSize', 3);
        plot(0:1:10, x(1)+x(2)*(0:1:10), 'r--');
        plot(0:1:10, x(1)+x(2)*(0:1:10)+1, 'k--');
        plot(0:1:10, x(1)+x(2)*(0:1:10)-1, 'k--');
        plot(r_SNR(i), b_Hs_interp(i), 'o', 'Color', 'm', 'LineWidth', 2)
        hold off;

        set(gcf, 'Position', [0, 0, 1820, 980]);
        xlim([0.4, 2]);
        ylim([0, 5]);
        title("SNR-Hs");
        xlabel("SNR"); ylabel("Hs [m]");


        savename = ['./error/', num2str(i), '_', 'Total_', datestr(r_Date(i), 'yyyymmdd_HHMM'), '.png'];
        saveas(t, savename);


        originalname = [d(i).folder, '/', d(i).name];
        copyname = ['./error/', num2str(i), '_', d(i).name];
        copyfile(originalname, copyname);

    end
end

%% 나이스
nice_over_1m = b_Hs_interp - r_Hs;
abs_nice_over_1m = abs(nice_over_1m);
idx_nice_over_1m = (abs_nice_over_1m < 0.05) & (b_Hs_interp > 2);
date_nice_over_1m = r_Date(idx_nice_over_1m);

y_mean = mean(b_Hs_interp);
SS_tot = sum((b_Hs_interp - y_mean).^2);
SS_res = sum((b_Hs_interp - r_Hs).^2);
R_squared = 1 - (SS_res / SS_tot);

figure(5);
hold on;
plot(b_Hs_interp, r_Hs, '.', 'MarkerSize', 3);
plot(0:0.5:5, 0:0.5:5, '--', 'Color', 'r', 'LineWidth', 1);
plot(0:0.5:5, 1:0.5:6, '--', 'Color', 'k', 'LineWidth', 1);
plot(0:0.5:5, -1:0.5:4, '--', 'Color', 'k', 'LineWidth', 1);
for i = 1 : length(idx_nice_over_1m)
    if idx_nice_over_1m(i) == 1
        plot(b_Hs_interp(i), r_Hs(i), 'o', 'Color', 'm')
    end
end
hold off;

set(gcf, 'Position', [0, 0, 980, 980]);
xlim([0, 5]);
ylim([0, 5]);
title(['R^2 = ', num2str(R_squared)], "FontSize", 15);
xlabel("Bouy Hs [m]"); ylabel("Rader Hs [m]");

figure(6);
hold on;
plot(r_SNR, b_Hs_interp, '.', 'MarkerSize', 3);
plot(0:1:10, x(1)+x(2)*(0:1:10), 'r--');
plot(0:1:10, x(1)+x(2)*(0:1:10)+1, 'k--');
plot(0:1:10, x(1)+x(2)*(0:1:10)-1, 'k--');
for i = 1 : length(idx_nice_over_1m)
    if idx_nice_over_1m(i) == 1
        plot(r_SNR(i), b_Hs_interp(i), 'o', 'Color', 'm')
    end
end
hold off;

set(gcf, 'Position', [0, 0, 1820, 980]);
xlim([0.4, 2]);
ylim([0, 5]);
title("SNR-Hs", "FontSize", 15);
xlabel("SNR"); ylabel("Hs [m]");

%% 이미지 분류
close all;

d = dir('./img/*.png');

for i = 1 : length(d)
    if idx_nice_over_1m(i) == 1


        figure(7);
        set(gcf, 'Position', [0, 0, 1820, 500]);
        t = tiledlayout(2,2);
        maintitle = [datestr(r_Date(i), 'yyyy/mm/dd HH:MM'), '  &  Bouy = ', num2str(round(b_Hs_interp(i), 2)), 'm  &  ', 'Radar = ', num2str(round(r_Hs(i), 2)), 'm'];
        title(t, maintitle, "FontSize", 15);


        nexttile;
        hold on;
        plot(r_Date, movmean(r_Hs, 1));
        plot(r_Date, movmean(b_Hs_interp, 1));
        highlight_start = r_Date(i-1); % 시작 시간
        highlight_end = r_Date(i+1); % 끝나는 시간
        y_limits = ylim; % y축 범위 가져오기
        fill([highlight_start highlight_end highlight_end highlight_start], ...
            [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
            [1 1 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', 'm'); % 노란색 형광펜 효과
        hold off;

        if i < 60
            xlim([r_Date(1), r_Date(i+60)]);
        elseif i > length(d) - 60
            xlim([r_Date(i-60), r_Date(length(d))]);
        else
            xlim([r_Date(i-60), r_Date(i+60)]);
        end
        title("Significant Wave Height");
        xlabel("Date [mm dd]"); ylabel("Hs [m]");
        legend('Rader', 'Bouy');


        nexttile(2, [2, 1]);
        hold on;
        plot(b_Hs_interp, r_Hs, '.', 'MarkerSize', 3);
        plot(0:0.5:5, 0:0.5:5, '--', 'Color', 'r', 'LineWidth', 1);
        plot(0:0.5:5, 1:0.5:6, '--', 'Color', 'k', 'LineWidth', 1);
        plot(0:0.5:5, -1:0.5:4, '--', 'Color', 'k', 'LineWidth', 1);
        plot(b_Hs_interp(i), r_Hs(i), 'o', 'Color', 'm', 'LineWidth', 2)
        hold off;

        xlim([0, 5]);
        ylim([0, 5]);
        title(['R^2 = ', num2str(R_squared)]);
        xlabel("Bouy Hs [m]"); ylabel("Rader Hs [m]");


        nexttile;
        hold on;
        plot(r_SNR, b_Hs_interp, '.', 'MarkerSize', 3);
        plot(0:1:10, x(1)+x(2)*(0:1:10), 'r--');
        plot(0:1:10, x(1)+x(2)*(0:1:10)+1, 'k--');
        plot(0:1:10, x(1)+x(2)*(0:1:10)-1, 'k--');
        plot(r_SNR(i), b_Hs_interp(i), 'o', 'Color', 'm', 'LineWidth', 2)
        hold off;

        set(gcf, 'Position', [0, 0, 1820, 980]);
        xlim([0.4, 2]);
        ylim([0, 5]);
        title("SNR-Hs");
        xlabel("SNR"); ylabel("Hs [m]");


        savename = ['./nice/', num2str(i), '_', 'Total_', datestr(r_Date(i), 'yyyymmdd_HHMM'), '.png'];
        saveas(t, savename);


        originalname = [d(i).folder, '/', d(i).name];
        copyname = ['./nice/', num2str(i), '_', d(i).name];
        copyfile(originalname, copyname);

    end
end

%% top 추출
clc;
E = b_Hs_interp - r_Hs;
abs_E = abs(E);
[sort_E, sort_E_idx] = sort(abs_E);
%nice top
disp('nice------------------------------------');
for i = 1 : 100
    if b_Hs_interp(sort_E_idx(i)) > 2
        disp(num2str(sort_E_idx(i)));
        disp(datestr(r_Date(sort_E_idx(i)), 'yyyymmdd_HHMM'));
    end
end
%bad top
disp('bad------------------------------------');
for i = length(sort_E) - 100 : -1 : length(sort_E) - 120
    disp(num2str(sort_E_idx(i)));
    disp(datestr(r_Date(sort_E_idx(i)), 'yyyymmdd_HHMM'));
end