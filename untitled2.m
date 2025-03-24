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
load("snr_y1910_wave.mat");
r_Date = Date;
r_wave_SNR = SNR;
r_wave_Ux = Ux;
r_wave_Uy = Uy;

load("snr_y1911_wave.mat");
r_Date = [r_Date ; Date];
r_wave_SNR = [r_wave_SNR ; SNR];
r_wave_Ux = [r_wave_Ux ; Ux];
r_wave_Uy = [r_wave_Uy ; Uy];

load("snr_y1912_wave.mat");
r_Date = [r_Date ; Date];
r_wave_SNR = [r_wave_SNR ; SNR];
r_wave_Ux = [r_wave_Ux ; Ux];
r_wave_Uy = [r_wave_Uy ; Uy];

r_wave_SNR = sqrt(r_wave_SNR);

% Radar surf
load("snr_y1910_surf.mat");
r_Date = Date;
r_surf_SNR = SNR;
r_surf_Ux = Ux;
r_surf_Uy = Uy;

load("snr_y1911_surf.mat");
r_Date = [r_Date ; Date];
r_surf_SNR = [r_surf_SNR ; SNR];
r_surf_Ux = [r_surf_Ux ; Ux];
r_surf_Uy = [r_surf_Uy ; Uy];

load("snr_y1912_surf.mat");
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
figure(1);
t = tiledlayout(2,1);
a = nexttile;
hold on;
% yyaxis left
plot(b2r_Date, movmean(b2r_Hs, 3), 'Color', [0, 0, 0, 0.8], 'LineStyle', '-');
plot(r2b_Date, movmean(r2b_wave_Hs, 3), 'Color', [1, 0, 0, 0.8], 'LineStyle', '-');
plot(r2b_Date, movmean(r2b_surf_Hs, 3), 'Color', [0, 0, 1, 0.8], 'LineStyle', '-');
ylabel("Hs [m]");
% yyaxis right
% plot(b2r_Date, movmean(b2r_Wind, 3), 'Color', [0.9, 0.2, 0.9, 0.7], 'LineStyle', '-.');
% ylabel("Wind Velocity [m/s]");
hold off;

set(gcf, 'Position', [0, 0, 1820, 980]);
xlim([datetime(2019, 10, 1), datetime(2019, 10, 4)]);
title("Significant Wave Height", "FontSize", 15);
xlabel("Date [mm dd]");
legend('Bouy', 'Wave', 'Surf');
% legend('Bouy', 'Wave', 'Surf', 'Wind Velocity');



% 데이터 커서 모드 활성화
dcm = datacursormode(gcf);
datacursormode on;

% 추가 인수 (예: 추가 정보를 담은 변수)
extraParam1 = b2r_Date;
extraParam2 = b2r_Wind;

% 익명 함수로 추가 인수 전달
set(dcm, 'UpdateFcn', @(obj, event_obj) cursorCallback(obj, event_obj, extraParam1, extraParam2));

b = nexttile;

% 클릭 시 좌표를 출력하는 함수
function txt = cursorCallback(~, event_obj, extraParam1, extraParam2)
    pos = event_obj.Position; % [x, y] 값
    clickedDate = datetime(pos(1), 'ConvertFrom', 'datenum') + calyears(2019) + calmonths(9) + caldays(1); % x값을 datetime으로 변환
    
    idx = extraParam1 == clickedDate;
    windspeed = extraParam2(idx);

    % 'yyyy-mm-dd HH:MM' 포맷으로 변환하여 출력
    txt = {sprintf('Date = %s', datestr(clickedDate, 'yyyy-mm-dd HH:MM')), ...
           sprintf('Hs = %.2f [m]', pos(2)), ...
           sprintf('Wind Speed = %.2f m/s', windspeed)}; % 추가 인수 사용

    % 두 번째 타일에 이미지 업데이트
    nexttile(2);
    imshow(['C:\Users\Hyeonjong Im\Documents\새 폴더\image\Image_', datestr(clickedDate, 'yyyymmdd_HHMM'), '.png']);
end
